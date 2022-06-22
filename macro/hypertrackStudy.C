#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"

#include <gsl/gsl>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "CommonDataFormat/RangeReference.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "StrangenessTracking/HyperTracker.h"

#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = 211;

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcV0alpha(const V0 &v0);
double recomputeV0Pt(const V0 &v0);
std::vector<std::array<int, 2>> matchV0stoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec);
std::array<int, 2> matchV0DautoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, o2::dataformats::V0::GIndex dauID);
std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx);

void hypertrackStudy()
{
    // Output Histograms
    TH1D *hChi2Sgn = new TH1D("Chi2 Signal", "; #chi^{2}; Counts", 102, -2, 100);
    TH1D *hChi2Bkg = new TH1D("Chi2 Fake assoc", "; #chi^{2}; Counts", 102, -2, 100);

    TH1D *hSigBkg = new TH1D("Hypertracker eff", "; ; Efficiency", 2, 0, 2);
    TH2D *hMChisto = new TH2D("histo mc", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH2D *hV0histo = new TH2D("histo V0", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH2D *hHyperhisto = new TH2D("histo hyperV0", ";#it{p}_{T} (GeV/#it{c}); Radius^2 (cm) ; Counts", 20, 1, 10, 30, 1, 900);
    TH1D *hResV0histo = new TH1D("pT resolution before hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResHyperhisto = new TH1D("pT resolution after hypertracking", ";(#it{p}_{T}^{gen} - #it{p}_{T}^{rec})/#it{p}_{T}^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResV0histoR2 = new TH1D("R2 resolution before hypertracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hResHyperhistoR2 = new TH1D("R2 resolution after hypertracking", ";(R2^{gen} - R2^{rec})/R2^{gen}; Counts", 20, -0.2, 0.2);
    TH1D *hV0Counter = new TH1D("V0 counter", ";V0 counter; Counts", 1, 0.5, 1.5);

    TH1D *hRecHypCounter = new TH1D("Rec V0 hyp counter", "; ; Counts", 1, 0.5, 1.5);
    TH1D *hHypertrackerStats = new TH1D("hypertracker_stats", "; ; Counts", 3, 0.5, 3.5);
    TH1D *hHyperCounter = new TH1D("Hypertrack counter", ";Hypertrack counter; Counts", 1, 0.5, 1.5);
    TH1D *hFakeAssocCounter = new TH1D("Fake assoc counter", ";Fake assoc counter; Counts", 1, 0.5, 1.5);

    std::string path = "/data/fmazzasc/its_data/sim/hyp/";
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 4) == "tf16")
        {
            dirs.push_back(path + file);
            auto innerdir = (TSystemDirectory *)fileObj;
            auto innerfiles = innerdir->GetListOfFiles();
            for (auto innerfileObj : *innerfiles)
            {
                TString innerfile = ((TSystemFile *)innerfileObj)->GetName();
                if (innerfile.EndsWith("Kine.root") && innerfile.Contains("sgn"))
                {
                    kine_files.push_back(innerfile);
                }
            }
        }
    }
    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        LOG(info) << "kine file " << kine_file;
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        auto fHyperTracks = TFile::Open((dir + "/o2_hypertrack.root").data());
        auto fSecondaries = TFile::Open((dir + "/o2_secondary_vertex.root").data());
        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((dir + "/o2match_tof_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((dir + "/o2match_tof_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());

        // Geometry
        o2::base::GeometryManager::loadGeometry(dir + "/o2sim_geometry.root");

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeHypertracks = (TTree *)fHyperTracks->Get("o2sim");
        auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");
        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;
        // Hypertracks
        std::vector<V0> *hyperV0vec = nullptr;
        std::vector<int> *hyperITSvec = nullptr;
        std::vector<o2::track::TrackParametrizationWithError<float>> *hypertrackVec = nullptr;
        std::vector<float> *hyperChi2vec = nullptr;
        std::vector<o2::strangeness_tracking::ClusAttachments> *nAttachments = nullptr;

        // Secondary Vertices
        std::vector<V0> *v0vec = nullptr;
        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        // Setting branches
        treeHypertracks->SetBranchAddress("V0s", &hyperV0vec);
        treeHypertracks->SetBranchAddress("ITSTrackRefs", &hyperITSvec);
        treeHypertracks->SetBranchAddress("Hypertracks", &hypertrackVec);
        treeHypertracks->SetBranchAddress("ITSV0Chi2", &hyperChi2vec);
        treeHypertracks->SetBranchAddress("ClusUpdates", &nAttachments);

        treeSecondaries->SetBranchAddress("V0s", &v0vec);
        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);

        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // define detector map
        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}};

        // load geometry
        auto gman = o2::its::GeometryTGeo::Instance();
        gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

        // fill MC matrix
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            mcTracksMatrix[n].resize(MCtracks->size());
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
                if (abs(MCtracks->at(mcI).GetPdgCode()) == motherPDG)
                {
                    auto &mcTrack = mcTracksMatrix[n][mcI];
                    hMChisto->Fill(mcTrack.GetPt(), calcDecLength(MCtracks, mcTrack, firstDaughterPDG));
                }
            }
        }

        // Starting matching  loop
        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) ||
                !treeHypertracks->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame))
                continue;

            std::vector<std::array<int, 2>> V0sMCref = matchV0stoMC(mcTracksMatrix, map, v0vec);

            std::vector<bool> matchedHypertrackVec;
            matchedHypertrackVec.resize(hyperV0vec->size());
            for (auto &&isMatched : matchedHypertrackVec)
                isMatched = false;

            for (unsigned int iV0vec = 0; iV0vec < v0vec->size(); iV0vec++)
            {
                hV0Counter->Fill(1);
                auto &v0MCref = V0sMCref[iV0vec];
                auto &v0 = (*v0vec)[iV0vec];
                if (v0MCref[0] == -1 || v0MCref[1] == -1)
                    continue;
                hRecHypCounter->Fill(1);

                auto &mcTrack = mcTracksMatrix[v0MCref[0]][v0MCref[1]];
                hV0histo->Fill(recomputeV0Pt(v0), v0.calcR2());

                // Matching ITS tracks to MC tracks and V0
                std::array<int, 2> ITSref = {-1, 1};
                o2::its::TrackITS ITStrack;
                std::array<std::array<int, 2>, 7> clsRef;

                int iTrack = -1;
                bool isMatched = false;

                for (unsigned int iITStrack = 0; iITStrack < ITStracks->size(); iITStrack++)
                {
                    auto &labITS = (*labITSvec)[iITStrack];
                    auto &trackIdx = (*ITSTrackClusIdx)[iITStrack];

                    ITSref = matchITStracktoMC(mcTracksMatrix, labITS);

                    if (ITSref[0] == V0sMCref[iV0vec][0] && ITSref[1] == V0sMCref[iV0vec][1])
                    {
                        ITStrack = (*ITStracks)[iITStrack];
                        auto firstClus = ITStrack.getFirstClusterEntry();
                        auto ncl = ITStrack.getNumberOfClusters();
                        for (int icl = 0; icl < ncl; icl++)
                        {
                            auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                            auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                            auto layer = gman->getLayer(clus.getSensorID());
                            clsRef[layer] = matchCompLabelToMC(mcTracksMatrix, labCls);
                        }
                        hHypertrackerStats->Fill(1);
                        isMatched = true;
                        iTrack = iITStrack;
                        break;
                    }
                }

                if (!isMatched)
                    continue;

                // Matching hypertracks to MC tracks, V0s and ITS tracks
                bool isHypertracked = false;
                for (unsigned int iHyperVec = 0; iHyperVec < hyperV0vec->size(); iHyperVec++)
                {
                    auto &hyperV0 = hyperV0vec->at(iHyperVec);
                    auto &hyperTrack = hypertrackVec->at(iHyperVec);
                    auto &hyperChi2 = hyperChi2vec->at(iHyperVec);
                    auto &hyperITSref = hyperITSvec->at(iHyperVec);
                    if (hyperV0.getProngID(0) == v0.getProngID(0) && hyperV0.getProngID(1) == v0.getProngID(1) && hyperITSref == iTrack)
                    {
                        LOG(info) << "++++++++++++++++++++++++";
                        LOG(info) << "Hypertrack found!: ITS track ref: " << hyperITSref;
                        isHypertracked = true;
                        matchedHypertrackVec[iHyperVec] = true;
                        hHyperCounter->Fill(1);
                        hHypertrackerStats->Fill(2);
                        hHyperhisto->Fill(hyperTrack.getPt(), hyperV0.calcR2());
                        hResHyperhisto->Fill((hyperTrack.getPt() - mcTrack.GetPt()) / mcTrack.GetPt());
                        hResV0histo->Fill((recomputeV0Pt(v0) - mcTrack.GetPt()) / mcTrack.GetPt());
                        hResV0histoR2->Fill((v0.calcR2() - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                        hResHyperhistoR2->Fill((hyperV0.calcR2() - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                        hChi2Sgn->Fill(hyperChi2);
                        break;
                    }
                }

                if (!isHypertracked)
                {
                    LOG(info) << "------------------";
                    LOG(info) << "No hypertrack found, but ITS track found";
                    LOG(info) << "processing frame " << frame << ", dir: " << dir;
                    LOG(info) << "V0 pos: " << v0.getProngID(0) << " V0 neg: " << v0.getProngID(1) << " V0pt: " << v0.getPt() << " ITSpt: " << ITStrack.getPt();
                    LOG(info) << "V0 Eta: " << v0.getEta() << " V0 phi" << v0.getPhi() << " ITS eta: " << ITStrack.getEta() << " ITS phi: " << ITStrack.getPhi();

                    LOG(info) << "Number of hits: " << ITStrack.getNClusters();
                    LOG(info) << "ITS Track ref: " << ITSref[0] << " " << ITSref[1] << " , PDG: " << mcTracksMatrix[ITSref[0]][ITSref[1]].GetPdgCode();
                    LOG(info) << "+++++++";
                    hHypertrackerStats->Fill(3);

                    for (unsigned int i{0}; i < 7; i++)
                    {
                        if (ITStrack.hasHitOnLayer(i))
                        {
                            LOG(info) << "ITS track has hit on layer " << i << ", is fake: " << ITStrack.isFakeOnLayer(i);
                            if (clsRef[i][0] != -1 && clsRef[i][1] != -1)
                                LOG(info) << "Cluster ref: " << clsRef[i][0] << " " << clsRef[i][1] << " , PDG: " << mcTracksMatrix[clsRef[i][0]][clsRef[i][1]].GetPdgCode();
                            if (clsRef[i][0] != ITSref[0])
                            {
                                LOG(info) << "EvID mismatch: " << clsRef[i][0] << " " << ITSref[0];
                                continue;
                            }
                            if (clsRef[i][1] != ITSref[1])
                            {

                                auto motherID = mcTracksMatrix[clsRef[i][0]][clsRef[i][1]].getMotherTrackId();
                                if (motherID != -1)
                                    LOG(info) << "Mother cluster ref: " << clsRef[i][0] << " " << motherID << " , PDG: " << mcTracksMatrix[clsRef[i][0]][motherID].GetPdgCode();
                            }
                        }
                    }
                }
            }

            LOG(info) << " Start studying fake associations";
            for (unsigned int iHyperVec = 0; iHyperVec < hyperV0vec->size(); iHyperVec++)
            {

                auto &hyperChi2 = hyperChi2vec->at(iHyperVec);
                auto &clusAttachments = nAttachments->at(iHyperVec);
                auto &ITStrack = ITStracks->at(hyperITSvec->at(iHyperVec));
                auto &ITStrackLab = labITSvec->at(hyperITSvec->at(iHyperVec));

                auto clusAttArr = clusAttachments.arr;
                auto isMatched = matchedHypertrackVec[iHyperVec];
                auto &hyperV0 = hyperV0vec->at(iHyperVec);

                if (!isMatched)
                {
                    LOG(info) << "**************";
                    LOG(info) << "V0 pos: " << hyperV0.getProngID(0) << " V0 neg: " << hyperV0.getProngID(1) << " V0pt: " << hyperV0.getPt() << " ITSpt: " << ITStrack.getPt();
                    auto v0PosRef = matchV0DautoMC(mcTracksMatrix, map, hyperV0.getProngID(0));
                    auto v0NegRef = matchV0DautoMC(mcTracksMatrix, map, hyperV0.getProngID(1));
                    auto ITStrackRef = matchCompLabelToMC(mcTracksMatrix, ITStrackLab);
                    LOG(info) << "V0 pos lab: " << v0PosRef[0] << " " << v0PosRef[1] << ", V0 neg lab: " << v0NegRef[0] << " " << v0NegRef[1] << ", ITS ref: " << ITStrackRef[0] << " " << ITStrackRef[1];  


                    hChi2Bkg->Fill(hyperChi2);
                    hFakeAssocCounter->Fill(1);
                    auto firstClus = ITStrack.getFirstClusterEntry();
                    auto ncl = ITStrack.getNumberOfClusters();
                    for (int icl = 0; icl < ncl; icl++)
                    {
                        auto &labCls = (clusLabArr->getLabels(ITSTrackClusIdx->at(firstClus + icl)))[0];
                        auto &clus = (*ITSclus)[(*ITSTrackClusIdx)[firstClus + icl]];
                        auto layer = gman->getLayer(clus.getSensorID());
                        std::array<int, 2> clsRef = matchCompLabelToMC(mcTracksMatrix, labCls);
                        if (clsRef[0] > -1 && clsRef[1] > -1)
                            LOG(info) << "Layer: " << layer << "PDG: " << mcTracksMatrix[clsRef[0]][clsRef[1]].GetPdgCode() << ", Attached to: " << clusAttArr[layer];
                        else
                            LOG(info) << "Layer: " << layer << ", No valid cluster ref";
                    }
                }
            }
        }
        auto outFile = TFile("hypertrack_study.root", "recreate");
        hChi2Sgn->Write();
        hChi2Bkg->Write();

        hV0histo->Write();
        hResV0histo->Write();
        hResV0histoR2->Write();
        hHyperhisto->Write();
        hResHyperhisto->Write();
        hResHyperhistoR2->Write();
        hMChisto->Write();

        hV0Counter->Write();
        hHypertrackerStats->Write();
        hHyperCounter->Write();
        hFakeAssocCounter->Write();
        hRecHypCounter->Write();
        outFile.Close();
    }
}
std::vector<std::array<int, 2>> matchV0stoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec)
{
    std::vector<std::array<int, 2>> outArray;
    outArray.resize(v0vec->size());
    int count_V0 = 0;
    for (unsigned int iV0vec = 0; iV0vec < v0vec->size(); iV0vec++)
    {
        std::vector<int> motherIDvec;
        std::vector<int> daughterIDvec;
        std::vector<int> evIDvec;

        outArray[iV0vec] = {-1, -1};
        auto &v0 = (*v0vec)[iV0vec];

        for (unsigned int iV0 = 0; iV0 < 2; iV0++)
        {
            if (map[v0.getProngID(iV0).getSourceName()])
            {
                auto labTrackType = map[v0.getProngID(iV0).getSourceName()];
                auto lab = labTrackType->at(v0.getProngID(iV0).getIndex());

                int trackID, evID, srcID;
                bool fake;
                lab.get(trackID, evID, srcID, fake);
                if (lab.isValid())
                {
                    auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
                    motherIDvec.push_back(motherID);
                    daughterIDvec.push_back(trackID);
                    evIDvec.push_back(evID);
                }
            }
        }

        if (motherIDvec.size() < 2)
            continue;
        if (motherIDvec[0] != motherIDvec[1] || evIDvec[0] != evIDvec[1])
            continue;

        if (motherIDvec[0] <= 0 || motherIDvec[0] > 10000)
            continue;

        int pdg0 = mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPdgCode();
        int pdg1 = mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPdgCode();

        if (!(std::abs(pdg0) == firstDaughterPDG && std::abs(pdg1) == secondDaughterPDG) && !(std::abs(pdg0) == secondDaughterPDG && std::abs(pdg1) == firstDaughterPDG))
            continue;

        // std::cout << "Mother PDG: " << mcTracksMatrix[evIDvec[0]][motherIDvec[0]].GetPt() << std::endl;
        outArray[iV0vec] = {evIDvec[0], motherIDvec[0]};
        count_V0++;
    }
    std::cout << "Number of V0s: " << count_V0 << std::endl;
    return outArray;
}

std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel)

{
    std::array<int, 2> outArray = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    ITSlabel.get(trackID, evID, srcID, fake);
    if (!ITSlabel.isNoise() && ITSlabel.isValid() && std::abs(mcTracksMatrix[evID][trackID].GetPdgCode()) == motherPDG)
    {
        outArray = {evID, trackID};
    }

    return outArray;
}

std::array<int, 2> matchCompLabelToMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel compLabel)
{
    std::array<int, 2> compRef = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    compLabel.get(trackID, evID, srcID, fake);
    if (compLabel.isValid())
    {
        compRef = {evID, trackID};
    }
    return compRef;
}

std::array<int, 2> matchV0DautoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, o2::dataformats::V0::GIndex dauID)
{
    std::array<int, 2> outArray{-1, -1};

    if (map[dauID.getSourceName()])
    {
        auto labTrackType = map[dauID.getSourceName()];
        auto lab = labTrackType->at(dauID.getIndex());

        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        if (lab.isValid())
        {
            auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
            outArray = {evID, motherID};
        }
    }

    return outArray;
}

std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx)
{

    std::vector<ITSCluster> outVec;
    auto firstClus = ITStrack.getFirstClusterEntry();
    auto ncl = ITStrack.getNumberOfClusters();
    for (int icl = 0; icl < ncl; icl++)
    {
        outVec.push_back(ITSClustersArray[(*ITSTrackClusIdx)[firstClus + icl]]);
    }
    return outVec;
}

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (std::abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto decLength = (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) *
                                 (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) +
                             (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) *
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY());
            return decLength;
        }
    }
    return -1;
}

double calcV0alpha(const V0 &v0)
{
    std::array<float, 3> fV0mom, fPmom, fNmom = {0, 0, 0};
    v0.getProng(0).getPxPyPzGlo(fPmom);
    v0.getProng(1).getPxPyPzGlo(fNmom);
    v0.getPxPyPzGlo(fV0mom);

    TVector3 momNeg(fNmom[0], fNmom[1], fNmom[2]);
    TVector3 momPos(fPmom[0], fPmom[1], fPmom[2]);
    TVector3 momTot(fV0mom[0], fV0mom[1], fV0mom[2]);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

double recomputeV0Pt(const V0 &v0)
{
    double alpha = calcV0alpha(v0);
    std::array<float, 3> fPmom, fNmom = {0, 0, 0};
    v0.getProng(0).getPxPyPzGlo(fPmom);
    v0.getProng(1).getPxPyPzGlo(fNmom);

    if (alpha > 0)
        return sqrt((2 * fPmom[0] + fNmom[0]) * (2 * fPmom[0] + fNmom[0]) + (2 * fPmom[1] + fNmom[1]) * (2 * fPmom[1] + fNmom[1]));

    return sqrt((fPmom[0] + 2 * fNmom[0]) * (fPmom[0] + 2 * fNmom[0]) + (fPmom[1] + 2 * fNmom[1]) * (fPmom[1] + 2 * fNmom[1]));
}