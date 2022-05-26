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
const int secondDaughterPDG = -211;

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcV0alpha(const V0 &v0);
double recomputeV0Pt(const V0 &v0);
std::vector<std::array<int, 2>> matchV0stoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, std::map<std::string, std::vector<o2::MCCompLabel> *> &map, std::vector<V0> *v0vec);
std::array<int, 2> matchITStracktoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::array<int, 2> matchITSclustoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel);
std::vector<ITSCluster> getTrackClusters(const o2::its::TrackITS &ITStrack, const std::vector<ITSCluster> &ITSClustersArray, std::vector<int> *ITSTrackClusIdx);
std::string get_str_between_two_str(const std::string &s,
                                    const std::string &start_delim,
                                    const std::string &stop_delim)
{
    unsigned first_delim_pos = s.find(start_delim);
    unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
    unsigned last_delim_pos = s.find(stop_delim);

    return s.substr(end_pos_of_first_delim,
                    last_delim_pos - end_pos_of_first_delim);
}

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

    TH1D *hRecHypCounter = new TH1D("Rec V0 hyp counter", ";Rec V0 hyp counter; Counts", 1, 0.5, 1.5);
    TH1D *hV0wTrackCounter = new TH1D("V0 with track counter", ";V0 with track counter; Counts", 1, 0.5, 1.5);
    TH1D *hHyperCounter = new TH1D("Hypertrack counter", ";Hypertrack counter; Counts", 1, 0.5, 1.5);
    TH1D *hFakeAssocCounter = new TH1D("Fake assoc counter", ";Fake assoc counter; Counts", 1, 0.5, 1.5);

    std::string path = "/home/fmazzasc/alice/run_sim/";
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            dirs.push_back(file);
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
        o2::base::GeometryManager::loadGeometry(dir + "/o2sim_geometry-aligned.root");

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

        // Topology dictionary
        o2::itsmft::TopologyDictionary mdict;
        mdict.readFromFile(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::ITS));

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;
        // Hypertracks
        std::vector<V0> *hyperV0vec = nullptr;
        std::vector<int> *hyperITSvec = nullptr;
        std::vector<o2::track::TrackParametrizationWithError<float>> *hypertrackVec = nullptr;
        std::vector<float> *hyperChi2vec = nullptr;
        std::vector<o2::strangeness_tracking::He3Attachments> *nUpdates = nullptr;

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
        treeHypertracks->SetBranchAddress("He3Updates", &nUpdates);

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
                if (MCtracks->at(mcI).GetPdgCode() == motherPDG)
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
                            clsRef[layer] = matchITSclustoMC(mcTracksMatrix, labCls);
                        }
                        hV0wTrackCounter->Fill(1);
                        iTrack = iITStrack;
                        break;
                    }
                }

                bool isMatched = true;
                if (!(ITSref[0] == V0sMCref[iV0vec][0]) || !(ITSref[1] == V0sMCref[iV0vec][1]))
                    isMatched = false;

                // Matching hypertracks to MC tracks, V0s and ITS tracks
                int iHypertrack = -1;
                for (unsigned int iHyperVec = 0; iHyperVec < hyperV0vec->size(); iHyperVec++)
                {
                    auto &hyperV0 = (*hyperV0vec)[iHyperVec];
                    if (hyperV0.getProngID(0) == v0.getProngID(0) || hyperV0.getProngID(1) == v0.getProngID(1))
                    {
                        iHypertrack = iHyperVec;
                        break;
                    }
                }

                if (iHypertrack == -1 && isMatched)
                {
                    LOG(info) << "------------------";
                    LOG(info) << "No hypertrack found, but ITS track found";
                    LOG(info) << "processing frame " << frame << ", dir: " << dir;
                    LOG(info) << "V0 pos: " << v0.getProngID(0) << " V0 neg: " << v0.getProngID(1) << " ITS: " << iTrack;
                    LOG(info) << "Number of hits: " << ITStrack.getNClusters();
                    LOG(info) << "ITS Track ref: " << ITSref[0] << " " << ITSref[1] << " , PDG: " << mcTracksMatrix[ITSref[0]][ITSref[1]].GetPdgCode();
                    LOG(info) << "+++++++";
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

                if (iHypertrack == -1)
                {
                    continue;
                }

                auto &hyperTrack = (*hypertrackVec)[iHypertrack];
                auto &hyperV0 = (*hyperV0vec)[iHypertrack];
                auto &hyperChi2 = (*hyperChi2vec)[iHypertrack];

                std::array<unsigned int, 7> &isUpdated = (*nUpdates)[iHypertrack].arr;
                if (isMatched)
                {
                    hHyperCounter->Fill(1);
                    hHyperhisto->Fill(hyperTrack.getPt(), hyperV0.calcR2());
                    hResHyperhisto->Fill((hyperTrack.getPt() - mcTrack.GetPt()) / mcTrack.GetPt());
                    hResV0histo->Fill((recomputeV0Pt(v0) - mcTrack.GetPt()) / mcTrack.GetPt());
                    hResV0histoR2->Fill((v0.calcR2() - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                    hResHyperhistoR2->Fill((hyperV0.calcR2() - calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG)) / calcDecLength(&mcTracksMatrix[v0MCref[0]], mcTrack, firstDaughterPDG));
                    hChi2Sgn->Fill(hyperChi2);
                }
                else
                {
                    hChi2Bkg->Fill(hyperChi2);
                    hFakeAssocCounter->Fill(1);
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
    hV0wTrackCounter->Write();
    hHyperCounter->Write();
    hFakeAssocCounter->Write();
    hRecHypCounter->Write();
    outFile.Close();
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
                if (!lab.isNoise() && lab.isValid() && lab.isCorrect() && srcID)
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

        if (pdg0 != firstDaughterPDG && pdg0 != secondDaughterPDG)
            continue;
        if (pdg1 != firstDaughterPDG && pdg1 != secondDaughterPDG)
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
    if (!ITSlabel.isNoise() && ITSlabel.isValid() && srcID && mcTracksMatrix[evID][trackID].GetPdgCode() == motherPDG)
    {
        outArray = {evID, trackID};
    }

    return outArray;
}

std::array<int, 2> matchITSclustoMC(const std::vector<std::vector<o2::MCTrack>> &mcTracksMatrix, o2::MCCompLabel ITSlabel)
{
    std::array<int, 2> clusRef = {-1, -1};
    int trackID, evID, srcID;
    bool fake;
    ITSlabel.get(trackID, evID, srcID, fake);
    if (ITSlabel.isValid() && srcID)
    {
        clusRef = {evID, trackID};
    }
    return clusRef;
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
        if (dauTrack.GetPdgCode() == dauPDG)
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