import os
import sys
from joblib import Parallel, delayed


def run_hypertracking(dire):
    os.chdir(dire)
    os.system('o2-hypertracking-workflow -b')

tf_paths = []
dirs = os.listdir()
base_path = os.getcwd()
for dire in dirs:
    if not dire.startswith('tf'):
        continue
    path = base_path + '/' + dire
    tf_paths.append(path)
    os.chdir(path)
    if not os.path.islink('ITSdictionary.bin'):
        os.symlink('/home/fmazzasc/alice/run_sim/ITSdictionary.bin', 'ITSdictionary.bin')
    if not os.path.islink('o2sim_geometry.root'):
        os.symlink('/home/fmazzasc/alice/run_sim/bkg_geometry.root', 'o2sim_geometry.root')
    if not os.path.islink('o2sim_geometry-aligned.root'):
        os.symlink('/home/fmazzasc/alice/run_sim/bkg_geometry-aligned.root', 'o2sim_geometry-aligned.root')
    
    if not os.path.islink('o2sim_grp.root'):
        os.symlink('/home/fmazzasc/alice/run_sim/o2sim_grp.root', 'o2sim_grp.root')

# run_hypertracking(tf_paths[0])
results = Parallel(n_jobs=len(tf_paths))(delayed(run_hypertracking)(dire) for dire in tf_paths)
