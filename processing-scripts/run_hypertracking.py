from distutils.log import debug
import os
import sys
from joblib import Parallel, delayed


debug_level = "debug"
base_path = "/data/fmazzasc/its_data/sim/hyp"

def run_hypertracking(dire, debug_level='info'):
    os.chdir(dire)
    os.system(f'o2-hypertracking-workflow -b --severity {debug_level}')

tf_paths = []
dirs = os.listdir(base_path)
for dire in dirs:
    if not dire.startswith('tf16'):
        continue
    path = base_path + '/' + dire
    tf_paths.append(path)
    os.chdir(path)
    files_list = os.listdir(path)
    if os.path.islink('o2sim_geometry.root'):
        os.unlink('o2sim_geometry.root')
    if os.path.islink('o2sim_grp.root'):
        os.unlink('o2sim_grp.root')
    geom_file = [f for f in files_list if (f.endswith('_geometry.root') and f.startswith('sgn'))][0]
    grp_file = [f for f in files_list if (f.endswith('_grp.root') and f.startswith('sgn'))][0]
    os.symlink(geom_file, 'o2sim_geometry.root')
    os.symlink(grp_file, 'o2sim_grp.root')

    # print(path)
    # print(geom_file, ' ', grp_file)

# run_hypertracking(tf_paths[0])
results = Parallel(n_jobs=len(tf_paths))(delayed(run_hypertracking)(dire, debug_level) for dire in tf_paths)
