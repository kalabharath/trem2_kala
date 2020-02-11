# Analysis using PMI analysis
# (1) initalize the analysis class - AnalysisTrajectories (run_analysis_trajectories.py)
# (2) Add restraints you want to be analyzed
# (3) Read stat files to obtain the scores, nuisance parameters, and info about the rmf files
# (4) Obtain the statstics of the XL restraints Psi parameter
# (5) Do HBDSCAN clustering for selected scores and/or nuisance parameters
# (6) Get info about XL satisfaction
# (7) Re-run clustering without having to re-read the stat files
# (8) Extract models from the rmf  (run_extract_models.py)
# (9) Test for convergence and do structural clustering (run_clustering.py)


import numpy as np
import pandas as pd
import math
import glob
import sys
import os
sys.path.append('/home/kalabharath/Prion/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 4
#top_dir =  sys.argv[1]
top_dir = "../modeling/testsystem/"
#analys_dir = top_dir+'/analys/'
analys_dir = "analys/testsystem/"
#analys_dir = "analys/splitXLs/"
# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run'
#out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')
out_dirs = glob.glob(top_dir+ dir_head+ '*/output/')
print(out_dirs)

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'NHSF':30.0, 'BS3':34.0} #, 'intra':30.0}
#XLs_cutoffs = {'NHSF_KST':20.0, 'NHSF_KKY':24.0, 'NHSF_Nter':18.0, 'NHSF_KH':22.0, 'BS3':28.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc,
                          nskip=20)
                          #th = 200)
AT.ambiguous_XLs_restraint= True
AT.Multiple_XLs_restraints= True

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
                             ambiguous_XLs_restraint = True,
                             Multiple_XLs_restraints = True,
                             get_nuisances = False)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
# AT.set_analyze_score_only_restraint('barrier')
# AT.set_analyze_EM_restraint()
AT.set_select_by_Total_score(50.0)

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.read_models_info(XLs_cutoffs)
#AT.get_Psi_stats()

AT.hdbscan_clustering(['XLs_NHSF', 'XLs_BS3']) #, 'XLs_intra'])
# AT.hdbscan_clustering(['EV_sum', 'XLs_NHSF_KST', 'XLs_NHSF_KKY', 'XLs_NHSF_Nter', 'XLs_NHSF_KH', 'XLs_BS3'])
# AT.summarize_XLs_info(Multiple_XLs_restraints = True, ambiguous_XLs_restraint = True)

AT.summarize_sampling_info()
exit()
