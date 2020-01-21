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
top_dir = "../modeling/testsystem/"
#analys_dir = top_dir+'/analys/'
analys_dir = "testsystem/"
analys_dir = '/home/kalabharath/PycharmProjects/trem2_kala/analys/analys/testsystem/'
# How are the trajectories dir names
dir_head = 'run'
out_dirs = glob.glob(top_dir+dir_head+'*/output')

################################
# Extract frames
################################

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Create dir
gsms_A_dir = analys_dir+'GSMs_cl-1/sample_A/'
gsms_B_dir = analys_dir+'GSMs_cl-1/sample_B/'

AT.create_gsms_dir(gsms_A_dir)
AT.create_gsms_dir(gsms_B_dir)

HA1 = AT.get_models_to_extract(analys_dir +'/selected_models_A_cluster3_detailed.csv')
HB1 = AT.get_models_to_extract(analys_dir +'/selected_models_B_cluster3_detailed.csv')
AT.do_extract_models(HA1, 'h1', gsms_A_dir)
AT.do_extract_models(HB1, 'h2', gsms_B_dir)


exit()
