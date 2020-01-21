from multiprocessing import Process
import time
import resource
import numpy as np
import scipy
import random
import math
import sys, os, glob
import time

from sklearn.cluster import KMeans
import scipy as sp
from scipy import spatial

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi
import IMP.pmi.analysis
import os
import sys
from math import sqrt


rmf_ref = "cluster_center_model.rmf3"

m = IMP.Model()
h_ref = IMP.pmi.analysis.get_hiers_from_rmf(m,0,rmf_ref)[0]
print(h_ref)
s0 = IMP.atom.Selection(h_ref, resolution=1)
print(s0)

pdb_name = "cluster_center_model.pdb"

o = IMP.pmi.output.Output()
o.init_pdb(pdb_name,h_ref)
o.write_pdb(pdb_name)
del o

