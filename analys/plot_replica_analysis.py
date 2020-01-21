import glob
import numpy as np
import argparse

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import matplotlib.colors as colors
mpl.rcParams.update({'font.size': 8})

################################################    
files = glob.glob('/home/rakesh/WORK/Projects/4-SPOTSComplex/IMPModeling/testRun/results/spots_1/output/replica_temp_acc.*.dat')
th = 500

# Read replicas information
S = dict()
temp = set()
for f in files:
    temp_vals = []
    replica = f.split('.')[1]
    print(replica)
    for line in open(f):
        vals=line.split()
        if vals[0]=='>':
            temp_vals.append([vals[1],vals[2],vals[3],vals[4],vals[5]])
    S[replica] = np.array(temp_vals,dtype='float')
    if len(S[replica])>th:
        ttemp = set(S[replica][:,0])
        temp = temp.union(ttemp)
        
temp_ord =  np.sort([np.float(t) for t in temp])
for i in range(7):
    print(i,i+1,temp_ord[i],temp_ord[i+1],temp_ord[i]/temp_ord[i+1])
    

# Accumulate scores per-temperature
Sc=dict()
for t in temp_ord:
    Sc[t] = list()
    print(t)
    
for i in range(len(S)):
    for t in temp_ord:
        if np.shape(S[str(i)])[0]>0:
            M = S[str(i)][th:,:]
            Sc[t] += list(M[M[:,0]==t,4])
        
# Compute heat capacity
Cv = np.zeros((len(temp_ord),2))
for i,t in enumerate(temp_ord):
    Cv[i,0] = t
    Cv[i,1] = (np.mean((np.array(Sc[t])*np.array(Sc[t])))-np.mean(Sc[t])**2)/t


# Read movers information
D = np.zeros((5000,4))
mF = glob.glob('output_wEM10.0_8/score_movacc.*.dat')
for f in mF:
    for line in open(f):
        vals=line.split()
        if vals[0]=='>':
            D[int(vals[1]),0] = int(vals[1])
            D[int(vals[1]),1] = float(vals[2])
            D[int(vals[1]),2] = float(vals[3])
            D[int(vals[1]),3] = float(vals[4])
            
#########################   
# Plot stuff
#########################   
jet = pl.get_cmap('jet', len(S))
nr = len(S)
fig = pl.figure(figsize=(10,6))
gs = gridspec.GridSpec(3, 2,
                       width_ratios = [0.5,0.5],
                       height_ratios = [0.33,0.33,0.33])


ax = pl.subplot(gs[0])
ax.plot(temp_ord,marker='o')
ax.set_xlabel('Replica', fontsize=10)
ax.set_ylabel('Temperature', fontsize=10)
ax.set_title("Replica's temperatures",fontsize=11)

# Get path for each replica
ax = pl.subplot(gs[1])
for i in range(len(S)):
    if np.shape(S[str(i)])[0]>0:
        ax.plot(S[str(i)][:,0],c=pl.cm.bwr(i/float(nr)))
ax.set_xlim([4800,5000])
print(temp_ord)
ax.set_ylim([0.9*temp_ord[0],1.1*temp_ord[-1]])
ax.set_xlabel('Step', fontsize=10)
ax.set_ylabel('Temperature', fontsize=10)
ax.set_title("Core's temperatures",fontsize=10)

# Histograms of scores
max_sc = 0
min_sc = 9999999
ax = pl.subplot(gs[2])
i=0
for t in Sc.keys():
    ax.hist(Sc[t],histtype='step', stacked=True, fill=False, color=pl.cm.nipy_spectral(i/float(nr)))
    if np.max(Sc[t]) > max_sc:
        max_sc = np.max(Sc[t])
    if np.min(Sc[t]) < min_sc:
        min_sc = np.min(Sc[t])
    i += 1

ax.set_xlim([min_sc,max_sc])
ax.set_ylim([0,5000])
ax.set_xlabel('Score', fontsize=10)
ax.set_ylabel('Frequency', fontsize=10)
ax.set_title("Scores distributions",fontsize=10)


# Plot heat capacity
ax = pl.subplot(gs[3])
ax.plot(Cv[:,0],Cv[:,1],marker='o')
ax.set_xlabel('Temperature', fontsize=10)
ax.set_ylabel('Cv', fontsize=10)

# Plot movers acc. probability
ax = pl.subplot(gs[4])
ax.hist(D[500:,2],histtype='step', stacked=True, fill=False,color='red')
ax.hist(D[500:,3],histtype='step', stacked=True, fill=False,color='blue')
ax.set_xlabel('Acc. probability', fontsize=10)
ax.set_ylabel('Freq.', fontsize=10)

ax = pl.subplot(gs[5])
ax.plot(D[1000:,2], c='red')
ax.plot(D[1000:,3], c='blue') 
ax.set_xlabel('Step', fontsize=10)
ax.set_ylabel('Acc. probability', fontsize=10)

fig.tight_layout()
fig.savefig('replica_analys_w10_r8.pdf') 
pl.close()
    



