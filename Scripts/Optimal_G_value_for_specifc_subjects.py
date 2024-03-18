# Importage 

import warnings
warnings.filterwarnings('ignore')

import os, sys, scipy.io, numpy as np, pandas as pd, seaborn as sns
import time
from matplotlib import pyplot as plt

from scipy.io import loadmat, savemat

import nibabel as nib
from nilearn.plotting import plot_surf, plot_surf_stat_map, plot_roi, plot_anat, plot_surf_roi
from nilearn.image import index_img

import scipy
import networkx as nx

tvb_lib_dir = '/external/rprshnas01/kcni/hharita/Code/shreyVB/tvb-library'
tvb_dat_dir = '/external/rprshnas01/kcni/hharita/Code/shreyVB/tvb-data'
sys.path += [tvb_lib_dir,tvb_dat_dir]
from tvb.simulator.lab import *
from tvb.datatypes.time_series import TimeSeriesRegion
# import tvb.analyzers.correlation_coefficient as corr_coeff


sns.set()

### ========================================================================================================================================


start_time = time.time()

sub_id = sys.argv[1]

sub_id = int(sub_id)

print("Imports, path setup complete!")

### ========================================================================================================================================

# Initialise the Deco 2014 Model.

rww_2014 = models.ReducedWongWangExcInh()

# -------------------------------------------------------------------------------------------------------------------------------------

class MyRWW(models.ReducedWongWangExcInh):
    def dfun(self, state, coupling, local_coupling=0.0):
        # save the x and H value as attribute on object
        S = state
        c_0 = coupling
        lc_0 = local_coupling * S
        self.x = self.w * self.J_N * S + self.I_o + self.J_N * c_0 + self.J_N * lc_0
        self.H = (self.a*self.x - self.b) / (1 - np.exp(-self.d*(self.a*self.x - self.b)))
        # call the default implementation
        return super(MyRWW, self).dfun(state, coupling, local_coupling=local_coupling)
    
# -------------------------------------------------------------------------------------------------------------------------------------
    
def __(num): return np.array([num])

# -------------------------------------------------------------------------------------------------------------------------------------

def run_rww_sim_pcc_test(con, G, regime, dt, simlen, noise_term):
    
    # put regime dict vals into np arrays
    regime = {k: __(v) for k,v in regime.items()}
    
    # Initialise Simulator.
    sim = simulator.Simulator(
        model=models.ReducedWongWangExcInh(**regime),
        connectivity=con,
        coupling=coupling.Scaling(a=__(G)),
#         mynoise = noise.Additive(nsig=__((D**2)/2),noise_seed=np.random.randint(0,1000)),
        integrator=integrators.HeunStochastic(dt=dt, noise=noise_term),
        monitors=(monitors.TemporalAverage(period=1.),
                  monitors.Bold(period=200.0),)
    )
    sim.initial_conditions = (0.3)*np.ones((1, 2, nregions, 1)) #np.random.rand(nregions,1) #(0.001)*np.ones((1, 1, nregions, 1))
    sim.configure()
    
    
    # Launch simulation
        
    res = sim.run(simulation_length=simlen)
    (Sts, Ss),(Bts,Bs) = res

    # Bts and Sts contain 2 sets of data -- one excitatory and one inhibitory? ask JG. 
    # Regardless, you need split the 2 sets of data and store them seprately, as this model uses 2 neuronal popn. types.
    
    test_Bs = np.squeeze(Bs)
    test_Ss = np.squeeze(Ss)
    
    df_B_e = pd.DataFrame(np.squeeze(test_Bs[:,0,:]),index=Bts)
    df_S_e = pd.DataFrame(np.squeeze(test_Ss[:,0,:]),index=Sts)
    
    df_B_i = pd.DataFrame(np.squeeze(test_Bs[:,1,:]),index=Bts)
    df_S_i = pd.DataFrame(np.squeeze(test_Ss[:,1,:]),index=Sts)
    
    # Remove transient time
    
    _test_Bs_e = Bs[10:int(simlen/200),0,:,:]
    _test_Bs_i = Bs[10:int(simlen/200),1,:,:]

    _test_Bts = Bts[10:int(simlen/200)]
    
    # Build a TimeSeries Datatype
    tsr_e = TimeSeriesRegion(connectivity=con,
                       data=_test_Bs_e, 
                       sample_period=sim.monitors[0].period)

    tsr_i = TimeSeriesRegion(connectivity=con,
                       data=_test_Bs_i, 
                       sample_period=sim.monitors[0].period)
    tsr_e.configure()
    tsr_i.configure()
    
    # Compute FC
    FC_e = np.corrcoef(np.squeeze(tsr_e.data).T)
    FC_i = np.corrcoef(np.squeeze(tsr_i.data).T) 
#     savemat('FC_' + str(G) + '_' + str(simlen) + '.mat', {'B': Bs, 'FC': FC})

    # Take triangular upper part of connectivity matrices and compute pearson correlations
    pcc_FC_e = np.corrcoef(HCP_FC[mask], FC_e[mask])[0, 1]
    pcc_SC_e = np.corrcoef(HCP_SC[mask], FC_e[mask])[0, 1]

    pcc_FC_i = np.corrcoef(HCP_FC[mask], FC_i[mask])[0, 1]
    pcc_SC_i = np.corrcoef(HCP_SC[mask], FC_i[mask])[0, 1]
    
    #return pcc
    return pcc_FC_e, pcc_SC_e, df_B_e, df_S_e, tsr_e, pcc_FC_i, pcc_SC_i, df_B_i, df_S_i, tsr_i


print('Model setup complete')

### ========================================================================================================================================

Davide_path = '/external/rprshnas01/netdata_kcni/jglab/Data/Shrey/Neg_Corrs_Jam_Session_June_2022/Davide_HCP_Data_Matrix'

ShreyVB_path = '/external/rprshnas01/netdata_kcni/jglab/Data/Shrey/Neg_Corrs_Jam_Session_June_2022'

Wts_Path = '/external/rprshnas01/netdata_kcni/jglab/Data/Shrey/Improved_WWD_HCP_model_runs/Sub_Specific_SC_Wts_2'

pconn_path = '/external/rprshnas01/netdata_kcni/jglab/Data/Shrey/Shrey_SS_parcellated_Func_Conns_II/'

# _sub_list = [200614, 199958, 177746, 164131, 141826, 130619, 127933, 116726, 100307, 100206]
# sub_list = list(reversed(_sub_list))

parcs = np.arange(0,200,1)

# -------------------------------------------------------------------------------------------------------------------------------------

# 5 Hz
regime = {'J_N': 0.18, 'J_i': 1.03, 'W_e': 1, 'W_i': 0.9} 

# -------------------------------------------------------------------------------------------------------------------------------------

# D1 = 0.01
# mynoise = noise.Additive(nsig=__((D1**2)/2),noise_seed=np.random.randint(0,100))

print('Params set for 5 Hz r_e/r_i')

### ========================================================================================================================================

# SC_wts

_gg = np.loadtxt(Wts_Path + "/{0}_SC_wts.txt".format(sub_id)) 

# Structural Connectivity
HCP_SC1 = _gg[parcs][:,parcs] #_step2a.copy() 
HCP_SC1 = HCP_SC1 + HCP_SC1.T # --> Symmetric
HCP_SC = HCP_SC1.copy()

print('SC setup complete')

### ========================================================================================================================================

# Get pconns
pconn1LR = pconn_path + '{0}_rfMRI_REST1_LR_Schaefer200_cifti_correlated.pconn.nii'.format(sub_id)
pconn1RL = pconn_path + '{0}_rfMRI_REST1_RL_Schaefer200_cifti_correlated.pconn.nii'.format(sub_id)
pconn2LR = pconn_path + '{0}_rfMRI_REST2_LR_Schaefer200_cifti_correlated.pconn.nii'.format(sub_id)
pconn2RL = pconn_path + '{0}_rfMRI_REST2_RL_Schaefer200_cifti_correlated.pconn.nii'.format(sub_id)

# Load pconns
_pconn_img1LR = nib.load(pconn1LR)
_pconn_dat1LR = _pconn_img1LR.get_data()
_pconn_dat1LR = _pconn_dat1LR/1

_pconn_img1RL = nib.load(pconn1RL)
_pconn_dat1RL = _pconn_img1RL.get_data()
_pconn_img2LR = nib.load(pconn2LR)
_pconn_dat2LR = _pconn_img2LR.get_data()
_pconn_img2RL = nib.load(pconn2RL)
_pconn_dat2RL = _pconn_img2RL.get_data()
_test_pconn = (_pconn_dat1LR + _pconn_dat1LR + _pconn_dat2LR + _pconn_dat2RL)/4

# An alternate version of this is to use only 1 resting-state run for each subject ... more neg corrs. Averaging like we've done above reduces overall number of neg corrs. 

HCP_FC1 = _test_pconn.copy()
HCP_FC = HCP_FC1[parcs][:,parcs]

# for a single rs scan ... decided to go with REST1_LR ...

HCP_FC_single_rs = _pconn_dat1LR.copy()


print('FC setup done!')

### ========================================================================================================================================

# Tract Lengths -->  irrelevant ... doing it so it doesn't give error, lazy to correct, but not lazy to type this comment apparently!

tract_lengths = np.loadtxt(ShreyVB_path + "/tract_lengths.txt")
tract_lengths = tract_lengths[parcs][:,parcs]

# Centroids
centroids = np.loadtxt(ShreyVB_path + "/all_centroids.txt")
centroids = centroids[parcs,:]


# Labels
f = open(Davide_path + "/labels.txt", "r")
q = f.read()
q = q.rsplit('\n')
del q[-32:]
labels = q.copy()
labels = np.array(labels)
labels = labels[parcs]

nregions = len(parcs)

mask = np.tril_indices(len(parcs), -1)

print("Loaded other Stuff!")

### ========================================================================================================================================

# Build HCP_con

HCP_con = connectivity.Connectivity()

HCP_con.tract_lengths = tract_lengths.copy()
HCP_con.speed = np.array(np.inf)
HCP_con.weights = HCP_SC.copy()
HCP_con.centres = centroids.copy()
HCP_con.areas = np.empty(0)
HCP_con.cortical = np.empty(0,dtype=bool)
HCP_con.orientations = np.empty(0)
HCP_con.hemispheres = np.concatenate((np.ones(100, dtype=bool),np.zeros(100,dtype=bool)))
HCP_con.region_labels = np.array(labels)
# HCP_con.delays = conduction_delays.copy()

HCP_con.weights = -np.diag((HCP_con.weights/np.linalg.norm(HCP_con.weights)).sum(0)) + HCP_con.weights/np.linalg.norm(HCP_con.weights)

HCP_con.configure()

print('Built HCP_con!')

### ========================================================================================================================================

print('Running model ... ')

Gs = np.arange(20, 60, 0.5)
# Gs = np.arange(40, 45, 0.1)
# Gs = np.arange(37, 39, 1)

test_pcc_FC_e = [None]*len(Gs)
test_pcc_SC_e = [None]*len(Gs)
test_df_B_e = [None]*len(Gs)
test_df_S_e = [None]*len(Gs)
test_tsr_e = [None]*len(Gs)

test_pcc_FC_i = [None]*len(Gs)
test_pcc_SC_i = [None]*len(Gs)
test_df_B_i = [None]*len(Gs)
test_df_S_i = [None]*len(Gs)
test_tsr_i = [None]*len(Gs)

for iG, G in enumerate(Gs):
    D1 = 0.01
    mynoise = noise.Additive(nsig=__((D1**2)/2),noise_seed=np.random.randint(0,100))
    if(iG == 0):
        print("started ...")
    elif(iG == round((len(Gs)/2))):
        print("Half-way ...")        
    test_pcc_FC_e[iG], test_pcc_SC_e[iG], test_df_B_e[iG], test_df_S_e[iG], test_tsr_e[iG], test_pcc_FC_i[iG], test_pcc_SC_i[iG], test_df_B_i[iG], test_df_S_i[iG], test_tsr_i[iG] = run_rww_sim_pcc_test(HCP_con, G, regime, 0.5, 300000, mynoise)
    
print("Finished Running Gs sweep !")    

### ========================================================================================================================================

# PCC for df_S_e vs. emp FC

FC_sim_best_test = []
FC_sim_best_test_single_rs = []

for i in range(0,len(Gs)):
    
    B = test_df_S_e[i].copy()
    FC_sim = np.corrcoef(B.T)
    
    pcc_FC_test = np.corrcoef(FC_sim[mask], HCP_FC[mask])[0,1]
    FC_sim_best_test.append(pcc_FC_test)
    
    pcc_FC_test_single_rs = np.corrcoef(FC_sim[mask], HCP_FC_single_rs[mask])[0,1]
    FC_sim_best_test_single_rs.append(pcc_FC_test_single_rs)
    
    
# get max value of PCC

max_pcc = 0
idx = 0

for i in range(len(FC_sim_best_test)):
    if max_pcc < FC_sim_best_test[i]:
        max_pcc = FC_sim_best_test[i].copy()
        idx = i

print(max_pcc)
print ("max pcc (with df_S_e) for subj {0} is {1}, with Gs = {2}".format(sub_id, max_pcc, Gs[idx]))

print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")


max_pcc_single_rs = 0
idx_single_rs = 0

for i in range(len(FC_sim_best_test_single_rs)):
    if max_pcc_single_rs < FC_sim_best_test_single_rs[i]:
        max_pcc_single_rs = FC_sim_best_test_single_rs[i].copy()
        idx_single_rs = i

print(max_pcc_single_rs)
print ("max pcc for REST1_LR ONLY(with df_S_e) for subj {0} is {1}, with Gs = {2}".format(sub_id, max_pcc_single_rs, Gs[idx_single_rs]))

### ========================================================================================================================================


### ========================================================================================================================================

print("Time taken to complete : --- %s seconds ---" % (time.time() - start_time))

### ========================================================================================================================================
