# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 07:49:50 2023

@author: CAG329
"""

# Original code from
# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html

# Packages ####################################################################

from brian2 import *
# %matplotlib qt
matplotlib.rcParams.update({'font.size': 13})

# Parameters ##################################################################

seed(4321)
np.random.seed(4321)


duration        = 1*second
num_neurons     = 4000

area = 20000*umetre**2
Cm = (1*ufarad*cm**-2) * area
gl = (5e-5*siemens*cm**-2) * area
El = -60*mV
EK = -90*mV
ENa = 50*mV
g_na = (100*msiemens*cm**-2) * area
g_kd = (30*msiemens*cm**-2) * area
VT = -63*mV
# Time constants
taue = 5*ms
taui = 10*ms
# Reversal potentials
Ee = 0*mV
Ei = -80*mV
we = 6*nS  # excitatory synaptic weight
wi = 67*nS  # inhibitory synaptic weight

# lambda_values   = np.linspace(1,2,5)
EI_ratio_values = np.repeat(0.8,30)
lambda_values   = np.linspace(0.5,1.5,7)
connect_prob = 0.02

results = np.zeros((num_neurons,len(lambda_values),len(EI_ratio_values)))
results_bursts = np.zeros((num_neurons,len(lambda_values),len(EI_ratio_values)))

# Loop over the values of lambda ##############################################

for idx_lambd, lambd in enumerate(lambda_values):
    for idx_EI_ratio, EI_ratio in enumerate(EI_ratio_values):
        print(idx_EI_ratio)
    
    # The model ###############################################################
    
        eqs = Equations('''
        dv/dt = (gl*(El-v)+ge*(Ee-v)+gi*(Ei-v)-
                 g_na*(m*m*m)*h*(v-ENa)-
                 g_kd*(n*n*n*n)*(v-EK))/Cm : volt
        dm/dt = alpha_m*(1-m)-beta_m*m : 1
        dn/dt = alpha_n*(1-n)-beta_n*n : 1
        dh/dt = alpha_h*(1-h)-beta_h*h : 1
        dge/dt = -ge*(1./taue) : siemens
        dgi/dt = -gi*(1./taui) : siemens
        alpha_m = lambd*(0.32*(mV**-1)*4*mV/exprel((13*mV-v+VT)/(4*mV))/ms) : Hz
        beta_m = lambd*(0.28*(mV**-1)*5*mV/exprel((v-VT-40*mV)/(5*mV))/ms) : Hz
        alpha_h = lambd*(0.128*exp((17*mV-v+VT)/(18*mV))/ms) : Hz
        beta_h = lambd*(4./(1+exp((40*mV-v+VT)/(5*mV)))/ms) : Hz
        alpha_n = lambd*(0.032*(mV**-1)*5*mV/exprel((15*mV-v+VT)/(5*mV))/ms) : Hz
        beta_n = lambd*(.5*exp((10*mV-v+VT)/(40*mV))/ms) : Hz
        ''')
        # Threshold and refractoriness are only used for spike counting
        group = NeuronGroup(num_neurons, eqs,
                            threshold='v > -20*mV',
                            refractory=3*ms,
                            method='exponential_euler')
        
        Pe = group[:round(EI_ratio*num_neurons)]
        Pi = group[round(EI_ratio*num_neurons):]
        Ce = Synapses(Pe, group, on_pre='ge+=we')
        Ci = Synapses(Pi, group, on_pre='gi+=wi')
        Ce.connect(p=connect_prob)
        Ci.connect(p=connect_prob)
        
        group.ge = '(randn() * 1.5 + 4) * 10.*nS'
        group.gi = '(randn() * 12 + 20) * 10.*nS'
    
        group.v = 'El + (randn() * 5 - 5)*mV'
        # group.I = '0.0007*nA'
        
        monitor = SpikeMonitor(group)
        M = StateMonitor(group, 'v', record=True)
        
        run(duration)
        
        results[:,idx_lambd,idx_EI_ratio] = monitor.count / duration
        
        
        spike_trains = monitor.spike_trains()
        for i in range(num_neurons):
            spike_trains_ = spike_trains[i]/msecond
            is_in_burst_idx = np.where(spike_trains_[1:] - spike_trains_[:-1] < 20)[0]
            if len(is_in_burst_idx)>0:
                results_bursts[i,idx_lambd,idx_EI_ratio] = (len(np.where(np.diff(is_in_burst_idx)>1)[0])+1) / duration
            
fig,axs=plt.subplots(nrows=1,ncols=2,sharex=True,figsize=(6.4,4))

axs[0].plot(lambda_values,np.mean(results,axis=(0,2)),c='tab:blue')
axs[0].fill_between(lambda_values,
                      np.mean(results,axis=(0,2))+np.std(results,axis=(0,2))/np.sqrt(len(EI_ratio_values)), 
                      np.mean(results,axis=(0,2))-np.std(results,axis=(0,2))/np.sqrt(len(EI_ratio_values)), alpha=0.2,color='tab:blue')

axs[1].plot(lambda_values,np.mean(results_bursts,axis=(0,2)),c='tab:blue')
axs[1].fill_between(lambda_values,
                      np.mean(results_bursts,axis=(0,2))+np.std(results_bursts,axis=(0,2))/np.sqrt(len(EI_ratio_values)), 
                      np.mean(results_bursts,axis=(0,2))-np.std(results_bursts,axis=(0,2))/np.sqrt(len(EI_ratio_values)), alpha=0.2,color='tab:blue')


axs[0].set_xlabel(r'$\lambda$ [-]' )
axs[1].set_xlabel(r'$\lambda$ [-]' )
axs[0].set_title('Firing rates' )
axs[1].set_title('Burst rates' )
axs[0].set_ylabel(r'[$s^{-1}$]')
# axs[1].set_ylabel(r'[$s^{-1}$]')

# axs[0].set_title('p = ' + str(connect_prob))
# axs[0].legend()
axs[0].grid()
# axs[1].set_title('p = ' + str(connect_prob))
# axs[1].legend()
axs[1].grid()
axs[0].axvline(x=1.0, color='k', linestyle='--')
axs[1].axvline(x=1.0, color='k', linestyle='--')
fig.tight_layout()
savefig("Fig1D.svg", dpi=300) 
    

