#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 16:07:14 2025

@author: camille
"""

# Packages ####################################################################

from brian2 import *
matplotlib.rcParams.update({'font.size': 12})
# %matplotlib qt

# Parameters ##################################################################

n_samples = 10
we_lambda = [-0.1*mV,-0.05*mV,0.0*mV,0.05*mV,0.1*mV]

seed(4321)
np.random.seed(4321)

duration = 50*second
# Neuronal population
N_E = 1000
# Total capacitance
C = 200*pF 
# Total leak capacitance
gL = 10*nS
# Effective rest potential
EL = -60*mV
# Threshold slope factor
DeltaT = 2*mV
# Effective threshold potential
VT = -50*mV
# Conductance
a = 2*nS
# Time constant
tauw = 800*msecond
# Adaptation
b = 30*pA
# Reset potential
Vr = -46*mV
# External input
I_values = 102*pA

# Synaptic parameters
connect_prob = 0.02
we_values = 0.5*mV
sigma = 1*mV
tau_E = 0.01*second

res_FR = np.zeros((N_E,n_samples,len(we_lambda)))
res_bursts = np.zeros((N_E,n_samples,len(we_lambda)))

# Example voltages
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('viridis',len(we_lambda)).colors)
fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True,sharey=True)

for n in range(n_samples):
    print(n)
    for we_idx,we_lamb in enumerate(we_lambda):

# Equations ###################################################################
    
        eqs_E = Equations('''
        dV/dt = (-gL*(V-EL) + gL*DeltaT*exp((V-VT)/DeltaT) + I - w)/C  + sigma*sqrt(2/tau_E)*xi : volt
        dw/dt = (a*(V-EL) - w)/tauw : amp
        I : amp
        we: volt
        ''')
        
        G_E = NeuronGroup(N_E, 
                        model=eqs_E,
                        threshold='V > 0*mV', 
                        reset='V=Vr; w+=b',
                        )
        
        G_E.I = I_values
        G_E.we = we_values + we_lamb
        
        S_EE = Synapses(G_E, G_E, on_pre='V += we')
        
        S_EE.connect(p=connect_prob)
        
        spikemon = SpikeMonitor(G_E)
        M = StateMonitor(G_E, 'V', record=True)
        P = StateMonitor(G_E, 'w', record=True)
        
        G_E.V = 'EL+EL*rand() * 0.1 '
        
        run(duration)
                
        res_FR[:,n,we_idx] = spikemon.count / duration
        
        spike_trains = spikemon.spike_trains()
        for i in range(N_E):
            spike_trains_ = spike_trains[i]/msecond
            is_in_burst_idx = np.where(spike_trains_[1:] - spike_trains_[:-1] < 20)[0]
            if len(is_in_burst_idx)>0:
                res_bursts[i,n,we_idx] = (len(np.where(np.diff(is_in_burst_idx)>1)[0])+1) / duration
                
        # Results #####################################################################
        
        if n == 0:
            axs[0].plot(0.001*M.t/ms,M.V[0]/mV, label = r'$w_{e}^{\lambda} = $' + str(we_lamb))
            axs[0].axhline(y=VT/mV,linestyle='--',c='k')
            axs[1].axhline(y=VT/mV,linestyle='--',c='k')
    
            axs[1].plot(0.001*M.t/ms,M.V[1]/mV)
            axs[1].set_xlabel('Time [s]')
            axs[0].set_ylabel('V [mV]')
            axs[1].set_ylabel('V [mV]')
            axs[0].grid()
            axs[1].grid()
            axs[0].legend(loc='lower left')
            axs[0].set_xlim([2.0, 2.84])
            axs[0].set_ylim([-70, 5])
            fig.tight_layout()
            savefig("Fig6B.svg", dpi=300) 
            
# Results #####################################################################

fig,axs=plt.subplots(nrows=1,ncols=2,sharex=True)
axs[0].plot(we_lambda/mV,np.mean(res_FR,axis=(0,1)), c='tab:blue')
axs[0].fill_between(we_lambda/mV,
                      np.mean(res_FR,axis=(0,1))+np.std(res_FR,axis=(0,1))/np.sqrt((n_samples)), 
                      np.mean(res_FR,axis=(0,1))-np.std(res_FR,axis=(0,1))/np.sqrt((n_samples)), alpha=0.2, color='tab:blue')

axs[1].plot(we_lambda/mV,np.mean(res_bursts,axis=(0,1)), c='tab:blue')
axs[1].fill_between(we_lambda/mV,
                      np.mean(res_bursts,axis=(0,1))+np.std(res_bursts,axis=(0,1))/np.sqrt((n_samples)), 
                      np.mean(res_bursts,axis=(0,1))-np.std(res_bursts,axis=(0,1))/np.sqrt((n_samples)), alpha=0.2, color='tab:blue')


axs[0].set_xlabel(r'$w_{e}^{\lambda}$ [mV]' )
axs[1].set_xlabel(r'$w_{e}^{\lambda}$ [mV]' )
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
fig.tight_layout()
savefig("Fig6A.svg", dpi=300) 

