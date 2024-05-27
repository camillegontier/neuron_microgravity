# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:28:36 2024

@author: CAG329
"""

# Packages ####################################################################

from brian2 import *
# %matplotlib qt

# Parameters ##################################################################

n_samples = 10
I_lambda = [0*pA,10*pA,20*pA,30*pA]

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
we = 0.5*mV
sigma = 1*mV
tau_E = 0.01*second

res_FR = np.zeros((N_E,n_samples,len(I_lambda)))
res_bursts = np.zeros((N_E,n_samples,len(I_lambda)))

for n in range(n_samples):
    print(n)
    for lambda_idx,I_lamb in enumerate(I_lambda):

# Equations ###################################################################
    
        eqs_E = Equations('''
        dV/dt = (-gL*(V-EL) + gL*DeltaT*exp((V-VT)/DeltaT) + I - w)/C  + sigma*sqrt(2/tau_E)*xi : volt
        dw/dt = (a*(V-EL) - w)/tauw : amp
        I : amp
        ''')
        
        G_E = NeuronGroup(N_E, 
                        model=eqs_E,
                        threshold='V > 0*mV', 
                        reset='V=Vr; w+=b',
                        )
        
        G_E.I = I_values+I_lamb
        
        S_EE = Synapses(G_E, G_E, on_pre='V += we')
        
        S_EE.connect(p=connect_prob)
        
        spikemon = SpikeMonitor(G_E)
        M = StateMonitor(G_E, 'V', record=True)
        P = StateMonitor(G_E, 'w', record=True)
        
        G_E.V = 'EL+EL*rand() * 0.1 '
        
        run(duration)
    
            
        res_FR[:,n,lambda_idx] = spikemon.count / duration
        
        
        spike_trains = spikemon.spike_trains()
        for i in range(N_E):
            spike_trains_ = spike_trains[i]/msecond
            is_in_burst_idx = np.where(spike_trains_[1:] - spike_trains_[:-1] < 20)[0]
            if len(is_in_burst_idx)>0:
                res_bursts[i,n,lambda_idx] = (len(np.where(np.diff(is_in_burst_idx)>1)[0])+1) / duration
            

# Results #####################################################################

fig,axs=plt.subplots(nrows=1,ncols=2,sharex=True)
axs[0].plot(I_lambda/pA,np.mean(res_FR,axis=(0,1)))
axs[0].fill_between(I_lambda/pA,
                      np.mean(res_FR,axis=(0,1))+np.std(res_FR,axis=(0,1))/np.sqrt((n_samples)), 
                      np.mean(res_FR,axis=(0,1))-np.std(res_FR,axis=(0,1))/np.sqrt((n_samples)), alpha=0.2)

axs[1].plot(I_lambda/pA,np.mean(res_bursts,axis=(0,1)))
axs[1].fill_between(I_lambda/pA,
                      np.mean(res_bursts,axis=(0,1))+np.std(res_bursts,axis=(0,1))/np.sqrt((n_samples)), 
                      np.mean(res_bursts,axis=(0,1))-np.std(res_bursts,axis=(0,1))/np.sqrt((n_samples)), alpha=0.2)


axs[0].set_xlabel(r'$I_{\lambda}$ [pA]' )
axs[1].set_xlabel(r'$I_{\lambda}$ [pA]' )
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
savefig("Fig4B.png", dpi=200) 