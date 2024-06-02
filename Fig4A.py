# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:28:36 2024

@author: CAG329
"""

# Packages ####################################################################

from brian2 import *
# %matplotlib qt

# Parameters ##################################################################

I_lambda = [0*pA,10*pA]

seed(4321)
np.random.seed(4321)

duration = 20*second
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

fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True)

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

        
    # Results #####################################################################
    
    axs[lambda_idx].plot(0.001*spikemon.t/ms, spikemon.i, '.k', markersize=0.5)
    
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Units')
    axs[0].set_ylabel('Units')
    
fig.tight_layout()
savefig("Fig4A.png", dpi=200) 