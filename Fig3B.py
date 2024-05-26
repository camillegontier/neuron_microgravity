# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 18:42:59 2023

@author: CAG329
"""

# Packages ####################################################################

from brian2 import *
# %matplotlib qt

# Parameters ##################################################################

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

# Equations ###################################################################

eqs_E = Equations('''
dV/dt = (-gL*(V-EL) + gL*DeltaT*exp((V-VT)/DeltaT) + I - w)/C : volt
dw/dt = (a*(V-EL) - w)/tauw : amp
I : amp
''')

G_E = NeuronGroup(N_E, 
                model=eqs_E,
                threshold='V > 0*mV', 
                reset='V=Vr; w+=b',
                )

G_E.I = I_values

S_EE = Synapses(G_E, G_E, on_pre='V += we')

S_EE.connect(p=connect_prob)

spikemon = SpikeMonitor(G_E)
M = StateMonitor(G_E, 'V', record=True)
P = StateMonitor(G_E, 'w', record=True)

G_E.V = 'EL+EL*rand() * 0.1 '

run(duration)

print(spikemon.count[0] / duration)
spike_trains = spikemon.spike_trains()
spike_trains_ = spike_trains[0]/msecond
is_in_burst_idx = np.where(spike_trains_[1:] - spike_trains_[:-1] < 20)[0]
if len(is_in_burst_idx)>0:
    print((len(np.where(np.diff(is_in_burst_idx)>1)[0])+1) / duration)
    
print(spikemon.count[1] / duration)
spike_trains = spikemon.spike_trains()
spike_trains_ = spike_trains[1]/msecond
is_in_burst_idx = np.where(spike_trains_[1:] - spike_trains_[:-1] < 20)[0]
if len(is_in_burst_idx)>0:
    print((len(np.where(np.diff(is_in_burst_idx)>1)[0])+1) / duration)
    
# Results #####################################################################

fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True)
axs[0].plot(M.t/ms,M.V[0]/mV)
axs[0].axhline(y=VT/mV,linestyle='--',c='k')
axs[1].axhline(y=VT/mV,linestyle='--',c='k')

axs[1].plot(M.t/ms,M.V[1]/mV)
axs[1].set_xlabel('Time [ms]')
axs[0].set_ylabel('V [mV]')
axs[1].set_ylabel('V [mV]')
axs[0].grid()
axs[1].grid()
fig.tight_layout()
savefig("Fig3B.png", dpi=200) 

