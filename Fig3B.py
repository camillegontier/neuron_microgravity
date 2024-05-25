# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 18:42:59 2023

@author: CAG329
"""

# Packages ####################################################################

from brian2 import *
# %matplotlib qt

# Parameters ##################################################################

duration = 4*second

# Neuronal population
N_E = 100

# Strength
we = (0.67*2)*mV # excitatory synaptic weight (voltage)

# Membrane time constants
tau_E = 0.01*second
tau_A = 0.15*second

g_A = 10*mV

threshold = 20*mV
reset = 0*volt

# External input
I_E_values = 21*mV

sigma = 1*mV

plot_offset = 70

# Equations ###################################################################

eqs_E = Equations('''
dV_E/dt = (I_E - V_E - A_E)/tau_E + sigma*sqrt(2/tau_E)*xi : volt
dA_E/dt = - A_E/tau_A : volt
I_E : volt
''')

G_E = NeuronGroup(N_E, 
                model=eqs_E,
                threshold='V_E > threshold', 
                reset='V_E=reset; A_E+=g_A',
                refractory=3*ms)
G_E.I_E = I_E_values

S_EE = Synapses(G_E, G_E, on_pre='V_E += we')

C = np.ones(N_E,dtype=int) - np.identity(N_E,dtype=int)  # The connection matrix as a numpy array of 0's and 1's
sources, targets = C.nonzero()
# S_EE.connect(i=sources, j=targets)
S_EE.connect(p=0.5)

spikemon = SpikeMonitor(G_E)
M = StateMonitor(G_E, 'V_E', record=True)
P = StateMonitor(G_E, 'A_E', record=True)

G_E.V_E = 'rand() * 20*mV '

run(duration)

# Results #####################################################################

fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True)
axs[0].plot(M.t/ms,M.V_E[0]/mV-plot_offset)
axs[0].axhline(y=threshold/mV-plot_offset,linestyle='--',c='k')
axs[1].axhline(y=threshold/mV-plot_offset,linestyle='--',c='k')

axs[1].plot(M.t/ms,M.V_E[1]/mV-plot_offset)
axs[1].set_xlabel('Time [ms]')
axs[0].set_ylabel('V [mV]')
axs[1].set_ylabel('V [mV]')
axs[0].grid()
axs[1].grid()
fig.tight_layout()
savefig("Fig3B.png", dpi=300) 

# plt.figure()
# plot(spikemon.t/ms, spikemon.i, '.k')
# xlabel('Time (ms)')
# ylabel('Neuron index')

# plt.figure()
# plot(M.t/ms, M.V_E[1]/mV)
# plot(M.t/ms, M.V_E[2]/mV)
# plot(M.t/ms, M.V_E[3]/mV)
# plot(M.t/ms, M.V_E[4]/mV)
# plot(M.t/ms, M.V_E[5]/mV,linestyle='--')
# xlabel('t (ms)')
# ylabel('v (mV)')
# show() 

# plt.figure()
# plot(M.t/ms, P.A_E[1]/mV)
# xlabel('t (ms)')
# ylabel('v (mV)')
# show() 
