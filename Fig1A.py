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

# Parameters ##################################################################

num_neurons     = 100

duration        = 2*second

area            = 20000*umetre**2
Cm              = 1*ufarad*cm**-2 * area
gl              = 5e-5*siemens*cm**-2 * area
El              = -60*mV
EK              = -90*mV
ENa             = 50*mV
g_na            = 100*msiemens*cm**-2 * area
g_kd            = 30*msiemens*cm**-2 * area
VT              = -63*mV

lambda_values = np.linspace(0.75,2,6)
fr_res = np.zeros(len(lambda_values))

# Loop over the values of lambda ##############################################
from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('viridis',len(lambda_values)).colors)
fig,axs=plt.subplots(nrows=1,ncols=1,sharex=True,figsize=(6.4,4))

for idx, lambd in enumerate(lambda_values):
    
    # The model ###############################################################

    eqs = Equations('''
    dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + I)/Cm : volt
    dm/dt = lambd*(0.32*(mV**-1)*4*mV/exprel((13.*mV-v+VT)/(4*mV))/ms*(1-m)-0.28*(mV**-1)*5*mV/exprel((v-VT-40.*mV)/(5*mV))/ms*m) : 1
    dn/dt = lambd*(0.032*(mV**-1)*5*mV/exprel((15.*mV-v+VT)/(5*mV))/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n) : 1
    dh/dt = lambd*(0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h) : 1
    P_K = n*n*n*n : 1
    P_Na = (m*m*m)*h : 1
    I : amp
    ''')
    # Threshold and refractoriness are only used for spike counting
    group = NeuronGroup(num_neurons, eqs,
                        threshold='v > -50*mV',
                        refractory='v > -50*mV',
                        method='exponential_euler')
    group.v = El
    group.I = '0.7*nA * i / num_neurons'
    
    monitor = SpikeMonitor(group)
    
    run(duration)
    
    # Verification that, for a given current, the FR increases with lambda ####
    
    
    axs.plot(group.I/nA, monitor.count / duration, label = r'$\lambda = $' + str(lambd))
    fr_res[idx] = (monitor.count[-1] / duration)*second
    
left, bottom, width, height = [0.65, 0.25, 0.28, 0.28]
ax2  = fig.add_axes([left, bottom, width, height])
ax2.plot(lambda_values,fr_res, color='k')
ax2.axvline(x=1.0, color='k', linestyle='--')
ax2.set_xlabel(r'$\lambda$ [-]')
ax2.set_ylabel(r'Firing rate [$s^{-1}$]')
ax2.grid()

axs.set_xlabel('I [nA]')
axs.set_ylabel(r'Firing rate [$s^{-1}$]')
axs.legend()
axs.grid()
    # show()

fig.tight_layout()
savefig("Fig2A.png", dpi=300)   
