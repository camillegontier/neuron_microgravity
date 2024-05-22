# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 07:50:01 2023

@author: CAG329
"""

# Original code from
# https://brian2.readthedocs.io/en/stable/examples/IF_curve_Hodgkin_Huxley.html

# Packages ####################################################################

from brian2 import *
# %matplotlib qt

# Parameters ##################################################################

duration        = 2*second
num_neurons     = 100

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

# Loop over the values of lambda ##############################################

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
    M_K = StateMonitor(group, 'P_K', record=True)
    M_Na = StateMonitor(group, 'P_Na', record=True)
    
    run(duration)
    
    # Effect of lambda on channels opening probabilities ######################
    
    jj = round(num_neurons/2)
    
    from cycler import cycler
    plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('viridis',len(lambda_values)).colors)
    
    tt = np.arange(0,len(M_K.P_K[jj,:]))*1e-4

    
    subplot(221)
    plot(tt,M_K.P_K[jj,:],label= r'$\lambda = $' + str(lambd))
    xlabel("Time [s]")
    xlim(0, 0.07)
    ylabel(r'$P_K$ [-]')
    title('I = ' + str(0.7*(jj/num_neurons)) + ' [nA]')
    # legend()
    grid()
    
    subplot(222)
    plot(tt,M_Na.P_Na[jj,:],label= r'$\lambda = $' + str(lambd))
    xlabel("Time [s]")
    xlim(0, 0.07)
    ylabel(r'$P_{Na}$ [-]')
    title('I = ' + str(0.7*(jj/num_neurons)) + ' [nA]')
    grid()
    
    subplot(223)
    plot(group.I/nA,np.mean(M_K.P_K,axis=1)/monitor.count)
    yticks([])
    xlabel('I [nA]')
    ylabel(r'<$P_{K}>$ [-]')
    grid()
    
    subplot(224)
    plot(group.I/nA,np.mean(M_Na.P_Na,axis=1)/monitor.count)
    yticks([])
    xlabel('I [nA]')
    ylabel(r'<$P_{Na}>$ [-]')
    grid()
    
    tight_layout()
    
savefig("Fig2C.png", dpi=300)     