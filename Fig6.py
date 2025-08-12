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

seed(4321)
np.random.seed(4321)

duration        = 2*second
num_neurons     = 5

area            = 20000*umetre**2
Cm              = 1*ufarad*cm**-2 * area
gl              = 5e-5*siemens*cm**-2 * area
El              = -60*mV
EK              = -90*mV
ENa             = 50*mV
g_na            = 100*msiemens*cm**-2 * area
g_kd            = 30*msiemens*cm**-2 * area
VT              = -63*mV

lambda_values   = np.linspace(0.75,2,6)

nb_samples = 20

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler('color', plt.get_cmap('viridis',len(lambda_values)).colors)
fig,axs=plt.subplots(nrows=1,ncols=1,sharex=True,figsize=(6.4,4))

results = np.zeros((nb_samples,len(lambda_values),num_neurons))

for idx_sample in range(nb_samples):
    print(idx_sample)
    
    rand_values = uniform(0.9,1.1,9)
    
    area            = area*rand_values[0]
    Cm              = Cm*rand_values[1]
    gl              = gl*rand_values[2]
    El              = El*rand_values[3]
    EK              = EK*rand_values[4]
    ENa             = ENa*rand_values[5]
    g_na            = g_na*rand_values[6]
    g_kd            = g_kd*rand_values[7]
    VT              = VT*rand_values[8]

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
        group.I = '0.7*nA * (i+1) / num_neurons'
        
        monitor = SpikeMonitor(group)
        
        run(duration)
        
        results[idx_sample,idx,:] = monitor.count / duration
    
# Plots #######################################################################

barWidth = 0.15
br = arange(num_neurons)
for i in range(len(lambda_values)):
           
    bar(br, np.mean(results[:,i,:],axis=0), width = barWidth,
            edgecolor ='grey', label = r'$\lambda = $' + str(lambda_values[i]),
            yerr = np.std(results[:,i,:],axis=0)/np.sqrt(nb_samples))

    br_prev = br
    br = [x + barWidth for x in br_prev]
    
    xlabel('I [nA]')
    ylabel(r'Firing rate [$s^{-1}$]')
    plt.xticks(np.arange(num_neurons)+2.5*barWidth, [round(0.7 * (i+1) / num_neurons,2) for i in range(num_neurons)])
    legend()
    grid()
savefig("Fig5.png", dpi=300)  

# plt.figure()
# for i in range(len(lambda_values)):
#     plt.plot(results[16,i,:])
 
