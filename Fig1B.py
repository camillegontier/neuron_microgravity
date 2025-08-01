# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 19:49:52 2023

@author: CAG329
"""

# From https://brian2.readthedocs.io/en/stable/user/plotting_functions.html

# Packages ####################################################################

from brian2 import *
matplotlib.rcParams.update({'font.size': 14})

# Parameters ##################################################################

VT = -63*mV
lambda_values = [1,1.5]

# Results #####################################################################

fig,axs=plt.subplots(nrows=1,ncols=2,sharex=True, figsize=(8,4))
linestyles = ['-','--']
linewidths = [1,3]

for lambda_idx, lambda_value in enumerate(lambda_values):

    eq = Equations("""
    alpha_m = lambda_value*(0.32*(mV**-1)*4*mV/exprel((13*mV-v+VT)/(4*mV))/ms) : Hz
    beta_m = lambda_value*(0.28*(mV**-1)*5*mV/exprel((v-VT-40*mV)/(5*mV))/ms) : Hz
    alpha_h = lambda_value*(0.128*exp((17*mV-v+VT)/(18*mV))/ms) : Hz
    beta_h = lambda_value*(4./(1+exp((40*mV-v+VT)/(5*mV)))/ms) : Hz
    alpha_n = lambda_value*(0.032*(mV**-1)*5*mV/exprel((15*mV-v+VT)/(5*mV))/ms) : Hz
    beta_n = lambda_value*(.5*exp((10*mV-v+VT)/(40*mV))/ms) : Hz
    tau_n = 1/(alpha_n + beta_n) : second
    tau_m = 1/(alpha_m + beta_m) : second
    tau_h = 1/(alpha_h + beta_h) : second
    n_infty = alpha_n/(alpha_n + beta_n) : 1
    m_infty = alpha_m/(alpha_m + beta_m) : 1
    h_infty = alpha_h/(alpha_h + beta_h) : 1
    """)
    
    group = NeuronGroup(100, eq + Equations("v : volt"))
    group.v = np.linspace(-100, 100, len(group))*mV
    
    linestyle = linestyles[lambda_idx]
    linewidth = linewidths[lambda_idx]
    
    axs[0].plot(group.v/mV, group.tau_m[:]/ms, label=r"$\tau_m (V)$",linestyle=linestyle,c='tab:blue',linewidth = linewidth)
    axs[0].plot(group.v/mV, group.tau_n[:]/ms, label=r"$\tau_n (V)$",linestyle=linestyle,c='tab:orange',linewidth = linewidth)
    axs[0].plot(group.v/mV, group.tau_h[:]/ms, label=r"$\tau_h (V)$",linestyle=linestyle,c='tab:green',linewidth = linewidth)
    axs[0].set_xlabel('[mV]')
    axs[0].set_ylabel('[ms]')
    if lambda_idx == 0:
        axs[0].legend()
    
    axs[1].plot(group.v/mV, group.m_infty[:], label=r"$m_{\infty} (V)$",linestyle=linestyle,c='tab:blue',linewidth = linewidth)
    axs[1].plot(group.v/mV, group.n_infty[:], label=r"$n_{\infty} (V)$",linestyle=linestyle,c='tab:orange',linewidth = linewidth)
    axs[1].plot(group.v/mV, group.h_infty[:], label=r"$h_{\infty} (V)$",linestyle=linestyle,c='tab:green',linewidth = linewidth)
    axs[1].set_xlabel('[mV]')
    if lambda_idx == 0:
        axs[1].legend()
    
axs[0].grid()
axs[1].grid()
fig.tight_layout()
plt.savefig('Fig1B.png',dpi=300)