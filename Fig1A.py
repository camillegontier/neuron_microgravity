# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 09:01:57 2023

@author: CAG329
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.array([23,24,25,26,27,28,29,30,31,32,33,34,35,36,37])
ya = np.array([61.4,64.96,73.64,89.88,92.55,103.0,103.67,107.01,
      106.56,107.01,108.34,109.45,111.01,113.01,116.35])
stda = np.array([66.96,69.41,78.09,93.44,97.44,107.23,109.01,112.12,
        111.46,111.01,110.34,112.12,113.24,115.46,118.35])
origin_a = 121.02
ratio_a = 4/np.abs(49.61-origin_a)
stda = np.abs(ya-stda)*ratio_a

ya = np.abs(origin_a-ya)*ratio_a
yb = np.array([141.94,146.16,151.5,152.84,158.4,160.18,163.96,167.07,170.19,
      174.19,176.86,179.75,183.09,188.43,189.54])
stdb = np.array([145.49,149.94,154.39,156.17,160.4,162.85,165.74,169.52,
        172.41,175.75,179.09,182.65,185.76,189.77,191.55])
origin_b = 208.68
ratio_b = 5/np.abs(137.26-origin_b)
stdb = np.abs(yb-stdb)*ratio_b

yb = np.abs(origin_b-yb)*ratio_b


fig,axs=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(4.8,4.8))
axs[0].plot(x,ya)
axs[0].fill_between(x,ya-stda,
                      ya+stda,
                      alpha = 0.5)
axs[0].set_ylabel(r"$\tau_m$ [ms]")
axs[0].grid()
axs[1].plot(x,yb)
axs[1].fill_between(x,yb-stdb,
                      yb+stdb,
                      alpha = 0.5)
axs[1].set_ylabel(r"$\tau_h$ [ms]")
axs[1].set_xlabel(r"Temperature [$^{\circ}C$]")
axs[1].grid()
fig.tight_layout()
plt.savefig('Fig1A.png',dpi=300)


