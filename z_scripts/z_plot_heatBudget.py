import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import base
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

ld = r"i_data_literature/"
dd = r"i_data_processed/"
fd = r"o_finalPlots/"

wt_ls = [8.2, 20.07] # based off liquidus alignment - 0.5K
mass_ls = [4.1943,3.7107]

## plot energy decomp

fig, ax = plt.subplots(1,2,sharey=True)
idx = 0
wt = wt_ls[idx]
m = mass_ls[idx]
data = pd.read_csv(dd+ rf'{wt}wt%_hb_{m}g.csv', header=0)

ax[0].stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])

ax[0].set_ylim([0, 1])
ax[0].set_xlim(184.5,233.5)

ax[0].text(0.5, 0.10, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.5, 0.41, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.5, 0.81, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.72, 0.687, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax[0].transAxes)
ax[0].arrow(220,0.685,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax[0].text(0.05, 0.94, '8.2 wt%', fontsize=4, horizontalalignment='left',
     verticalalignment='center', transform=ax[0].transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))




# handles, labels = ax[0].get_legend_handles_labels()
# legend = ax[0].legend(handles[::-1], labels[::-1],fontsize=4, loc='upper left',frameon=True,fancybox=False,framealpha=1,edgecolor='k')
# legend.get_frame().set_linewidth(.3)

idx = 1
wt = wt_ls[idx]
m = mass_ls[idx]
data = pd.read_csv(dd+ rf'{wt}wt%_hb_{m}g.csv', header=0)

ax[1].stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])
ax[1].set_ylim([0, 1])
ax[1].set_xlim(184.5,227.5)

ax[1].text(0.5, 0.145, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.352, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.68, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.72, 0.487, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax[1].transAxes)
ax[1].arrow(216,0.485,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax[1].text(0.05, 0.94, '20.07 wt%', fontsize=4, horizontalalignment='left',
     verticalalignment='center', transform=ax[1].transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))


base.splt_axis_label(fig,'Temperature (K)','Heat (Normalized)')
plt.tight_layout()

# base.show_plot_max()
saveFig = fd + rf'energy_budget.png'
plt.savefig(saveFig)
plt.close()