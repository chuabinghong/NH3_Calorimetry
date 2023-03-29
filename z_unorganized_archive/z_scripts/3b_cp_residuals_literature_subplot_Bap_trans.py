"""
created by Bing Hong CHUA 10Nov22

script objective:
plot deviations from Tillner-Roth and Friend's Helmholtz Equation (1998) and specific heat data by Giauque and others
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import base
from scipy.optimize import curve_fit

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
ccycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

ld = r"i_data_literature/"
dd = r"i_data_processed/"
od = r"i_data_processed/"
fd = r"t_all_cp/"

## SCRIPT ARGUMENTS
TRF_arg = 0

ramp_rate = 0.1
wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912]
mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]# based off liquidus alignment - 0.5K
colour_ls = ['#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84','#081d58']
colour_ls = ['#d73027','#f46d43','#fdae61','#e8e884','#abd9e9','#74add1','#4575b4']

## LOAD LITERATURE DATA

# df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
# df_TRF = pd.read_csv(ld + 'Tillner-Roth_Friend_liq.csv', header=0)
# df_St = pd.read_csv(ld + 'Steve_cp_old.csv', header=0)
# CG = {'33': pd.read_csv(ld + 'Chan_Giauque_0.33.csv', header=0)}
# HG = { '49': pd.read_csv(ld + 'Hildenbrand_Giauque_0.49.csv', header=0),
#            '59': pd.read_csv(ld + 'Hildenbrand_Giauque_0.59.csv', header=0),
#            '61': pd.read_csv(ld + 'Hildenbrand_Giauque_0.61.csv', header=0),
#            '64': pd.read_csv(ld + 'Hildenbrand_Giauque_0.64.csv', header=0),
#            '65': pd.read_csv(ld + 'Hildenbrand_Giauque_0.65.csv', header=0)}
# df_WK = pd.read_csv(ld + 'Wrewsky_Kaigorodoff.csv', header=0)

def import_TRF():
    array = np.loadtxt(ld + 'Baptiste_TRF.csv', delimiter=',')
    H2O = 18.01528
    NH3 = 17.03026
    T = np.arange(180,340+1,1)
    molfrac = np.arange(0,1.01,.01)
    mol_in_g = molfrac*NH3+(1-molfrac)*H2O
    wt = molfrac*NH3/mol_in_g
    interp = RegularGridInterpolator((T, wt), array,bounds_error=False,fill_value=np.NaN)
    return(interp)

TRF = import_TRF()

df_chan = pd.read_csv(ld + 'Baptiste/' + 'chan_1964.csv', header=0)
df_allred = pd.read_csv(ld + 'Baptiste/' + 'allred_1981.csv', header=0)
df_chernenkaya = pd.read_csv(ld + 'Baptiste/' + 'chernenkaya_1971.csv', header=0)
df_fujita = pd.read_csv(ld + 'Baptiste/' + 'fujita_2008.csv', header=0)
df_hildenbrand = pd.read_csv(ld + 'Baptiste/' + 'hildenbrand_1953.csv', header=0)
df_wrewsky = pd.read_csv(ld + 'Baptiste/' + 'wrewsky_1924.csv', header=0)

## PLOT FIG 9 OF TRF 1998

def calc_dev(T,wt,cp):
    TRF_val = TRF((T,wt))
    # dev = 100 * (1-TRF_val/cp)
    dev = 100 * (cp-TRF_val)/TRF_val
    return(dev)

plt.close('all')

w = 3.3
fig, ax = plt.subplots(2,1,sharey=True,figsize=(w,w*3/3))
plt.tight_layout()
base.splt_axis_label(fig,"",'$100(C_p{}^{data}$ - $C_p{}^{TF98}$) / $C_p{}^{TF98})$')

plt.subplots_adjust(hspace=0.5)
plt.sca(ax[0])
plt.ylim(-20, 20)
plt.axhline(0, linewidth=.5, color='k')
plt.xlim(0,100)
ax[0].set_xlabel('$NH_3$ Mass Fraction (wt%)')
plt.sca(ax[1])
plt.ylim(-20, 20)
plt.axhline(0, linewidth=.5, color='k')
plt.xlim(180,340)
plt.xlabel('Temperature (K)')

# plot experimental data
for idx,wt in enumerate(wt_ls):
    # idx=6
    # wt=26.912

    colour = colour_ls[idx]
    m = mass_ls[idx]
    data = base.DataCP(m,wt)
    data.import_data(dd, mean=0)
    dev = calc_dev(data.df_p['T(K)'], data.df_p['X'], data.df_p['cp(J/gK)']*1000)
    plt.sca(ax[0])
    plt.scatter([float(wt)] * len(dev), dev, 3,marker='o', facecolors="None",edgecolors=colour,linewidths=.5,label=f'{wt} wt%',zorder=20)
    plt.sca(ax[1])
    plt.scatter(data.df_p['T(K)'], dev, 3,marker='o',facecolors="None",edgecolors=colour,linewidths=.5,label=f'{wt} wt%',zorder=20)

ax[0].scatter(0,250,alpha=0)
ax[1].scatter(0,250,alpha=0)

# plot chan giauque 33 wt %
dev = calc_dev(df_chan['T(K)'],df_chan['massfrac'],df_chan['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_chan['massfrac']*100, dev, 3, marker='s',facecolors="None",edgecolors=ccycle[1],linewidths=.5,label='Chan & Giauque 1964')
plt.sca(ax[1])
plt.scatter(df_chan['T(K)'], dev, 3, marker='s', facecolors="None",edgecolors=ccycle[1],linewidths=.5,label='Chan & Giauque 1964')

dev = calc_dev(df_hildenbrand['T(K)'],df_hildenbrand['massfrac'],df_hildenbrand['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_hildenbrand['massfrac']*100, dev, 3, marker='d', facecolors="None",edgecolors=ccycle[2],linewidths=.5,label='Hildenbrand & Giauque 1953')
plt.sca(ax[1])
plt.scatter(df_hildenbrand['T(K)'], dev, 3, marker='d', facecolors="None",edgecolors=ccycle[2],linewidths=.5,label='Hildenbrand & Giauque 1953')

dev = calc_dev(df_wrewsky['T(K)'],df_wrewsky['massfrac'],df_wrewsky['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_wrewsky['massfrac']*100, dev, 3, marker='P', facecolors="None",edgecolors=ccycle[3],linewidths=.5,label='Wrewsky & Kaigorodoff 1924')
plt.sca(ax[1])
plt.scatter(df_wrewsky['T(K)'], dev, 3, marker='P', facecolors="None",edgecolors=ccycle[3],linewidths=.5,label='Wrewsky & Kaigorodoff 1924')

dev = calc_dev(df_chernenkaya['T(K)'],df_chernenkaya['massfrac'],df_chernenkaya['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_chernenkaya['massfrac']*100, dev, 3, marker='X', facecolors="None",edgecolors=ccycle[4],linewidths=.5,label='Chernenkaya 1971')
plt.sca(ax[1])
plt.scatter(df_chernenkaya['T(K)'], dev, 3, marker='X', facecolors="None",edgecolors=ccycle[4],linewidths=.5,label='Chernenkaya 1971')

dev = calc_dev(df_allred['T(K)'],df_allred['massfrac'],df_allred['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_allred['massfrac']*100, dev, 3, marker='^', facecolors="None",edgecolors=ccycle[5],linewidths=.5,label='Allred & Wolley 1981')
plt.sca(ax[1])
plt.scatter(df_allred['T(K)'], dev, 3, marker='^', facecolors="None",edgecolors=ccycle[5],linewidths=.5,label='Allred & Wolley 1981')

dev = calc_dev(df_fujita['T(K)'],df_fujita['massfrac'],df_fujita['cp(J/kgK)'])
plt.sca(ax[0])
plt.scatter(df_fujita['massfrac']*100, dev, 3, marker='v', facecolors="None",edgecolors=ccycle[6],linewidths=.5,label='Fujita et al. 2008')
plt.sca(ax[1])
plt.scatter(df_fujita['T(K)'], dev, 3, marker='v', facecolors="None",edgecolors=ccycle[6],linewidths=.5,label='Fujita et al. 2008')

ax[0].text(0.025, 0.98, 'a)', horizontalalignment='left', verticalalignment='top', transform=ax[0].transAxes)
ax[1].text(0.025, 0.98, 'b)', horizontalalignment='left', verticalalignment='top', transform=ax[1].transAxes)

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

for idx,wt in enumerate(wt_ls):
    plt.scatter(1,1,3, marker='o',facecolors="None",edgecolors=colour_ls[idx],linewidths=.5,label=f'{wt} wt%',alpha=0)
plt.scatter(0,250,alpha=0)
plt.scatter(1,1, 3, marker='s', facecolors="None",edgecolors=ccycle[1],linewidths=.5,label='CG64', alpha=0)
plt.scatter(1,1, 3, marker='d', facecolors="None",edgecolors=ccycle[2],linewidths=.5,label='HG53', alpha=0)
plt.scatter(1,1, 3, marker='P', facecolors="None",edgecolors=ccycle[3],linewidths=.5,label='WK24', alpha=0)
plt.scatter(1,1, 3, marker='X', facecolors="None",edgecolors=ccycle[4],linewidths=.5,label='C71', alpha=0)
plt.scatter(1,1, 3, marker='^', facecolors="None",edgecolors=ccycle[5],linewidths=.5,label='AW81', alpha=0)
plt.scatter(1,1, 3, marker='v', facecolors="None",edgecolors=ccycle[6],linewidths=.5,label='F08', alpha=0)

leg = plt.legend(bbox_to_anchor=(1.04, .5),loc='center left',prop={'size': 4},
                 frameon=True, edgecolor='k', framealpha=1)
leg.get_frame().set_linewidth(0.5)
for lh in leg.legendHandles:
    lh.set_alpha(1)

plt.tight_layout()

# base.show_plot_max()
plt.savefig(fd + 'z_deviations_all_small.png')

plt.close('all')
## debug end line
print('Finished running script')

