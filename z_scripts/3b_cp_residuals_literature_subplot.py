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

def calc_dev(wt,df):
    if TRF_arg == 0:
        df_dev = df_TRF
    elif TRF_arg == 1:
        df_dev = df_St
    TRF_interp = interp1d(df_dev['T(K)'], df_dev[str(wt)])
    TRF_val = TRF_interp(df['T(K)'])
    dev = 100 * (1-TRF_val/df['cp(J/gK)'])
    return(dev)

plt.close('all')

fig, ax = plt.subplots(2,1,sharey=True)
base.splt_axis_label(fig,"",'$100(1-C_p{}^{EoS}/C_p{}^{expt})$ (%)')
plt.tight_layout()
plt.subplots_adjust(hspace=0.5)
plt.sca(ax[0])
plt.ylim(-20, 20)
plt.axhline(0, linewidth=.3, color='k')
plt.xlim(0,100)
ax[0].set_xlabel('$w_{NH_3}$ (wt%)')
plt.sca(ax[1])
plt.ylim(-20, 20)
plt.axhline(0, linewidth=.3, color='k')
plt.xlim(180,340)
plt.xlabel('T (K)')

# plot experimental data
for idx,wt in enumerate(wt_ls):
    # idx=6
    # wt=26.912

    colour = colour_ls[idx]
    m = mass_ls[idx]
    data = base.DataCP(m,wt)
    data.import_data(dd, mean=0)
    dev = calc_dev(wt, data.df_p)
    plt.sca(ax[0])
    plt.scatter([float(wt)] * len(dev), dev, 2, marker='o', facecolors=colour, linewidths=.3,label=f'{wt} wt%')
    plt.sca(ax[1])
    plt.scatter(data.df_p['T(K)'], dev, 2, marker='o', facecolors=colour, linewidths=.3,label=f'{wt} wt%')

# plot chan giauque 33 wt %
dev = calc_dev('33',CG['33'])
plt.sca(ax[0])
plt.scatter([33] * len(dev), dev, 8, marker='o', facecolors='k', edgecolors='k', linewidths=.1,label='Chan & Giauque')
plt.sca(ax[1])
plt.scatter(CG['33']['T(K)'], dev, marker='o', facecolors='k', edgecolors='k', linewidths=.1,label='Chan & Giauque')

# plot hildenbrand giauque other wt %
for _, wt in enumerate(HG):
    dev = calc_dev(wt,HG[wt])
    plt.sca(ax[0])
    plt.scatter([float(wt)]*len(dev),dev,8,marker='o',facecolors='w', edgecolors='k',linewidths=.1)
    plt.sca(ax[1])
    plt.scatter(HG[wt]['T(K)'], dev, 8, marker='o', facecolors='w', edgecolors='k',linewidths=.1)

plt.sca(ax[0])
plt.scatter(-1,0, marker='o', facecolors='w', edgecolors='k',linewidths=.1,label='Hildenbrand & Giauque')
plt.sca(ax[1])
plt.scatter(-1,0, marker='o', facecolors='w', edgecolors='k',linewidths=.1,label='Hildenbrand & Giauque')

# plot wrewsky and kaigorodoff
plt.sca(ax[0])
plt.scatter(df_WK['wt'], df_WK['dev'], 4, marker='x', facecolors='k', linewidths=.3,label='Wrewsky & Kaigorodoff')
plt.sca(ax[1])
plt.scatter(df_WK['T(K)'], df_WK['dev'], 4, marker='x', facecolors='k', linewidths=.3,label='Wrewsky & Kaigorodoff')

ax[0].legend(markerscale=10,fontsize=5)
ax[0].legend(bbox_to_anchor=(1, 1),loc='upper left')
base.show_plot_max()
plt.savefig(fd + 'deviations_all.png')

plt.close('all')
## debug end line
print('Finished running script')

