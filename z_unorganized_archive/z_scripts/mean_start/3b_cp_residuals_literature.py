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
from scipy.interpolate import UnivariateSpline
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


ramp_rate = 0.1
wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912]
mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]# based off liquidus alignment - 0.5K
colour_ls = ['tab:red','tab:orange','tab:olive','tab:green','tab:blue','tab:purple','tab:pink']

## LOAD LITERATURE DATA

# df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
df_TRF = pd.read_csv(ld + 'Tillner-Roth_Friend_liq.csv', header=0)
CG = {'33': pd.read_csv(ld + 'Chan_Giauque_0.33.csv', header=0)}
HG = { '49': pd.read_csv(ld + 'Hildenbrand_Giauque_0.49.csv', header=0),
           '59': pd.read_csv(ld + 'Hildenbrand_Giauque_0.59.csv', header=0),
           '61': pd.read_csv(ld + 'Hildenbrand_Giauque_0.61.csv', header=0),
           '64': pd.read_csv(ld + 'Hildenbrand_Giauque_0.64.csv', header=0),
           '65': pd.read_csv(ld + 'Hildenbrand_Giauque_0.65.csv', header=0)}
df_WK = pd.read_csv(ld + 'Wrewsky_Kaigorodoff.csv', header=0)

## PLOT FIG 9 OF TRF 1998

def calc_dev(wt,df):
    if wt == 4.91:
        TRF_interp = interp1d(df_TRF['T(K)'], df_TRF[str(4.9)])
    else:
        TRF_interp = interp1d(df_TRF['T(K)'], df_TRF[str(wt)])
    TRF_val = TRF_interp(df['T(K)'])
    dev = 100 * (1-TRF_val/df['cp(J/gK)'])
    return(dev)

def create_fig(name):
    plt.figure(name)
    plt.ylim(-20, 20)
    plt.ylabel('$100(1-C_p{}^{calc}/C_p{}^{exp})$')
    plt.axhline(0, linewidth=.3, color='k')

plt.close('all')
create_fig('dev_wt')
create_fig('dev_T')

# plot experimental data
for idx,wt in enumerate(wt_ls):
    # idx=5
    # wt=20.1

    colour = colour_ls[idx]
    m = mass_ls[idx]
    data = base.DataCP(m,wt)
    data.import_data(dd, mean=0)
    dev = calc_dev(wt, data.df_p)
    plt.figure('dev_wt')
    plt.scatter([float(wt) / 100] * len(dev), dev, 2, marker='o', facecolors=colour, linewidths=.3,label=f'{wt} wt%')
    plt.figure('dev_T')
    plt.scatter(data.df_p['T(K)'], dev, 2, marker='o', facecolors=colour, linewidths=.3,label=f'{wt} wt%')

# plot chan giauque 33 wt %
dev = calc_dev('33',CG['33'])
plt.figure('dev_wt')
plt.scatter([33 / 100] * len(dev), dev, 8, marker='o', facecolors='k', edgecolors='k', linewidths=.1,label='Chan & Giauque')
plt.figure('dev_T')
plt.scatter(CG['33']['T(K)'], dev, marker='o', facecolors='k', edgecolors='k', linewidths=.1,label='Chan & Giauque')

# plot hildenbrand giauque other wt %
for _, wt in enumerate(HG):
    dev = calc_dev(wt,HG[wt])
    plt.figure('dev_wt')
    plt.scatter([float(wt)/100]*len(dev),dev,8,marker='o',facecolors='w', edgecolors='k',linewidths=.1)
    plt.figure('dev_T')
    plt.scatter(HG[wt]['T(K)'], dev, 8, marker='o', facecolors='w', edgecolors='k',linewidths=.1)

plt.figure('dev_wt')
plt.scatter(-1,0, marker='o', facecolors='w', edgecolors='k',linewidths=.1,label='Hildenbrand & Giauque')
plt.figure('dev_T')
plt.scatter(-1,0, marker='o', facecolors='w', edgecolors='k',linewidths=.1,label='Hildenbrand & Giauque')

# plot wrewsky and kaigorodoff
plt.figure('dev_wt')
plt.scatter(df_WK['wt']/100, df_WK['dev'], 4, marker='x', facecolors='k', linewidths=.3,label='Wrewsky & Kaigorodoff')
plt.figure('dev_T')
plt.scatter(df_WK['T(K)'], df_WK['dev'], 4, marker='x', facecolors='k', linewidths=.3,label='Wrewsky & Kaigorodoff')


plt.figure('dev_wt')
plt.xlim(0,1)
plt.xlabel('$w_{NH3}$')
plt.legend(markerscale=10,fontsize=5)
plt.legend(bbox_to_anchor=(1,1))
plt.savefig(fd + 'deviations_wt.png')
# base.show_plot_max()

plt.figure('dev_T')
plt.xlim(180,340)
plt.xlabel('$T(K)$')
plt.legend(markerscale=10,fontsize=5)
plt.legend(bbox_to_anchor=(1,1))
plt.savefig(fd + 'deviations_T.png')
# base.show_plot_max()

plt.close('all')
## debug end line
print('Finished running script')

