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
colour_ls = ['#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84','#081d58']

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
df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
df_St = pd.read_csv(ld + 'Steve_cp.csv', header=0)

## PLOT FIG 9 OF TRF 1998

fig, ax = plt.subplots()
# plot experimental data
for idx,wt in enumerate(wt_ls):
    # idx=6
    # wt=wt_ls[idx]

    colour = colour_ls[idx]
    m = mass_ls[idx]
    data = base.DataCP(m,wt)
    data.import_data(dd, mean=0)
    data.calc_shomate('all')
    TRF_interp = interp1d(df_TRF['T(K)'], df_TRF[str(wt)])
    St_interp = interp1d(df_St['T (K)'], df_St[str(wt)])
    # SF_interp = interp1d(df_SF['T(K)'], df_SF[str(wt)])

    # plt.scatter(data.df_p['T(K)'], data.df_p['cp(J/gK)'], 1, c=colour, marker='o',label=f'{wt} %')
    # plt.plot(data.df_p['T(K)'], base.shomate_eqn(data.df_p['T(K)'], *data.shomate_p), c=colour,label=f'{wt} %')

    plt.plot(data.df_p['T(K)'],TRF_interp(data.df_p['T(K)']),c=colour,linestyle=(0, (5, 1)))
    plt.plot(data.df_p['T(K)'],St_interp(data.df_p['T(K)']),c=colour)
    # plt.plot(df_SF['T(K)'],df_SF[str(wt)],c=colour,linestyle=(0, (5, 1)))
    # plt.title('Cp wrt SF')
    plt.xlabel('Temperature (K)');
    plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')

    # plt.savefig(fd + '{}wt%_Cp_TRF.png'.format(wt))


# plt.scatter(CG['33']['T(K)'], CG['33']['cp(J/gK)'], 8, marker='o', facecolors='k', edgecolors='k', linewidths=.1,
#             label='CG 33 wt%')
# plt.ylim([4,4.6])
# plt.xlim([260,320])
plt.legend(prop={'size': 5})
# plt.legend(bbox_to_anchor=(1, .464))
base.show_plot_max()
plt.savefig(fd + 'all_Cp_SF_shomate.png')

## debug end line
print('Finished running script')

