"""
created by Bing Hong CHUA 29Sep22

script objective:
from cut specific heat data, apply various spline and moving average fits for further analysis
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import base

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
comp_arg = 0
melt_arg = 0
pure_arg = 1
lit_arg = 0
best_fit_arg = 0
water_arg = 1
inset_arg = 1

ramp_rate = 0.1
wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K

mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

# colour_ls = ['#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84','#081d58']
# colour_ls = ['#d73027','#fc8d59','#fee090','#f2f26d','#e0f3f8','#91bfdb','#4575b4']
colour_ls = ['#d73027','#f46d43','#fdae61','#e8e884','#abd9e9','#74add1','#4575b4']
marker_ls = ['o','s','d','P','X','^','v']
## PLOT ALL CP ON SAME PLOT. OPTIONS FOR SF AND TRF EOS
plt.close()
fig = plt.figure('combined')
ax = fig.add_subplot(111)

if water_arg ==1:
    df_water = pd.read_csv(ld + 'water_Cp.csv',names=['T(K)','Cp'])
    plt.plot(df_water['T(K)'],df_water['Cp'],linewidth=0.75,c='k')
    plt.text(0.85, 0.58, '$H_{2}O$ (0 wt%)', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)


coef_num =5
coef_arr = np.empty((len(wt_ls),coef_num))
sd_arr = np.empty((len(wt_ls),coef_num))

for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]
    marker = marker_ls[idx]
    m = mass_ls[idx]

    data = base.DataCP(m,wt)
    data.import_data(dd,mean=0)

    data.calc_shomate('all')

    # save coef and sd of polynomial for best fit calculation of all fits
    coef_arr[idx] = data.shomate_m
    sd_arr[idx] = data.sd_m

    if lit_arg == 1:
        from scipy.interpolate import interp1d
        df_St = pd.read_csv(ld + 'Steve_cp.csv', header=0)
        TRF_interp = interp1d(df_St['T(K)'], df_St[str(wt)])
        plt.plot(data.df_p['T(K)'],TRF_interp(data.df_p['T(K)']),'--',c=colour,zorder=2,linewidth=0.9)
        # plt.plot(np.arange(200,300,1),TRF_interp(np.arange(200,300,1)),'--',c=colour,zorder=2,linewidth=0.9)

    plt.figure('combined')
    data.plot_data(colour, marker, idx, raw_data=1, shomate=1,errors=1,melt=melt_arg, pure=pure_arg)

if best_fit_arg==1:
    # calculate overall best fit
    weights_arr = (1 / sd_arr ** 2)
    coef_mean = np.sum(coef_arr * weights_arr, 0) / np.sum(weights_arr, 0)
    sd_mean = np.sqrt(1 / np.sum(weights_arr, 0))
    # plot the overall best fit
    T_range = np.arange(184.5, 234.5)
    plt.plot(T_range, base.shomate_eqn(T_range,*coef_mean),'k',linewidth=1, zorder=5,label = 'overall fit')

if lit_arg==1:
    plt.plot(0, 0, 'k--', label='EoS')

if comp_arg == 1 and pure_arg==0:
    # plot secondary axis of composition
    plt.xlim([180,240])
    Xl_calc_v = np.vectorize(base.Xl_calc)
    plt.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(Xl_calc_v, base.T_calc))
    secax.set_xlabel('Mass Fraction (wt%)')

if melt_arg == 1 and pure_arg == 1:
    figname = 'all'
    title = ''
    # plt.ylim([3.45, 4.65])
    plt.xlim([180, 315])
elif melt_arg == 1 and pure_arg != 1:
    figname = 'melt'
    title = 'Partial Melting Regime'
    # plt.ylim([2.5, 4.5])
elif melt_arg != 1 and pure_arg == 1:
    figname = 'pure'
    title = 'Liquid Phase'
    # plt.ylim([4.2,4.5])
    # plt.xlim([255,312])

    plt.ylim([3.45,4.65])
    plt.xlim([210,315])

# plot the rest
plt.xlabel('Temperature (K)');
plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
plt.legend(prop={'size': 5})
# plt.legend(bbox_to_anchor=(1,1))

# plt.title(title)

if best_fit_arg == 1 and pure_arg==0:
    figBF = '_BF'
else:
    figBF = ''

# base.show_plot_max()

plt.savefig(fd + '{}{}.png'.format(figname,figBF))
# plt.close()
print('Finished for loop {} / {}'.format(idx,len(wt_ls)-1))



## debug end line
print('Finished running script')
