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
comp_arg = 1
melt_arg = 0
pure_arg = 1
best_fit_arg = 1

ramp_rate = 0.1
wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K

mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

colour_ls = ['#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84','#081d58']

## PLOT ALL CP ON SAME PLOT. OPTIONS FOR SF AND TRF EOS
plt.close()
fig, ax = plt.subplots(1,2)
# fig = plt.figure('combined')
# ax = fig.add_subplot(111)

coef_num =5
coef_arr = np.empty((len(wt_ls),coef_num))
sd_arr = np.empty((len(wt_ls),coef_num))

for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]
    m = mass_ls[idx]

    data = base.DataCP(m,wt)
    data.import_data(dd,mean=0)

    data.calc_shomate('all')

    # save coef and sd of polynomial for best fit calculation of all fits
    coef_arr[idx] = data.shomate_m
    sd_arr[idx] = data.sd_m

    # plt.figure('combined')
    plt.sca(ax[0])
    data.plot_data(colour, idx, raw_data=1, shomate=1,errors=1,melt=1, pure=0)
    plt.sca(ax[1])
    data.plot_data(colour, idx, raw_data=1, shomate=1,errors=1,melt=0, pure=1)

    print('ye')
if best_fit_arg==1 and pure_arg==0:
    # calculate overall best fit
    weights_arr = (1 / sd_arr ** 2)
    coef_mean = np.sum(coef_arr * weights_arr, 0) / np.sum(weights_arr, 0)
    sd_mean = np.sqrt(1 / np.sum(weights_arr, 0))
    # plot the overall best fit
    T_range = np.arange(183.5, 234.5)
    plt.plot(T_range, base.shomate_eqn(T_range,*coef_mean), 'k',label = 'overall fit')

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
elif melt_arg == 1 and pure_arg != 1:
    figname = 'melt'
    title = 'Partial Melting Regime'
    plt.ylim([2.5, 4.5])
elif melt_arg != 1 and pure_arg == 1:
    figname = 'pure'
    title = 'Liquid Phase'
    # plt.ylim([3,5])
    # plt.ylim([4.2,4.6])
    # plt.xlim([260,320])

# plot the rest
plt.xlabel('Temperature (K)');
plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
plt.legend(prop={'size': 5})
plt.legend(bbox_to_anchor=(1,1))



plt.title(title)

if best_fit_arg == 1 and pure_arg==0:
    figBF = '_BF'
else:
    figBF = ''
plt.savefig(fd + '{}{}_zoom_test.png'.format(figname,figBF))
# base.show_plot_max()
# plt.close()
print('Finished for loop {} / {}'.format(idx,len(wt_ls)-1))



## debug end line
print('Finished running script')
