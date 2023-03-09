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
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
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
fit_type = 'shomate'
melt_arg = 0
pure_arg = 1
best_fit_arg = 0

ramp_rate = 0.1
# wt_ls = [1.42,2,2.5,2.9,4.9,7.9,8.1,9.7,14.2,20,26.9] # based off liquidus alignment
wt_ls = [1.42,2.5,4.9,8.1,9.7,26.9,2.0,2.9,7.9,14.2,20] # based off liquidus alignment
# wt_ls = [1.42,2.84,5.11,6.88,8.73,25.26] # based off lab measurements 26C
# wt_ls = [2.00, 4.03, 9.45, 16.00, 21.01]  # based off lab measurements 46C
# wt_ls = [1.42,2.84,5.11,6.88,8.73,25.26, 2.00, 4.03, 9.45, 16.00, 21.01]

colour_ls = ['tab:blue','yellow','tab:green','tab:cyan','tab:orange','tab:brown','tab:red','tab:purple','tab:pink','tab:olive','dimgray']
# colour_ls = ['tab:green','tab:orange','tab:brown','tab:red','tab:purple','tab:pink','tab:olive','dimgray']

## PLOT ALL CP ON SAME PLOT. OPTIONS FOR SF AND TRF EOS
# plt.close()
fig = plt.figure('combined')
ax = fig.add_subplot(111)

for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]
    data = base.DataCP(wt)

    # data.import_data_man(dd)
    data.import_data(dd,'savgol')

    plt.figure('combined')
    if melt_arg == 1 and pure_arg != 1:
        phase = 'melt'
        plt.plot(data.df_m['T(K)'],data.df_m['cp_err(%)'],linewidth=0.35,color=colour,label=f'{wt} wt%')
    if melt_arg != 1 and pure_arg == 1:
        phase = 'pure'
        plt.plot(data.df_p['T(K)'], data.df_p['cp_err(%)'],linewidth=0.35,color=colour, label=f'{wt} wt%')
    # data.calc_shomate(base.shomate_eqn)
    # # save coef and sd of polynomial for best fit calculation of all fits
    # coef_arr[idx] = data.shomate_m
    # sd_arr[idx] = data.sd_m
    # data.plot_fit(colour, melt=melt_arg, pure=pure_arg)


# plot the rest
plt.xlabel('Temperature (K)');
plt.ylabel('Specific Heat Error (%)')
plt.title(f'Cp {phase} Errors')

if comp_arg == 1:
    # plot secondary axis of composition
    plt.xlim([176,265])
    Xl_calc_v = np.vectorize(base.Xl_calc)
    plt.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(Xl_calc_v, base.T_calc))
    secax.set_xlabel('Mass Fraction (wt%)')


plt.legend(prop={'size': 5})
plt.savefig(fd + f'err_{phase}.png')
# base.show_plot_max()

print('Finished for loop {} / {}'.format(idx,len(wt_ls)-1))



## debug end line
print('Finished running script')
