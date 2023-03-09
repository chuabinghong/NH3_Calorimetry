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
melt_arg = 1
pure_arg = 0
best_fit_arg = 0

ramp_rate = 0.1
# wt_ls = [1.42,2,2.5,2.9,4.9,7.9,8.1,9.7,14.2,20,26.9] # based off liquidus alignment
wt_ls = [1.42,2.5,4.9,8.1,9.7,26.9,2.0,2.9,4.91,7.9,14.2,20] # based off liquidus alignment
# wt_ls = [1.42,2.84,5.11,6.88,8.73,25.26] # based off lab measurements 26C
# wt_ls = [2.00, 4.03, 9.45, 16.00, 21.01]  # based off lab measurements 46C
# wt_ls = [1.42,2.84,5.11,6.88,8.73,25.26, 2.00, 4.03, 9.45, 16.00, 21.01]

colour_ls = ['tab:blue','yellow','tab:green','tab:cyan','tab:orange','tab:brown','tab:red','tab:purple','lightpink','tab:pink','tab:olive','dimgray']
# colour_ls = ['tab:green','tab:orange','tab:brown','tab:red','tab:purple','tab:pink','tab:olive','dimgray']


# adjust script arguments according to best_fit
if best_fit_arg == 1:
    wt_ls = [7.9, 8.1, 9.7, 14.2, 20, 26.9]  # ONLY GOOD VIBE DATA
    colour_ls = ['tab:brown', 'tab:red', 'tab:purple', 'tab:pink', 'tab:olive', 'dimgray']  # ONLY GOOD VIBE COLOURS

## LOAD STEVE'S EOS

df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
df_TRF = pd.read_csv(ld + 'Tillner-Roth_Friend_liq.csv', header=0)

## PLOT ALL CP ON SAME PLOT. OPTIONS FOR SF AND TRF EOS
plt.close()
fig = plt.figure('combined')
ax = fig.add_subplot(111)

# to plot SF or TRF
for idx, wt in enumerate(wt_ls):
    plt.figure('combined')
    colour = colour_ls[idx]
    # plt.plot(df_TRF['T(K)'], df_TRF['{}'.format(wt)], '-.', color=colour)
    # plt.plot(df_SF['T(K)'],df_SF['{}'.format(wt)]/1000,'--',color=colour)

if fit_type == 'linear':
    coef_num = 2
elif fit_type == 'shomate':
    coef_num =5
coef_arr = np.empty((len(wt_ls),coef_num))
sd_arr = np.empty((len(wt_ls),coef_num))

for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]

    data = base.DataCP(wt)

    # data.import_data_man(dd)
    data.import_data(dd)

    plt.figure('combined')

    if fit_type == 'linear':
        data.calc_linear()
        # save coef and sd of polynomial for best fit calculation of all fits
        coef_arr[idx] = data.coef_m
        sd_arr[idx] = data.sd_m
        data.plot_fit(colour, melt=melt_arg, pure=pure_arg)

    elif fit_type == 'spline':
        data.calc_spline()
        data.plot_spl(colour, melt=melt_arg, pure=pure_arg)
    # data.plot_fit(colour,melt=1)

    elif fit_type =='shomate':
        data.calc_shomate('all')
        # save coef and sd of polynomial for best fit calculation of all fits
        coef_arr[idx] = data.shomate_m
        sd_arr[idx] = data.sd_m
        data.plot_fit(colour, melt=melt_arg, pure=pure_arg)

if fit_type == 'linear' or 'shomate' and best_fit_arg==1:
    # calculate overall best fit
    weights_arr = (1 / sd_arr ** 2)
    coef_mean = np.sum(coef_arr * weights_arr, 0) / np.sum(weights_arr, 0)
    sd_mean = np.sqrt(1 / np.sum(weights_arr, 0))
    # plot the overall best fit
    T_range = np.arange(182, 260)
    plt.plot(T_range, base.shomate_eqn(T_range,*coef_mean), 'k',label = 'overall fit')


# plot the rest
plt.xlabel('Temperature (K)');
plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
plt.title('Cp {} Fit'.format(fit_type.capitalize()))
plt.title('Cp'.format(fit_type.capitalize()))
plt.ylim([2,7])
# plt.ylim([4,4.5])
# plt.xlim([176,320])
# plt.xlim([250,320])
# plt.plot(df_SF['T(K)'][1],df_SF['{}'.format(wt)][1],'-.',color='k',label='Tillner-Roth Friend EOS')
# plt.plot(df_SF['T(K)'][1],df_SF['{}'.format(wt)][1],'--',color='k',label='SeaFreeze EOS')

if comp_arg == 1:
    # plot secondary axis of composition
    plt.xlim([176,265])
    Xl_calc_v = np.vectorize(base.Xl_calc)
    plt.tick_params(axis='x', which='both', top=False)
    secax = ax.secondary_xaxis('top', functions=(Xl_calc_v, base.T_calc))
    secax.set_xlabel('Mass Fraction (wt%)')


plt.legend(prop={'size': 5})
if melt_arg == 1 and pure_arg == 1:
    figname = 'all'
elif melt_arg == 1 and pure_arg != 1:
    figname = 'melt'
elif melt_arg != 1 and pure_arg == 1:
    figname = 'pure'
if best_fit_arg == 1:
    figBF = '_BF'
elif best_fit_arg == 0:
    figBF = ''
# plt.savefig(fd + '{}_{}{}.png'.format(figname,fit_type,figBF))
base.show_plot_max()

print('Finished for loop {} / {}'.format(idx,len(wt_ls)-1))



## debug end line
print('Finished running script')
