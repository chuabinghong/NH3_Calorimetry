"""
created by Bing Hong CHUA 11Nov22

script objective:

"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
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
fd = r"o_specificHeat/"

## SCRIPT ARGUMENTS
ice_arg = 1 # 1 if using ice cp from our own 0wt% experiments
range_arg = 0 #1 if using -196C to 46C
wt_arg = 2
plot_arg = 1 # 1 if want to plot
T_bin = 3 # temperature bin width (in K) used to smooth

mass_ls = [4.6293,4.5858,4.5202,3.7778]
ramp_rate = 0.1
if wt_arg == 0:
    wt_ls = [1.42, 2.84, 5.11, 6.88, 8.73, 25.26]  # based off lab measurements
elif wt_arg ==1:
    wt_ls = [1.42, 2.5, 4.9, 8.1, 9.7, 26.9]  # based off liquidus alignment
elif wt_arg ==2:
    wt_ls = [ 5.2, 8.4, 10.0, 26.93]  # based off liquidus alignment
# wt_ls = [1.42,2.84,5.11,6.88,8.73,25.26] # based off lab measurements

# adjust script arguments according to range_arg
if range_arg == 1:
    dd = r"i_data_processed/"
    od = r"i_data_processed/"
    fd = r"o_specificHeat/"
    mass_ls = [4.5386,4.1943,3.8153,3.7107]
    if wt_arg == 0:
        wt_ls = [2.00,4.03,9.45,16.00,21.01] # based off lab measurements
    elif wt_arg == 1:
        wt_ls = [2.0, 2.9, 4.91, 7.9, 14.2, 20]  # based off liquidus alignment
    elif wt_arg == 2:
        wt_ls = [5.2, 8.2, 14.4, 20.1]  # based off liquidus alignment

## DATA PROCESSING

for idx, wt in enumerate(wt_ls):
    # idx = 1 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging

    m = mass_ls[idx]

    data = base.DataCP(m,wt)
    # data.import_data_man(dd)
    data.import_data(dd,mean=0)


    # smooth Cp average
    data.mean_cp('melt',T_bin,od)
    data.mean_cp('pure',T_bin,od)

    # Cp moving average (based on T_bin)
    data.calc_SMA('melt',1,od)
    data.calc_SMA('pure',1,od)

## plot stuff

    if plot_arg == 1:
        coeff_dict = {}
        data.calc_shomate('mean')
        coeff_dict['m_mean'] = data.shomate_m
        coeff_dict['p_mean'] = data.shomate_p
        data.calc_shomate('SMA')
        coeff_dict['m_SMA'] = data.shomate_m
        coeff_dict['p_SMA'] = data.shomate_p
        data.calc_shomate('all')
        coeff_dict['m_all'] = data.shomate_m
        coeff_dict['p_all'] = data.shomate_p

        plt.close('all')
        plt.figure()
        plt.scatter(data.df_m['T(K)'], data.df_m['cp(J/gK)'], .1, c = 'tab:blue', marker='o',label='data',alpha=.5)
        plt.scatter(data.df_m['T(K)'], data.df_m['cp_SMA_1K'], 0.1, c = 'tab:orange', marker='^',label='moving mean',alpha=.5)
        plt.scatter(data.df_m_mean['T(K)'], data.df_m_mean['cp(J/gK)'], 5, c = 'k', marker='x',label=f'{T_bin} K average',alpha=1)
        plt.errorbar(data.df_m_mean['T(K)'], data.df_m_mean['cp(J/gK)'], yerr=data.df_m_mean['std'], linestyle='',
                     ecolor='k', capsize=1, capthick=0.5, elinewidth=0.5,errorevery=20)

        plt.plot(data.df_m['T(K)'], base.shomate_eqn(data.df_m['T(K)'], *coeff_dict['m_all']),c='blue', label='data shomate')
        plt.plot(data.df_m['T(K)'], base.shomate_eqn(data.df_m['T(K)'], *coeff_dict['m_SMA']),c='orange', label='moving mean shomate')
        plt.plot(data.df_m['T(K)'], base.shomate_eqn(data.df_m['T(K)'], *coeff_dict['m_mean']),c='k', label=f'{T_bin} K mean shomate')

        plt.xlabel('Temperature (K)');
        plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
        plt.title(f'Cp melt {wt} wt% Fit')
        plt.legend(fontsize=3)

        plt.savefig(fd + f'{wt}wt%_cp_melt_fit_{data.m}g.png')
        # base.show_plot_max()


        plt.close('all')
        plt.figure()
        plt.scatter(data.df_p['T(K)'], data.df_p['cp(J/gK)'], .1, c = 'tab:blue', marker='o',label='data',alpha=.5)
        plt.scatter(data.df_p['T(K)'], data.df_p['cp_SMA_1K'], 0.1, c = 'tab:orange', marker='^',label='moving mean',alpha=.5)
        plt.scatter(data.df_p_mean['T(K)'], data.df_p_mean['cp(J/gK)'], 5, c = 'k', marker='x',label='1 K average')
        plt.errorbar(data.df_p_mean['T(K)'], data.df_p_mean['cp(J/gK)'],yerr=data.df_p_mean['std'],linestyle='',ecolor='k',capsize=1,capthick=0.5,elinewidth=0.5,errorevery=20)

        plt.plot(data.df_p['T(K)'], base.shomate_eqn(data.df_p['T(K)'], *coeff_dict['p_all']),c='blue', label='data shomate')
        plt.plot(data.df_p['T(K)'], base.shomate_eqn(data.df_p['T(K)'], *coeff_dict['p_SMA']),c='orange', label='moving mean shomate')
        plt.plot(data.df_p['T(K)'], base.shomate_eqn(data.df_p['T(K)'], *coeff_dict['p_mean']),c='k', label='1 K mean shomate')

        plt.xlabel('Temperature (K)');
        plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
        plt.title(f'Cp liquid {wt} wt% Fit')
        plt.legend(fontsize=3)

        plt.savefig(fd + f'{wt}wt%_cp_pure_fit_{data.m}g.png')
    # base.show_plot_max()


    print('Finished for loop {} / {}'.format(idx+1,len(wt_ls)))

## debug end line
print('Finished running script')

## EXTRAS
# import numpy as np
#
# # Deviation plot for ice Cp
# T = np.linspace(240,273,100)
#
# ice_data = base.DataPrep(2.5108, 0)
# ice_data.import_data(dd)
# ice_data.calc_cp_pure()
# iceCp_EXP = UnivariateSpline(ice_data.df_cp_p['T(K)'], ice_data.df_cp_p['cp(J/gK)'], k=3, ext=0)
# data_EXP = iceCp_EXP(T)
#
# df_iceCp = pd.read_csv(ld + 'ice_Cp.csv', names=['T', 'Cp'])
# iceCp_FS = interp1d(df_iceCp['T'], df_iceCp['Cp'] / 1000)
# data_FS = iceCp_FS(T)
#
# dev = data_EXP - data_FS
# dev_data = ice_data.df_cp_p['cp(J/gK)'].iloc[6489:-1] - iceCp_FS(ice_data.df_cp_p['T(K)'].iloc[6489:-1])
# dev_data_pct = dev_data / iceCp_FS(ice_data.df_cp_p['T(K)'].iloc[6489:-1])
# plt.figure()
# # plt.plot(dev)
# plt.plot(dev_data_pct)
# plt.axhline(0,'k--')


# # plot Lh
# fig = plt.figure()
# plt.plot(np.arange(180,300),iceLh_f(np.arange(180,300)),label='Extrapolation')
# plt.plot(df_iceLh['T'],df_iceLh['Lh'],linewidth=1.5,label='SeaFreeze')
# plt.plot(273.15,333.55,'x',color='k')
# plt.xlabel('T (K)')
# plt.ylabel(r'$\Delta$H (J $g^{-1}$)')
# plt.title('Latent Heat of Melting (Pure Ice)')
# plt.legend()
# # plt.show()
# saveFig = fd + 'Lh_SF.png'
# plt.savefig(saveFig)
# plt.close()