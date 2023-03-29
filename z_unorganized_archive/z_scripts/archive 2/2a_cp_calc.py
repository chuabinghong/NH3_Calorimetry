"""
created by Bing Hong CHUA 28Sepl22

script objective:
calculate specific heat
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
range_arg = 1 #1 if using -196C to 46C
wt_arg = 2
T_bin = 1 # temperature bin width (in K) used to smooth

ramp_rate = 0.1
mass_ls = [4.6293,4.5858,4.5202,3.7778]
if wt_arg == 0:
    wt_ls = [1.42, 2.84, 5.11, 6.88, 8.73, 25.26]  # based off lab measurements
elif wt_arg ==1:
    wt_ls = [1.42, 2.5, 4.9, 8.1, 9.7, 26.9]  # based off liquidus alignment
elif wt_arg ==2:
    wt_ls = [5.2, 8.4, 10.0, 26.93]  # based off liquidus alignment


# m_lb_ls = [181.5, 181.5, 181.9, 181.5, 181.5, 184]
# m_ub_ls = [263, 262, 240, 240, 239, 235]
# p_lb_ls = [279.5,276.5,273.5,268,265,214] #OLD ONE
# p_lb_ls = [283.01,279.01,276.01,268,265,214]
# p_ub_ls = [310,310,310,310,310,310]

m_lb_ls = [181.9, 181.5, 181.5, 184]
m_ub_ls = [240, 240, 239, 235]
p_lb_ls = [276.01,268,265,214]
p_ub_ls = [310,310,310,310]

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

    # m_lb_ls = [181.5, 181.7, 181.5, 181.5, 182, 183]
    # m_ub_ls = [263, 263, 261.5, 260, 250, 230]
    # p_lb_ls = [279, 277, 273, 269, 257, 240]
    # p_ub_ls = [330, 330, 330, 330, 330, 330]

    m_lb_ls = [181.5, 181.5, 182, 183]
    m_ub_ls = [261.5, 260, 250, 230]
    p_lb_ls = [273, 269, 257, 240]
    p_ub_ls = [330, 330, 330, 330]

## LOAD DATA REQUIRED IN CALCULATIONS

# load ice latent heat of melting
# based on SeaFreeze derivation (>=240K) already in J/g
df_iceLh = pd.read_csv(ld+'SF_Lh.csv',skiprows=1, names=['T','Lh'])
iceLh_f = interp1d(df_iceLh['T'],df_iceLh['Lh'],fill_value='extrapolate')

# load ice specific heat
if ice_arg == 1: # if arg = 1, get ice cp from our own 0wt% experiments
    ice_data = base.DataPrep(2.5108,0)
    ice_data.import_data(dd)
    ice_data.calc_cp_pure()
    ice_data.save_cp_pure(od)
    iceCp_f = UnivariateSpline(ice_data.df_cp_p['T(K)'], ice_data.df_cp_p['cp(J/gK)'], k=3,ext=0)
    iceCp_err_interp = interp1d(ice_data.df_cp_p['T(K)'],(ice_data.df_cp_p['cp_err(%)']/100),fill_value='extrapolate')
else: # ice specific heat based on Feistel and Wagner 2006
    df_iceCp = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
    iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)
    #TODO interp error for FW2006

## DATA PROCESSING

for idx, wt in enumerate(wt_ls):
    # idx = 0 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging

    m = mass_ls[idx]
    m_lb = m_lb_ls[idx]
    m_ub = m_ub_ls[idx]
    p_lb = p_lb_ls[idx]
    p_ub = p_ub_ls[idx]

    data = base.DataPrep(m,wt)
    data.import_data(dd)

    # calculate and save Cp
    data.calc_cp_melt(iceLh_f,iceCp_f,iceCp_err_interp)
    data.calc_cp_pure()
    data.save_cp_melt(od)
    data.save_cp_pure(od)

    # plot Cp plots
    data.plot_cp('cp',data.df_cp_m,data.df_cp_p,fd,save=1,ylim=[0,12])
    data.plot_cp('cp',data.df_cp_m,data.df_cp_p,fd,save=1)
    # data.plot_cp('cpm',data.df_cp_m,data.df_cp_p,fd,save=1,ylim=[0,120])
    # data.plot_cp('cpm',data.df_cp_m,data.df_cp_p,fd,save=1)
    data.plot_hb(fd,save=1)

    # cut Cp
    data.cut_savgol(data.df_cp_m,m_lb,od,fd,plot=1)
    # data.cut_cp(data.df_cp_m, m_lb, m_ub, od, fd, test=1)
    data.cut_cp(data.df_cp_p, p_lb, p_ub, od, fd, test=0)


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