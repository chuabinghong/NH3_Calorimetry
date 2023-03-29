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
import numpy as np


plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FILE DIRECTORY SETTINGS

wd = r"../../"
os.chdir(wd)

dd = r"i_data/water_tests/"
sd = r"i_data/water_tests/"
fd = r"i_data/water_tests/"
ld = r"i_data_literature/"

## SCRIPT ARGUMENTS
ice_arg = 0 # 1 if using ice cp from our own 0wt% experiments
range_arg = 0 #1 if using -196C to 46C
calib_arg = 0 #1 if applying linear calibration offset, 2 if constant offset

ramp_rate = 0.1
mass_ls = [0.190,0.340,6.000]
mass_ls = [6.001]
wt_ls = [0, 0, 0]  # based off liquidus alignment offset
wt_ls = [0]  # based off liquidus alignment offset

## LOAD DATA REQUIRED IN CALCULATIONS

# load ice latent heat of melting
# based on SeaFreeze derivation (>=240K) already in J/g
df_iceLh = pd.read_csv(ld+'SF_Lh.csv',skiprows=1, names=['T','Lh'])

iceLh_coeff = np.polyfit(df_iceLh['T'],df_iceLh['Lh'],1) #POLYNOMIAL
Trange = np.arange(150,300,0.01)
iceLh_val = np.polyval(iceLh_coeff,Trange)
iceLh_f = interp1d(Trange,iceLh_val)

## PLOT ICE
# iceLh_2 = interp1d(df_iceLh['T'],df_iceLh['Lh'],fill_value='extrapolate') #OLD
#
# iceLh_coeff = np.polyfit(df_iceLh['T'],df_iceLh['Lh'],1) #POLYNOMIAL
# Trange = np.arange(150,300,0.01)
# iceLh_val = np.polyval(iceLh_coeff,Trange)
# iceLh_3 = interp1d(Trange,iceLh_val)

# plt.plot(df_iceLh['T'],df_iceLh['Lh'],linewidth=3, label='SeaFreeze output')
# plt.plot(Trange,iceLh_f(Trange),label='Quadratic Fit')
# plt.plot(Trange,iceLh_3(Trange),label='Linear Fit')
# plt.plot(Trange,iceLh_2(Trange),label='interp1d (old)')
# plt.xlim([176,310])
# plt.ylim([0,400])
# plt.xlabel('T (K)')
# plt.ylabel('Latent Heat (J/g)')
# plt.legend()

##
# load ice specific heat
if ice_arg == 1: # if arg = 1, get ice cp from our own 0wt% experiments
    ice_data = pd.read_csv(dd + '0wt%_cp_pure_2.5108g.csv', header=0)
    iceCp_f = UnivariateSpline(ice_data['T(K)'], ice_data['cp(J/gK)'], k=3,ext=0)
    iceCp_err_interp = interp1d(ice_data['T(K)'],(ice_data['cp_err(%)']/100),fill_value='extrapolate')
else: # ice specific heat based on Feistel and Wagner 2006
    df_iceCp = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
    iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)
    Trange = np.arange(176,273,1)
    iceCp_err_interp =interp1d(Trange,np.repeat(2/100,len(Trange)),fill_value='extrapolate')
    #TODO interp error for FW2006

# load mixing excess enthalpy
df_mix = pd.read_csv(ld+'Mixing.csv', header=0)
mix_coeff = np.polyfit(df_mix['wt'],df_mix['Hmix(J/g)'],2)
## DATA PROCESSING


for idx, wt in enumerate(wt_ls):
    # idx = 3 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging

    m = mass_ls[idx]

    data = base.DataPrep(m,wt)
    data.import_data(dd)

    # calculate and save Cp
    if calib_arg ==1:
        df_calib = pd.read_csv(dd + 'calibration.csv')
    elif calib_arg ==2:
        df_calib = pd.DataFrame()
        df_calib['T(K)'] = np.arange(100, 340, 0.5)
        df_calib['dev%'] = 5.33479154
    else:
        df_calib = pd.DataFrame()
        df_calib['T(K)'] = np.arange(100,340,0.5)
        df_calib['dev%'] = 0
    data.calc_cp_pure(range_arg,df_calib)
    data.save_cp_pure(sd)



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