"""
created by Bing Hong CHUA 12Oct22

script objective:
plot different ice Cp values
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import base
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

dd = r"i_data_processed/"
ld = r"i_data_literature/"
sd = r"t_ice_cp/"
fd = r"t_ice_cp/"

## SCRIPT ARGUMENTS


## DATA PROCESSING

# load ice specific heat
# if arg = 1, get ice cp from our own 0wt% experiments
# ice_data = pd.read_csv(dd + '0wt%_cp_pure_2.5108g.csv', header=0)
# ice_data = pd.read_csv(dd + 'ice_cp_calib.csv', header=0)
ice_data = pd.read_csv(dd + 'ice_cp_nocalib.csv', header=0)

ice_data = ice_data.copy().iloc[64:]


# ice_exp = UnivariateSpline(ice_data['T(K)'], ice_data['cp(J/gK)'], k=3, ext=0)
ice_exp_c, pcov = np.polyfit(ice_data['T(K)'], ice_data['cp(J/gK)'],1,cov=1)
iceCp_err_interp = interp1d(ice_data['T(K)'], (ice_data['cp_err(%)'] / 100), fill_value='extrapolate')

# ice specific heat based on Feistel and Wagner 2006
df_iceCp = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
ice_FW = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)

# ice specific heat based on SeaFreeze
df_iceCp_SF = pd.read_csv(ld+'SF_cp_ice.csv', names=['T','Cp'],header=0)
ice_SF = interp1d(df_iceCp_SF['T'],df_iceCp_SF['Cp'])

water_data = pd.read_csv(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\t_ice_cp\0wt%_cp_pure_2.5108g_all.csv', header=0)


##
dev = 100 * (ice_data['cp(J/gK)'] - ice_FW(ice_data['T(K)']))/ice_FW(ice_data['T(K)'])

def func(x, a, b):
    return a*x + b
# dev_coeff, _ = curve_fit(func,ice_data['T(K)'],dev)
dev_coeff, _ = curve_fit(func,ice_data['T(K)'],dev,bounds=([0,-np.inf, ], [0.000000001,np.inf]))

# dev_coeff = np.polyfit(ice_data['T(K)'],dev,1) #POLYNOMIAL
Trange = np.arange(184.5,312.5,0.01)
dev_val = np.polyval(dev_coeff,Trange)
dev_f = interp1d(Trange,dev_val)

##

plt.figure()
plt.scatter(ice_data['T(K)'],dev,1,color='k')
plt.plot(Trange,dev_val)

plt.axhline(y=0,color='k')
plt.ylabel('Deviation (%)')
plt.xlabel('T (K)')
# base.show_plot_max()
plt.savefig(fd + 'iceCp_dev.png')

plt.figure()
plt.scatter(ice_data['T(K)'],ice_data['cp(J/gK)'],1,'gray',label='3 K averaged data')
plt.plot(ice_data['T(K)'],ice_data['cp(J/gK)']*(1+ice_data['cp_err(%)']/100),linewidth=.5,color='tab:gray')
plt.plot(ice_data['T(K)'],ice_data['cp(J/gK)']*(1-ice_data['cp_err(%)']/100),linewidth=.5,color='tab:gray',label='1-sigma')
plt.plot(ice_data['T(K)'],ice_exp_c[0]*ice_data['T(K)']+ice_exp_c[1],label='linear fit')
plt.plot(ice_data['T(K)'],ice_FW(ice_data['T(K)']),label='Feistel & Wagner 2006')
plt.legend()
# plt.plot(ice_data['T(K)'],ice_SF(ice_data['T(K)']))
# base.show_plot_max()
plt.savefig(fd + 'iceCp_new.png')

fig = plt.figure()
plt.plot(water_data['T(K)'].iloc[174:],water_data['cp(J/gK)'].iloc[174:])
plt.savefig(fd + 'waterCp_new.png')

## debug end line
print('Finished running script')