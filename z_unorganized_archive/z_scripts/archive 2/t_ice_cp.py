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
from scipy.interpolate import UnivariateSpline

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
ice_data = base.DataPrep(2.5108,0)
ice_data.import_data(dd)
ice_data.calc_cp_pure()
ice_exp = UnivariateSpline(ice_data.df_cp_p['T(K)'], ice_data.df_cp_p['cp(J/gK)'], k=3,ext=0)

# ice specific heat based on Feistel and Wagner 2006
df_iceCp = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
ice_FW = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)

# ice specific heat based on SeaFreeze
df_iceCp_SF = pd.read_csv(ld+'SF_cp_ice.csv', names=['T','Cp'],header=0)
ice_SF = interp1d(df_iceCp_SF['T'],df_iceCp_SF['Cp'])
##
plt.figure()
plt.scatter(ice_data.df_cp_p['T(K)'],ice_data.df_cp_p['cp(J/gK)'],0.1,'gray')
plt.plot(ice_data.df_cp_p['T(K)'],ice_exp(ice_data.df_cp_p['T(K)']))
plt.plot(df_iceCp['T'],ice_FW(df_iceCp['T']))
plt.plot(df_iceCp_SF['T'],ice_SF(df_iceCp_SF['T']))
base.show_plot_max()
## debug end line
print('Finished running script')