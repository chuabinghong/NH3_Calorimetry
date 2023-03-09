"""
created by Bing Hong CHUA 27Sep22

script objective:
from raw data obtained from the calorimeter
output corrected heatflow data to be used in further analysis
plot thermograms, crystal fraction curves and identify peaks
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

dd = r"i_data/0.10Kmin/"
sd = r"t_vary_X/data/"
fd = r"t_vary_X/"

## SCRIPT ARGUMENTS
peak_arg = 0 #1 if we are getting peaks

ramp_rate = 0.1
m = 4.8402
wt_ls = np.arange(2,3,0.1) # based off liquidus alignment
m_lb = 181.5
m_ub = 262
p_lb = 273.5
p_ub = 310

# ## DATA PROCESSING
#
# # import blank
# blank = base.DataRaw(0,0,ramp_rate)
# blank.df = pd.read_csv(dd + 'Blank.csv', skiprows=8) #TODO for 0.1 K/min
#
# for idx, wt in enumerate(wt_ls):
#     # idx = 3 #TODO for single data debugging
#     # wt = wt_ls[idx] #TODO for single data debugging
#
#     # data processing
#     data = base.DataRaw(m, wt,ramp_rate)
#     data.import_raw(dd)
#     data.correct_HF(blank.df) # calibration with blank
#
#     data.plot_Q_T(fd,1)
#     data.calc_Fc() # calculating crystal fraction
#     data.plot_Fc_T(fd,1)
#
#     data.save_data(sd)
#
#     print('Finished loop {}/{}'.format(idx,len(wt_ls)-1))
#

##

ld = r"i_data_literature/"
dd = sd
od = sd

## SCRIPT ARGUMENTS
ice_arg = 1 # 1 if using ice cp from our own 0wt% experiments

## LOAD DATA REQUIRED IN CALCULATIONS

# load ice latent heat of melting
# based on SeaFreeze derivation (>=240K) already in J/g
df_iceLh = pd.read_csv(ld+'SF_Lh.csv',skiprows=1, names=['T','Lh'])
iceLh_f = interp1d(df_iceLh['T'],df_iceLh['Lh'],fill_value='extrapolate')

# load ice specific heat
if ice_arg == 1: # if arg = 1, get ice cp from our own 0wt% experiments
    ice_data = base.DataPrep(2.5108,0)
    ice_data.import_data(r'i_data_processed/')
    ice_data.calc_cp_pure()
    iceCp_f = UnivariateSpline(ice_data.df_cp_p['T(K)'], ice_data.df_cp_p['cp(J/gK)'], k=3,ext=0)
else: # ice specific heat based on Feistel and Wagner 2006
    df_iceCp = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
    iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)

##
for idx, wt in enumerate(wt_ls):

    data = base.DataPrep(m, wt)
    data.import_data(dd)
    data.calc_cp_melt(iceLh_f, iceCp_f)
    data.calc_cp_pure()

    data.plot_cp('cp', data.df_cp_m, data.df_cp_p, fd, save=1, ylim=[0, 12])
    data.plot_cp('cp', data.df_cp_m, data.df_cp_p, fd, save=1)
    data.plot_cp('cpm', data.df_cp_m, data.df_cp_p, fd, save=1, ylim=[0, 120])
    data.plot_cp('cpm', data.df_cp_m, data.df_cp_p, fd, save=1)
    data.plot_hb(fd, save=1)

    data.cut_cp(data.df_cp_m, m_lb, m_ub, od, fd, test=0)
    data.cut_cp(data.df_cp_p, p_lb, p_ub, od, fd, test=0)

    print('Finished for loop {} / {}'.format(idx, len(wt_ls) - 1))

    ## DATA PROCESSING
    plt.close()

for idx, wt in enumerate(wt_ls):

    data = base.DataCP(m, wt)
    data.import_data(dd)
    data.calc_spline()

    plt.figure('combined')
    data.plot_spl2(pure=0)

plt.xlabel('Temperature (K)');
plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
plt.title('Mass Fraction from Liquidus Alignment')
plt.ylim([2, 7])
plt.legend(prop={'size': 4})
plt.savefig('t_cp_vary_X.png')
# base.show_plot_max()

print('Finished for loop {} / {}'.format(idx, len(wt_ls) - 1))

## debug end line
print('Finished running script')

## debug end line

print('Finished running script')
