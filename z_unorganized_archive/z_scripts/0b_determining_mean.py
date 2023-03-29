"""
created by Bing Hong CHUA 27Sep22
edited 03Jan23

script objective:
plot the difference between the two blanks.
identify a reasonable T_averaged bin size based off of sigma of difference
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import base
from scipy.signal import argrelextrema
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

dd = r"i_data/blanks/"
sd = r"i_data/blanks/"
fd = r"i_data/blanks/"

## SCRIPT ARGUMENTS
ramp_rate = 0.1
T_bin = 1
T_bin3 = 5
## DATA PROCESSING

# import blanks
b26 = base.DataRaw(0,'blank',ramp_rate,12981)
b26.df = pd.read_csv(dd + 'Blank_26C.csv', skiprows=8)
b26.correct_HF(b26.df)
b26.df2 = b26.df2.reset_index(drop=True)


b46 = base.DataRaw(0,'blank',ramp_rate,13066)
b46.df = pd.read_csv(dd + 'Blank_46C.csv', skiprows=8)
b46.correct_HF(b46.df)
b46.df2 = b46.df2.reset_index(drop=True)

# get deviation

def get_dev(df26,df46):
    # take difference and get sigma from raw blanks
    int26 = interp1d(df26['Sample Temperature(K)'], df26['Q_Corrected'], bounds_error=False, fill_value=0)
    int46 = interp1d(df46['Sample Temperature(K)'], df46['Q_Corrected'], bounds_error=False, fill_value=0)

    df_sort = df26.iloc[(df26['Sample Temperature(K)'] - 176).abs().argsort()[:1]]
    start_idx = df_sort.index.values.astype(int)[0]

    dev = int46(df26['Sample Temperature(K)'].iloc[start_idx:-1]) - int26(
        df26['Sample Temperature(K)'].iloc[start_idx:-1])

    out = {"dev": dev,
           "start_idx":start_idx}

    return(out)

out = get_dev(b26.df2,b46.df2)
sigma = np.std(out['dev'])

plt.plot(b26.df2['Sample Temperature(K)'].iloc[out['start_idx']:-1],out['dev'],label='original')

# get mean and dev in mean
b26.mean_smoothing(b26.df2,T_bin,sd)
b46.mean_smoothing(b46.df2,T_bin,sd)
out2 = get_dev(b26.df_mean,b46.df_mean)
plt.plot(b26.df_mean['Sample Temperature(K)'].iloc[out2['start_idx']:-1],out2['dev'],label=f'{T_bin} K Bin')

# get mean and dev in mean
T_bin2 = 3
b26.mean_smoothing(b26.df2,T_bin2,sd)
b46.mean_smoothing(b46.df2,T_bin2,sd)
out2 = get_dev(b26.df_mean,b46.df_mean)
plt.plot(b26.df_mean['Sample Temperature(K)'].iloc[out2['start_idx']:-1],out2['dev'],label=f'{T_bin2} K Bin')

# get mean and dev in mean
b26.mean_smoothing(b26.df2,T_bin3,sd)
b46.mean_smoothing(b46.df2,T_bin3,sd)
out2 = get_dev(b26.df_mean,b46.df_mean)
plt.plot(b26.df_mean['Sample Temperature(K)'].iloc[out2['start_idx']:-1],out2['dev'],label=f'{T_bin3} K Bin')



plt.title("46C blank - 26C blank")
plt.xlabel('Heat Flow (mW)')
plt.ylabel('T (K)')
plt.axhline(y=sigma,color='k',linestyle='--')
plt.axhline(y=-sigma,color='k',linestyle='--')
plt.annotate(f'sigma = {sigma:.3f}', xy=(0.05, 0.90), xycoords='axes fraction')
plt.legend(prop={'size': 5})
# base.show_plot_max()
plt.savefig(fd + 'blank_dev.png')



## debug end line

print('Finished running script')
