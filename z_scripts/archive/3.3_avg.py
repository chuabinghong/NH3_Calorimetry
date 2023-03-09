"""
created by Bing Hong CHUA 21Sepl22

script objective:
calculate plot rolling averages
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns
import os
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

## SCRIPT SETTINGS

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE SETTINGS

# set working and save directories
dd = r"../i_data_processed/"
sd = r"../i_data_processed/"
fd = r"../t_SMA/"

# experimental values
wtpct_ls = [1.42,2.84,5.11,6.88,8.73,26.96]
mass_ls_in = [4887.3,4840.2,4629.3,4585.8,4520.2,3777.8]

# script input values
SMA_window = 1 #36 for 1K, 178 for 5K, 354 for 10K

# prepare file lists to load
file_ls = [dd+f'{i}wt%_cp_cut.csv' for i in wtpct_ls]
fileL_ls =[dd+f'{i}wt%_cp_cut_liq.csv' for i in wtpct_ls]

##

df1 = pd.read_csv(dd+'5.11wt%_cp_cut.csv')
df2 = pd.read_csv(dd+r'0.25Kmin/6wt%_cp.csv')
df1['cp_SMA_5K'] = df1['cp(J/gK)'].rolling(178,center=False).mean()
df2['cp_SMA_5K'] = df2['cp(J/gK)'].rolling(104,center=False).mean()
plt.plot(df1['T(K)'], df1['cp_SMA_5K'],label='0.1K/min')
plt.plot(df2['T(K)'], df2['cp_SMA_5K']+1,label='0.25K/min')
plt.ylim(3.9,4.4)
plt.show()

##
# iterate through all samples
for samp_id, samp in enumerate(wtpct_ls):
    # samp = 5.11
    # samp_id = 2
    df = pd.read_csv(fileL_ls[samp_id])

    ## CALC SIMPLE MOVING AVERAGE & SPLINES

    df['cp_SMA_1K'] = df['cp(J/gK)'].rolling(36,center=False).mean()
    df['cp_SMA_5K'] = df['cp(J/gK)'].rolling(178,center=False).mean()
    df['cp_SMA_10K'] = df['cp(J/gK)'].rolling(354,center=False).mean()
    spl_uni = UnivariateSpline(df['T(K)'], df['cp(J/gK)'], k=3)

    # # plot SMA with Data
    # plt.figure('indiv')
    # plt.plot(df['T(K)'], df['cp(J/gK)'], 'o', markersize=0.2, alpha=0.3, label='data')
    # plt.plot(df['T(K)'], df['cp_SMA_1K'],label='moving average (1K)')
    # plt.plot(df['T(K)'], df['cp_SMA_5K'],label='moving average (5K)')
    # plt.plot(df['T(K)'], df['cp_SMA_10K'],label='moving average (10K)')
    # plt.plot(df['T(K)'], spl_uni(df['T(K)']), label='cubic spline')
    # plt.xlabel('T (K)')
    # plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    # plt.title('{}wt% curve fit check'.format(samp))
    # plt.legend(markerscale=2,prop={'size': 5})
    # # plt.show()
    # saveFig = fd + '{}wt%_cp_fit.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close('indiv')

    ## CALCULATE DEVIATION FROM SPLINE

    # dev = df['cp(J/gK)']-spl_uni(df['T(K)'])
    #
    # sd = (np.sum((dev**2))/len(df))**.5
    #
    # plt.figure('dev')
    # plt.axhline(0, color='k')
    # plt.axhline(sd, color='k',linestyle='--')
    # plt.axhline(-sd, color='k',linestyle='--')
    # plt.plot(df['T(K)'],dev,'o',markersize=0.2)
    # plt.text(180,sd+0.05,'$\sigma$={:.3f}'.format(sd),fontsize=7)
    # plt.xlabel('T(K)')
    # plt.ylabel(r'cp deviation $(\frac{J}{g ⋅ K})$')
    # plt.title('{}wt% spline deviation'.format(samp))
    # # plt.show()
    # saveFig = fd + '{}wt%_curve_dev.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close('dev')
    #
    # plt.figure('dist')
    # sns.kdeplot(dev)
    # plt.axvline(sd,color='k',linestyle='--')
    # plt.axvline(-sd,color='k',linestyle='--')
    # plt.xlabel(r'cp deviation $(\frac{J}{g ⋅ K})$')
    # plt.title('{}wt% spline deviation distribution'.format(samp))
    # # plt.show()
    # saveFig = fd + '{}wt%_curve_dev_dist.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close('dist')

    ## PLOT SPLINES ON SAME GRAPH

    # plt.figure('spline')
    # plt.plot(df['T(K)'], spl_uni(df['T(K)']), label='cubic spline')
    if samp_id != 5:
        plt.figure('sma')
        plt.scatter(df['T(K)'], df['cp(J/gK)'],0.2,marker='x',alpha=0.2)
        plt.plot(df['T(K)'], df['cp_SMA_5K'],label='{}wt%'.format(samp))
plt.xlabel('T (K)')
plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
plt.title('Pure phase moving average (5k)')
plt.legend(markerscale=2,prop={'size': 5})
# plt.show()
saveFig = fd + 'purePhase_SMA.png'.format(samp)
plt.savefig(saveFig)
plt.close('sma')

print('end')