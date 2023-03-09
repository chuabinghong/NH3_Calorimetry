"""
created by Bing Hong CHUA 21Sepl22

script objective:
template
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns
import os
from scipy.interpolate import interp1d
from numpy.fft import fft, ifft
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

## SCRIPT SETTINGS

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE SETTINGS

# set working and save directories
wd = r"../"
dd = r"i_data_processed/"
sd = r"i_data_processed/"
fd = r""
os.chdir(wd)
# experimental values
wtpct_ls = [0]
mass_ls_in = [4887.3,4840.2,4629.3,4585.8,4520.2,3777.8]

# script input values
SMA_window = 1 #36 for 1K, 178 for 5K, 354 for 10K

# prepare file lists to load
file_ls = [dd+f'{i}wt%_cp_cut.csv' for i in wtpct_ls]
fileL_ls =[dd+f'{i}wt%_cp_cut_liq.csv' for i in wtpct_ls]


##
start_offset = 6000  # FOR 0.1KMIN
end_offset = 12900
start_offset = 8500  # FOR 0.25KMIN
end_offset = 12100

df = pd.read_csv(r'i_data/0.25Kmin/blank.csv', skiprows=8)

df2 = pd.read_csv(r'i_data/0.10Kmin/Blank.csv', skiprows=8)
df = df.iloc[8500:12100].copy()
df2 = df2.iloc[6000:12900].copy()
df.reset_index(inplace=True,drop=True)
df2.reset_index(inplace=True,drop=True)

a = fft(df['HeatFlow(mW)'])
tp = len(df['HeatFlow(mW)'])
values = np.arange(int(tp))
f = 1/(df['Time(s)'].iloc[1] - df['Time(s)'].iloc[0])
timePeriod  = tp/f

frequencies = values/timePeriod
plt.figure()
plt.plot(frequencies,abs(a))


plt.figure()
plt.plot(df.index,df['HeatFlow(mW)'])
plt.plot(df2.index,df2['HeatFlow(mW)'])



## ICE CP

col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
df = pd.read_csv(dd+r'0.25Kmin/0wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)
m = 4.9859
m = 4.985

Cp = np.zeros([len(df), 1]);
Cp[:] = np.nan
Cp_T = np.zeros([len(df), 1])
for i in range(len(df) - 1):
    Q_integral = -(df['Q_corrected(mW)'].iloc[i + 1] / 1000 + df['Q_corrected(mW)'].iloc[i] / 1000) * (
            df['time(s)'].iloc[i + 1] - df['time(s)'].iloc[i]) / 2
    dT = df['sampleT(K)'].iloc[i + 1] - df['sampleT(K)'].iloc[i]

    Cp[i] = Q_integral / (m * dT)
    Cp_T[i] = (df['sampleT(K)'].iloc[i + 1] + df['sampleT(K)'].iloc[i]) / 2

df_cp = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(Cp)})

df_cp = df_cp[300:-700].copy()
df_cp2 = df_cp[100:-300].copy()

df_iceCp = pd.read_csv(dd+'ice_Cp.csv', names=['T','Cp'])
iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000,fill_value='extrapolate')

spl_uni = UnivariateSpline(df_cp['T(K)'],df_cp['cp(J/gK)'], k=3)
spl_uni2 = UnivariateSpline(df_cp['T(K)'],df_cp['cp(J/gK)'], k=3)

fig = plt.figure()
plt.plot(df_cp['T(K)'],df_cp['cp(J/gK)'], 'o', markersize=0.2, label='data')
# plt.plot(Cp_T, Cp, label='ice')
plt.plot(df_cp['T(K)'],spl_uni(df_cp['T(K)']),'-',label='spline')
plt.plot(df_cp['T(K)'],iceCp_f(df_cp['T(K)']),'-',label='Feistel & Wagner 2006')
plt.xlabel('T (K)')
plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
plt.title('specific heat of ice (0.25K/min) m=4.986')
# plt.ylim([0,3])
# plt.xlim([])
plt.legend(markerscale=5)
plt.show()
# saveFig = fd + 'ice_cp_0.25Kmin_m=4.986g.png'
# plt.savefig(saveFig)
# plt.close()

##

col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
df = pd.read_csv(dd+'0wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)
m = 2.511
# m = 1

Cp = np.zeros([len(df), 1]);
Cp[:] = np.nan

CpU = np.zeros([len(df), 1]);
CpU[:] = np.nan
CpD = np.zeros([len(df), 1]);
CpD[:] = np.nan
Cp_T = np.zeros([len(df), 1])
for i in range(len(df) - 1):
    Q_integral = -(df['Q_corrected(mW)'].iloc[i + 1] / 1000 + df['Q_corrected(mW)'].iloc[i] / 1000) * (
            df['time(s)'].iloc[i + 1] - df['time(s)'].iloc[i]) / 2
    dT = df['sampleT(K)'].iloc[i + 1] - df['sampleT(K)'].iloc[i]

    Cp[i] = Q_integral / (m * dT)
    Cp_T[i] = (df['sampleT(K)'].iloc[i + 1] + df['sampleT(K)'].iloc[i]) / 2

    Q_up = -((df['Q_corrected(mW)'].iloc[i + 1] +0.1) / 1000 + (df['Q_corrected(mW)'].iloc[i]+0.1) / 1000) * (
            df['time(s)'].iloc[i + 1] - df['time(s)'].iloc[i]) / 2
    CpU[i] = Q_up / ((m+0.001) * dT)

    Q_d = -((df['Q_corrected(mW)'].iloc[i + 1] -0.1) / 1000 + (df['Q_corrected(mW)'].iloc[i]-0.1) / 1000) * (
            df['time(s)'].iloc[i + 1] - df['time(s)'].iloc[i]) / 2
    CpD[i] = Q_d / ((m-0.001) * dT)

df_cp = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(Cp)})
df_cp = df_cp[700:-750].copy()

df_cpU = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(CpU)})
df_cpU = df_cpU[700:-750].copy()

df_cpD = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(CpD)})
df_cpD = df_cpD[700:-750].copy()

df_iceCp = pd.read_csv(dd+'ice_Cp.csv', names=['T','Cp'])
iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000,fill_value='extrapolate')

spl_uni = UnivariateSpline(df_cp['T(K)'],df_cp['cp(J/gK)'], k=3)
spl_uniU = UnivariateSpline(df_cpU['T(K)'],df_cpU['cp(J/gK)'], k=3)
spl_uniD = UnivariateSpline(df_cpD['T(K)'],df_cpD['cp(J/gK)'], k=3)

fig = plt.figure()
plt.plot(df_cp['T(K)'],df_cp['cp(J/gK)'], 'o', markersize=0.2, label='data')
# plt.plot(Cp_T, Cp, label='ice')
plt.plot(df_cp['T(K)'],spl_uni(df_cp['T(K)']),'-',label='spline')
plt.plot(df_cp['T(K)'],spl_uniU(df_cp['T(K)']),'--',label='min')
plt.plot(df_cp['T(K)'],spl_uniD(df_cp['T(K)']),'--',label='max')
plt.plot(df_cp['T(K)'],iceCp_f(df_cp['T(K)']),'-',label='Feistel & Wagner 2006')
plt.xlabel('T (K)')
plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
plt.title('specific heat of ice (0.1K/min), m=2.511g')
# plt.ylim([0,3])
# plt.xlim([])
plt.legend(markerscale=5)
# plt.show()
# saveFig = fd + 'ice_cp_m=2.511_unc.png'
# plt.savefig(saveFig)
# plt.close()
##

plt.figure()
plt.plot(df_cp['T(K)'],df_cp['cp(J/gK)'], 'o', markersize=0.2, label='0.1 K/min')
plt.plot(df_cp2['T(K)'],df_cp2['cp(J/gK)'], 'o', markersize=0.2, label='0.25 K/min')
plt.xlabel('T (K)')
plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
plt.legend(markerscale=3)
saveFig = fd + 'ice_cp_both.png'
plt.savefig(saveFig)
plt.close()
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