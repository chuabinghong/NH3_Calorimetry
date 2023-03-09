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

## IMPORT ICE LH

# ice latent heat of melting based on SeaFreeze derivation (>=240K) already in J/g
df_iceLh = pd.read_csv(dd+'SF_Lh.csv',skiprows=1, names=['T','Lh'])
iceLh_f = interp1d(df_iceLh['T'],df_iceLh['Lh'],fill_value='extrapolate')

## IMPORT ICE CP

col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
df_0 = pd.read_csv(dd + '0wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)
m_0 = 2.511

Cp = np.zeros([len(df_0), 1]);
Cp[:] = np.nan
Cp_T = np.zeros([len(df_0), 1])

for i in range(len(df_0) - 1):
    Q_integral = -(df_0['Q_corrected(mW)'].iloc[i + 1] / 1000 + df_0['Q_corrected(mW)'].iloc[i] / 1000) * (
            df_0['time(s)'].iloc[i + 1] - df_0['time(s)'].iloc[i]) / 2
    dT = df_0['sampleT(K)'].iloc[i + 1] - df_0['sampleT(K)'].iloc[i]

    Cp[i] = Q_integral / (m_0 * dT)
    Cp_T[i] = (df_0['sampleT(K)'].iloc[i + 1] + df_0['sampleT(K)'].iloc[i]) / 2

df_cp = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(Cp)})
df_cp = df_cp[700:-750].copy()

df_iceCp = pd.read_csv(dd + 'ice_Cp.csv', names=['T', 'Cp'])
own_ice = UnivariateSpline(df_cp['T(K)'], df_cp['cp(J/gK)'], k=3, ext=0)

df_iceCp = pd.read_csv(dd + 'ice_Cp.csv', names=['T', 'Cp'])
lit_ice = interp1d(df_iceCp['T'], df_iceCp['Cp'] / 1000)

## LOAD TEST DATA

m = 4.5202
col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
df = pd.read_csv(dd+'8.73wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)

## DO CP CALCULATIONS W OWN ICE
dfB = df[df["Fc"].first_valid_index():df["Fc"].last_valid_index()].copy()

Cp = np.zeros([len(dfB), 1]);
Cp[:] = np.nan
Cp_T = np.zeros([len(dfB), 1])
for i in range(len(dfB) - 1):
    Q_integral = -(dfB['Q_corrected(mW)'].iloc[i + 1] / 1000 + dfB['Q_corrected(mW)'].iloc[i] / 1000) * (
            dfB['time(s)'].iloc[i + 1] - dfB['time(s)'].iloc[i]) / 2
    dfFc = dfB['Fc'].iloc[i] - dfB['Fc'].iloc[i + 1]
    dT = dfB['sampleT(K)'].iloc[i + 1] - dfB['sampleT(K)'].iloc[i]

    # Cp[i] = (Q_integral - dfFc * m * latent_heat - dfB['Fc'].iloc[i] * m * iceCp_f(dfB['sampleT(K)'].iloc[i+1]) * dT) / (m * (1 - dfB['Fc'].iloc[i+1]) * dT)
    Cp[i] = (Q_integral - dfFc * m * iceLh_f(dfB['sampleT(K)'].iloc[i + 1]) - dfB['Fc'].iloc[i] * m * own_ice(
        dfB['sampleT(K)'].iloc[i + 1]) * dT) / (m * (1 - dfB['Fc'].iloc[i + 1]) * dT)
    Cp_T[i] = (dfB['sampleT(K)'].iloc[i + 1] + dfB['sampleT(K)'].iloc[i]) / 2


## DO CP CALCULATIONS W LIT ICE
dfB = df[df["Fc"].first_valid_index():df["Fc"].last_valid_index()].copy()

Cp2 = np.zeros([len(dfB), 1]);
Cp2[:] = np.nan
Cp_T2 = np.zeros([len(dfB), 1])
for i in range(len(dfB) - 1):
    Q_integral = -(dfB['Q_corrected(mW)'].iloc[i + 1] / 1000 + dfB['Q_corrected(mW)'].iloc[i] / 1000) * (
            dfB['time(s)'].iloc[i + 1] - dfB['time(s)'].iloc[i]) / 2
    dfFc = dfB['Fc'].iloc[i] - dfB['Fc'].iloc[i + 1]
    dT = dfB['sampleT(K)'].iloc[i + 1] - dfB['sampleT(K)'].iloc[i]

    # Cp[i] = (Q_integral - dfFc * m * latent_heat - dfB['Fc'].iloc[i] * m * iceCp_f(dfB['sampleT(K)'].iloc[i+1]) * dT) / (m * (1 - dfB['Fc'].iloc[i+1]) * dT)
    Cp2[i] = (Q_integral - dfFc * m * iceLh_f(dfB['sampleT(K)'].iloc[i + 1]) - dfB['Fc'].iloc[i] * m * lit_ice(
        dfB['sampleT(K)'].iloc[i + 1]) * dT) / (m * (1 - dfB['Fc'].iloc[i + 1]) * dT)
    Cp_T2[i] = (dfB['sampleT(K)'].iloc[i + 1] + dfB['sampleT(K)'].iloc[i]) / 2

## PLOT BOTH

plt.figure()
plt.plot(Cp_T2, Cp2, 'o', markersize=0.2, label='Feistel & Wagner 2006')
plt.plot(Cp_T, Cp, 'o', markersize=0.2, label='Experimental')
plt.xlabel('T (K)')
plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
plt.title('Test of specific heat used for ice')
plt.ylim([2,7])
plt.legend(markerscale=5)
# plt.show()
saveFig = 't_ice_cp_test'
plt.savefig(saveFig)
plt.close()

print('end')
##

plt.plot(dfB['sampleT(K)'],lit_ice(dfB['sampleT(K)']),'-',label='spline')
##