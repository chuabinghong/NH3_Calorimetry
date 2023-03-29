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

#
# figure_size = plt.gcf().get_size_inches()
# factor = 0.8
# plt.gcf().set_size_inches(factor * figure_size)
## FILE DIRECTORY SETTINGS

wd = r"../../"
os.chdir(wd)

dd = r"i_data/water_tests/"
sd = r"i_data/water_tests/"
fd = r"i_data/water_tests/"
ld = r"i_data_literature/"

## SCRIPT ARGUMENTS

yr_ls = [2014,2015,2022]
m_ls = [6,6,2.5108]

a_range = np.linspace(-100,100,200)
b_range = np.linspace(-100,100,200)

## DATA PROCESSING

# Save FW06 cp for conversion later on
df_FW_raw = pd.read_csv(ld+'ice_Cp.csv', names=['T','Cp'])
f_FW_Cp = interp1d(df_FW_raw['T'],df_FW_raw['Cp']/1000)

df = { '2014': pd.read_csv(dd + f'2014_hf_prepped.csv', header=0,usecols=['sampleT(K)','Q_corrected(mW)']),
       '2015': pd.read_csv(dd + f'2015_hf_prepped.csv', header=0,usecols=['sampleT(K)','Q_corrected(mW)']),
       '2022': pd.read_csv(dd + f'2022_hf_prepped.csv', header=0,usecols=['sampleT(K)','Q_corrected(mW)'])}


df['2014']['m'] = 6
df['2015']['m'] = 6
df['2022']['m'] = 2.5108

df['2014']['yr'] = 2014
df['2015']['yr'] = 2015
df['2022']['yr'] = 2022


# get FW Q values
def f_get_FW_Q(df):
    cp = f_FW_Cp(df['sampleT(K)'])
    df['FW_Q(mW)'] = -(df['m'] * 1 * cp / (10 * 60) * 1000)

for yr in df:
    f_get_FW_Q(df[yr])
    df[yr]['dev'] = df[yr]['Q_corrected(mW)'] - df[yr]['FW_Q(mW)']
    df[yr]['rel_dev'] = df[yr]['dev']/df[yr]['FW_Q(mW)']*100


df_all = pd.concat([df['2014'],df['2015'],df['2022']])
dev_coeff = np.polyfit(df_all['FW_Q(mW)'],df_all['rel_dev'],1) #POLYNOMIAL
dev_val = np.polyval(dev_coeff,df_all['FW_Q(mW)'])


w = 3.3
fig,ax = plt.subplots(2,1,sharex=True,figsize=(w, w * 3 / 3))


# ax[0].set_aspect('equal', adjustable='box')
ax[0].scatter(df_all['FW_Q(mW)'],df_all['Q_corrected(mW)'], 1,color='#4575b4', label='Calibration Data')
ax[0].plot([-30,0],[-30,0],'k--',linewidth=.5)
ax[0].set_xlim(-22,-3)
ax[0].set_ylim(-22,-3)
ax[0].set_ylabel('Heat Flow Data (mW)')

# ax[1].set_aspect('equal', adjustable='box')
ax[1].scatter(df_all['FW_Q(mW)'], df_all['rel_dev'],1,color='#4575b4',label='Calibration Data')
ax[1].plot(df_all['FW_Q(mW)'], dev_val,color='k',label='Correction Function')
ax[1].axhline(0,linewidth=.5, color='k')
ax[1].set_ylabel(r'100($\frac{\mathrm{d}q_{data}}{\mathrm{d}t}$ - $\frac{\mathrm{d}q_{FW06}}{\mathrm{d}t}$)/$\frac{\mathrm{d}q_{FW06}}{\mathrm{d}t}$ (%)')

ax[0].text(0.025, 0.98, 'a)', horizontalalignment='left', verticalalignment='top', transform=ax[0].transAxes)
ax[1].text(0.025, 0.96, 'b)', horizontalalignment='left', verticalalignment='top', transform=ax[1].transAxes)

ax[0].legend(prop={'size':4},loc='lower right')
ax[1].legend(prop={'size':4},bbox_to_anchor=(0.995,.05),loc='lower right')

base.splt_axis_label(fig,'Heat Flow FW06 (mW)','')
plt.tight_layout()
# base.show_plot_max()
plt.savefig(fd + 'Correction_equal.png')

# np.savetxt(sd + "calib_coeff.csv", dev_coeff, delimiter=",")



##

# get FW Q values
def f_get_FW_Q(df):
    cp = f_FW_Cp(df['sampleT(K)'])
    df['FW_Q(mW)'] = -(df['m'] * 1 * cp / (10 * 60) * 1000)

fig,ax = plt.subplots(1,2)

for yr in df:
    f_get_FW_Q(df[yr])
    df[yr]['dev'] = df[yr]['Q_corrected(mW)'] - df[yr]['FW_Q(mW)']

    ax[0].scatter(df[yr]['sampleT(K)'], df[yr]['Q_corrected(mW)'],1,label=f'{yr}')
    ax[0].plot(df[yr]['sampleT(K)'], df[yr]['FW_Q(mW)'])

    ax[1].scatter(df[yr]['FW_Q(mW)'],df[yr]['Q_corrected(mW)'], 1, label=f'{yr}')

df_all = pd.concat([df['2014'],df['2015'],df['2022']])
dev_coeff = np.polyfit(df_all['FW_Q(mW)'],df_all['dev'],1) #POLYNOMIAL
dev_val = np.polyval(dev_coeff,df_all['FW_Q(mW)'])

ax[0].set_aspect('equal', adjustable='box')
ax[0].scatter(df_all['FW_Q(mW)'],df_all['Q_corrected(mW)'], 1, label=f'{yr}')



ax[0].plot([-30,0],[-30,0],'k')

ax[1].set_xlabel('Heat Flow FW06 (mW)')
ax[1].set_ylabel('Heat Flow Data (mW)')
corr_ori = df_all['Q_corrected(mW)'].corr(df_all['FW_Q(mW)'])
# plt.text(-20,-5,'correlation before correction = {:.5f}'.format(corr_ori),fontsize=3)

ax[1].set_xlim(-22,-3)
ax[1].set_ylim(-22,-3)

ax[0].plot(df[yr]['sampleT(K)'], df[yr]['FW_Q(mW)'],'k',alpha=0,label='corresponding FS06')
leg = ax[0].legend(prop={'size':3.5})
for lh in leg.legendHandles:
    lh.set_alpha(1)
ax[1].legend(prop={'size':4})

ax[0].set_xlabel('Temperature (K)')
ax[0].set_ylabel('Heat Flow (mW)')


plt.tight_layout()
base.show_plot_max()
plt.savefig(fd + 'HF_plots_2.png')

np.savetxt(sd + "calib_coeff.csv", dev_coeff, delimiter=",")