"""
created by Bing Hong CHUA 25Jul22

script objective:
plot custom plots of processed calorimetry data
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
plt.rcParams.update({'figure.dpi': '600'})

# set working, figure and output directories if required
wd = r".."
dd = r"i_data_processed/"
fd = r"o_specificHeat/"
od = r"i_data_processed/"
os.chdir(wd)

## IMPORT DATA FILES

# wtpct_ls = [2,4,6,8,10]
# file_ls = ["1.76wt%_climb_proc.txt", "3.52wt%_climb_proc.txt", "6wt%_climb_proc.txt", "8wt%_climb_proc.txt", "10wt%_climb_proc.txt"]
# mass_ls = [4.8389, 4.7459, 4.5983, 4.5500, 4.3987]

wtpct_ls = [1.42,2.84,5.11,6.88,8.73,26.96,2.842]
mass_ls = [4887.3,4840.2,4629.3,4585.8,4520.2,3777.8,0.9642]
# wtpct_ls = [1.76]
# mass_ls = [4887.3]

file_ls = [dd+f'{i}wt%_climb_proc.txt' for i in wtpct_ls]
mass_ls = [x*0.001 for x in mass_ls]

# --- import datafiles of samples --- #
col_names = ['index','time(s)','furnanceT(K)','sampleT(K)','Q(mW)','Q_corrected(mW)','Fc']
files = enumerate(file_ls, 1)
dfs= {idx: pd.read_csv(f,skiprows=1,names=col_names,delim_whitespace=True) for idx, f in files}


# col_names = ['index','time(s)','furnanceT(K)','sampleT(K)','Q(mW)','Q_corrected(mW)','Fc']
# df0 = pd.read_csv(dd+file_ls[0], skiprows=1, names=col_names, delim_whitespace=True)
# df2 = pd.read_csv(file_ls[1], skiprows=1, names=col_names, delim_whitespace=True)
# df4 = pd.read_csv(file_ls[2], skiprows=1, names=col_names, delim_whitespace=True)
# df6 = pd.read_csv(file_ls[3], skiprows=1, names=col_names, delim_whitespace=True)
# df8 = pd.read_csv(file_ls[4], skiprows=1, names=col_names, delim_whitespace=True)
# df10 = pd.read_csv(file_ls[5], skiprows=1, names=col_names, delim_whitespace=True)
# df0old = pd.read_csv(dd+r"0.25Kmin/0wt%_climb_proc.txt", skiprows=1, names=col_names, delim_whitespace=True)

# import debug line
print('Finished imports')

## ADH DOUBLE PEAK

plt.figure()
plt.plot(dfs[1]['sampleT(K)'].iloc[2788:3064],dfs[1]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[0]))
plt.plot(dfs[2]['sampleT(K)'].iloc[2788:3064],dfs[2]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[1]))
plt.plot(dfs[3]['sampleT(K)'].iloc[2788:3064],dfs[3]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[2]))
plt.plot(dfs[4]['sampleT(K)'].iloc[2788:3064],dfs[4]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[3]))
plt.plot(dfs[5]['sampleT(K)'].iloc[2788:3064],dfs[5]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[4]))
plt.plot(dfs[6]['sampleT(K)'].iloc[2788:3064],dfs[6]['Q_corrected(mW)'].iloc[2788:3064],label='{}wt%'.format(wtpct_ls[5]))
plt.plot(dfs[7]['sampleT(K)'].iloc[2400:2650],dfs[7]['Q_corrected(mW)'].iloc[2400:2650],label='2.84wt% 1mL')
plt.xlabel('Temperature (K)', fontweight='bold')
plt.ylabel('Heat Flow (mW)', fontweight='bold')
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'ADH_melt_double.png'
plt.savefig(saveFig)
plt.close()

## LIQUIDUS PEAKS

plt.close()
plt.figure()
plt.plot(dfs[1]['sampleT(K)'].iloc[5500:-400],dfs[1]['Q_corrected(mW)'].iloc[5500:-400],label='{}wt%'.format(wtpct_ls[0]))
plt.plot(dfs[2]['sampleT(K)'].iloc[5500:-400],dfs[2]['Q_corrected(mW)'].iloc[5500:-400],label='{}wt%'.format(wtpct_ls[1]))
plt.plot(dfs[3]['sampleT(K)'].iloc[5500:-400],dfs[3]['Q_corrected(mW)'].iloc[5500:-400],label='{}wt%'.format(wtpct_ls[2]))
plt.plot(dfs[4]['sampleT(K)'].iloc[5500:-400],dfs[4]['Q_corrected(mW)'].iloc[5500:-400],label='{}wt%'.format(wtpct_ls[3]))
plt.plot(dfs[5]['sampleT(K)'].iloc[5500:-400],dfs[5]['Q_corrected(mW)'].iloc[5500:-400],label='{}wt%'.format(wtpct_ls[4]))
plt.plot(dfs[7]['sampleT(K)'].iloc[5000:-990],dfs[7]['Q_corrected(mW)'].iloc[5000:-990],label='2.84wt% 1mL',color='tab:grey')
plt.xlabel('Temperature (K)', fontweight='bold')
plt.ylabel('Heat Flow (mW)', fontweight='bold')
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'liq_melt.png'
plt.savefig(saveFig)
plt.close()

## 3.52wt% SUBPLOTS


fig, ax = plt.subplots(2,1, sharex='col',sharey=True)
l1 = ax[0].plot(df4['sampleT(K)'], df4['Q_corrected(mW)'], 'b', label='3.52wt% 1')
l2 = ax[1].plot(df03['sampleT(K)'],df03['Q_corrected(mW)'],'b',label='3.52wt% 2')
ax[0].legend(fontsize=3,loc='lower left')
ax[1].legend(fontsize=3,loc='lower left')

fig.add_subplot(1, 1, 1, frame_on=False)
plt.tick_params(labelcolor="none", which='both',bottom=False, left=False, top=False, right=False)
plt.xlabel('Temperature (K)', fontweight='bold')
plt.ylabel('Heat Flow (mW)', fontweight='bold')
plt.show()
saveFig = 'all_3.52wt%_test_separatePlots.png'
plt.savefig(saveFig)
plt.close()

## 0wt% OVERLAP

fig = plt.figure()
plt.plot(df0['sampleT(K)'], df0['Q_corrected(mW)'], 'k', label='0.1K/min')
plt.plot(df0old['sampleT(K)'],df0old['Q_corrected(mW)']/2,'r',label='0.25K/min (cor)')
plt.legend()
plt.xlim([270,290])
plt.xlabel('T (K)', fontweight='bold')
plt.ylabel('Q (mW)', fontweight='bold')
# plt.show()
saveFig = 'x_0wt%compare3.png'
plt.savefig(saveFig)
plt.close()

## ALL HEALTFLOW SUBPLOT

fig, ax = plt.subplots(3,2, sharex='col',sharey=True)
ax[0,0].plot(df0['sampleT(K)'], df0['Q_corrected(mW)'], label='0wt%')
l1 = ax[1,0].plot(df2['sampleT(K)'], df2['Q_corrected(mW)'], 'y', label='1.76wt%')
l2 = ax[2,0].plot(df4['sampleT(K)'], df4['Q_corrected(mW)'], 'b', label='3.52wt%')
l3 = ax[0,1].plot(df6['sampleT(K)'], df6['Q_corrected(mW)'], 'r', label='6wt%')
l4 = ax[1,1].plot(df8['sampleT(K)'], df8['Q_corrected(mW)'], 'g', label='8wt%')
l5 = ax[2,1].plot(df10['sampleT(K)'],df10['Q_corrected(mW)'],'m',label='10wt%')
# ax[0,0].set_ylim([-700,100])
ax[0,0].legend(fontsize=3,loc='lower left')
ax[0,1].legend(fontsize=3,loc='lower left')
ax[1,0].legend(fontsize=3,loc='lower left')
ax[1,1].legend(fontsize=3,loc='lower left')
ax[2,0].legend(fontsize=3,loc='lower left')
ax[2,1].legend(fontsize=3,loc='lower left')
fig.legend(fontsize=6)
#
# Adding a plot in the figure which will encapsulate all the subplots with axis showing only
fig.add_subplot(1, 1, 1, frame_on=False)

# Hiding the axis ticks and tick labels of the bigger plot
plt.tick_params(labelcolor="none", which='both',bottom=False, left=False, top=False, right=False)

# Adding the x-axis and y-axis labels for the bigger plot
plt.xlabel('T (K)', fontweight='bold')
plt.ylabel('Q (mW)', fontweight='bold')

# plt.show()
saveFig = sd + 'all_3.52wt%_test.png'
plt.savefig(saveFig)
plt.close()

##
fig, ax = plt.subplots(3,2, sharex='col',sharey=True)
# ax[0,0].plot(df0['sampleT(K)'],df0['Q_corrected(mW)'],label='0wt%')
l1 = ax[1,0].plot(df2['sampleT(K)'], df2['Q_corrected(mW)'], 'y', label='1.76wt%')
l2 = ax[2,0].plot(df4['sampleT(K)'], df4['Q_corrected(mW)'], 'b', label='3.52wt%')
l3 = ax[0,1].plot(df6['sampleT(K)'], df6['Q_corrected(mW)'], 'r', label='6wt%')
l4 = ax[1,1].plot(df8['sampleT(K)'], df8['Q_corrected(mW)'], 'g', label='8wt%')
l5 = ax[2,1].plot(df10['sampleT(K)'],df10['Q_corrected(mW)'],'m',label='10wt%')
# ax[0,0].set_ylim([-700,100])
# ax[0,0].legend(fontsize=3,loc='lower left')
ax[0,1].legend(fontsize=3,loc='lower left')
ax[1,0].legend(fontsize=3,loc='lower left')
ax[1,1].legend(fontsize=3,loc='lower left')
ax[2,0].legend(fontsize=3,loc='lower left')
ax[2,1].legend(fontsize=3,loc='lower left')
fig.legend(fontsize=6)
#
# Adding a plot in the figure which will encapsulate all the subplots with axis showing only
fig.add_subplot(1, 1, 1, frame_on=False)

# Hiding the axis ticks and tick labels of the bigger plot
plt.tick_params(labelcolor="none", which='both',bottom=False, left=False, top=False, right=False)

# Adding the x-axis and y-axis labels for the bigger plot
plt.xlabel('T (K)', fontweight='bold')
plt.ylabel('Q (mW)', fontweight='bold')

plt.show()
# saveFig = sd + 'all_3.52wt%_test.png'
# plt.savefig(saveFig)
# plt.close()

## FC OVERLAP PLOT

fig = plt.figure(1)
plt.plot(df2['Fc'], df2['sampleT(K)'], 'b-', label='1.6wt%')
plt.plot(df4['Fc'], df4['sampleT(K)'], 'r-', label='3.2wt%')
plt.plot(df6['Fc'], df6['sampleT(K)'], 'y-', label='6.0wt%')
plt.plot(df8['Fc'], df8['sampleT(K)'], 'g-', label='8.0wt%')
plt.plot(df10['Fc'],df10['sampleT(K)'],'m-',label='10.0wt%')
plt.xlabel('Fraction of Crystal')
plt.ylabel('Temperature (K)')
plt.xlim([0,1])
plt.legend()
# plt.show()
saveFig = sd + 'FcPlot_T-Fc.png'
plt.savefig(saveFig)
plt.close()

fig = plt.figure(2)
plt.plot(df2['sampleT(K)'], df2['Fc'], 'b-', label='1.6wt%')
plt.plot(df4['sampleT(K)'], df4['Fc'], 'r-', label='3.2wt%')
plt.plot(df6['sampleT(K)'], df6['Fc'], 'y-', label='6.0wt%')
plt.plot(df8['sampleT(K)'], df8['Fc'], 'g-', label='8.0wt%')
plt.plot(df10['sampleT(K)'],df10['Fc'],'m-',label='10.0wt%')
plt.ylabel('Fraction of Crystal')
plt.xlabel('Temperature (K)')
plt.ylim([0,1])
plt.legend()
# plt.show()
saveFig = sd + 'FcPlot_Fc-T.png'
plt.savefig(saveFig)
plt.close()

# debug end line
print('Finished running script')