"""
created by Bing Hong CHUA 04Aug22

script objective:
plot representations of specific heat
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['scatter.edgecolors'] = 'k'
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode


# set working, figure and output directories if required
wd = r"../"
dd = r"i_data_processed/"
ddb = r"i_data_processed/0.25Kmin/"
fd = r"/"
os.chdir(wd)

## FUNCTIONS

## IMPORT DATA FILES

wtpct_ls = [1.6,3.52,6,8,10]
file_ls = ["1.6wt%_cp.csv", "3.2wt%_cp.csv", "6wt%_cp.csv", "8wt%_cp.csv", "10wt%_cp.csv"]

# 0.1Kmin data
df1 = pd.read_csv(dd+file_ls[0])
df1l = pd.read_csv(dd+"1.6wt%_cp_liquidPhase.csv")
df2 = pd.read_csv(dd+file_ls[1])
df2l = pd.read_csv(dd+"3.2wt%_cp_liquidPhase.csv")
# df3 = pd.read_csv(file_ls[2])
df4 = pd.read_csv(dd+file_ls[3])
df4l = pd.read_csv(dd+"8wt%_cp_liquidPhase.csv")
# df5 = pd.read_csv(file_ls[4])
df26 =  pd.read_csv(dd+'26.96wt%_cp.csv')
df26l =  pd.read_csv(dd+'26.96wt%_cp_liquidPhase.csv')

#0.25Kmin data
# df1b = pd.read_csv(ddb+file_ls[1])
# df1bl = pd.read_csv(ddb+"1.6wt%_cp_liquidPhase.csv")
# df4b = pd.read_csv(ddb+file_ls[3])
# df4bl = pd.read_csv(ddb+"8wt%_cp_liquidPhase.csv")

# literature data
df_33 = pd.read_csv(dd+'Chan_Giauque_0.33.csv')
df_49 = pd.read_csv(dd+'Hildenbrand_Giauque_0.49.csv')
df_59 = pd.read_csv(dd+'Hildenbrand_Giauque_0.59.csv')
df_61 = pd.read_csv(dd+'Hildenbrand_Giauque_0.61.csv')
df_64 = pd.read_csv(dd+'Hildenbrand_Giauque_0.64.csv')
df_65 = pd.read_csv(dd+'Hildenbrand_Giauque_0.65.csv')

# SF cp test data
df_SF8 = pd.read_csv(dd+'SF_cp_8wt%.csv')
df_SF27 = pd.read_csv(dd+'SF_cp_27wt%.csv')
## CUT TO REASONABLE REGION

dfCp = df26.copy()
# dfCp = dfCp.iloc[np.arange(190,2300)]

## PLOT CPM

fig = plt.figure(1)
plt.scatter(df4['T(K)'],df4['cpm(J/molK)'],0.5,linewidths=0.2,edgecolors='face',marker='x',label='8wt% (melt)')
plt.scatter(df4l['T(K)'],df4l['cpm(J/molK)'],0.5,linewidths=0.2,edgecolors='face',marker='x',label='8wt% (liq)')
plt.scatter(dfCp['T(K)'],dfCp['cpm(J/molK)'],1,linewidths=0.2,facecolor='grey',edgecolors='face',marker='x',label='newBottle (melt)')
plt.scatter(df26l['T(K)'],df26l['cpm(J/molK)'],1,linewidths=0.2,facecolor='coral',edgecolors='face',marker='x',label='newBottle (liq)')
plt.scatter(df_33['T(K)'],df_33['cpm(J/molK)'],10,linewidths=0.2,label='33wt%')
plt.scatter(df_49['T(K)'],df_49['cpm(J/molK)'],10,linewidths=0.2,marker='^',label='49wt%')
plt.scatter(df_59['T(K)'],df_59['cpm(J/molK)'],10,linewidths=0.2,marker='d',label='59wt%')
plt.scatter(df_61['T(K)'],df_61['cpm(J/molK)'],10,linewidths=0.2,marker='s',label='61wt%')
plt.scatter(df_64['T(K)'],df_64['cpm(J/molK)'],10,linewidths=0.2,marker='p',label='64wt%')
plt.scatter(df_65['T(K)'],df_65['cpm(J/molK)'],10,linewidths=0.2,marker='H',label='65wt%')
plt.xlabel('T (K)')
plt.ylabel('cpm (J $mol^{-1}$ $K^{-1}$)')
plt.ylim([20,90])
# plt.xlim([200,380])
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'CPM_test3.png'
plt.savefig(saveFig)
plt.close()

## PLOT CP

fig = plt.figure(1)
plt.scatter(df1['T(K)'],df1['cp(J/gK)'],1,linewidths=0.2,marker='x',label='1.6wt% (melt)')
plt.scatter(df2['T(K)'],df2['cp(J/gK)'],1,linewidths=0.2,marker='x',label='3.2wt% (melt)')
plt.scatter(df4['T(K)'],df4['cp(J/gK)'],1,linewidths=0.2,marker='x',label='8wt% (melt)')
plt.scatter(dfCp['T(K)'],dfCp['cp(J/gK)'],1,linewidths=0.2,marker='x',label='newBottle (melt)')
# plt.scatter(df1l['T(K)'],df1l['cp(J/gK)'],1,linewidths=0.2,marker='x',label='1.6wt% (liq)')
# plt.scatter(df2l['T(K)'],df2l['cp(J/gK)'],1,linewidths=0.2,marker='x',label='3.2wt% (liq)')
# plt.scatter(df4l['T(K)'],df4l['cp(J/gK)'],1,linewidths=0.2,marker='x',label='8wt% (liq)')
# plt.scatter(df26l['T(K)'],df26l['cp(J/gK)'],1,linewidths=0.2,marker='x',label='26.96wt% (liq)')
# plt.scatter(df_33['T(K)'],df_33['cp(J/gK)'],10,linewidths=0.2,marker='o',facecolor='dimgrey',label='33wt%')
# plt.scatter(df_49['T(K)'],df_49['cp(J/gK)'],10,linewidths=0.2,marker='^',facecolor='grey',label='49wt%')
# plt.scatter(df_59['T(K)'],df_59['cp(J/gK)'],10,linewidths=0.2,marker='d',facecolor='darkgrey',label='59wt%')
# plt.scatter(df_61['T(K)'],df_61['cp(J/gK)'],10,linewidths=0.2,marker='s',facecolor='lightgrey',label='61wt%')
# plt.scatter(df_64['T(K)'],df_64['cp(J/gK)'],10,linewidths=0.2,marker='p',facecolor='gainsboro',label='64wt%')
# plt.scatter(df_65['T(K)'],df_65['cp(J/gK)'],10,linewidths=0.2,marker='H',facecolor='white',label='65wt%')
# plt.plot(df_SF8['T(K)'],df_SF8['cp(J/gK)'],label='EOS 7.85wt% ',color='cadetblue')
# plt.plot(df_SF27['T(K)'],df_SF27['cp(J/gK)'],label='EOS 27.26wt% ',color='peru')

plt.xlabel('T (K)')
plt.ylabel('cp (J $g^{-1}$ $K^{-1}$)')
plt.ylim([1,5])
# plt.xlim([200,380])
plt.legend(prop={'size': 5})
plt.show()
# saveFig = 'CP_test_SF.png'
# plt.savefig(saveFig)
# plt.close()

## PLOT 0.25Kmin vs 0.1Kmin comparisons

plt.scatter(df1['T(K)'],df1['cp(J/gK)'],1.5,linewidths=0.1,marker='s',label='0.1K/min')
plt.scatter(df1l['T(K)'],df1l['cp(J/gK)'],1.5,linewidths=0.1,marker='s',label='0.1K/min')
plt.scatter(df1b['T(K)'],df1b['cp(J/gK)'],6,linewidths=0.1,marker='.',label='0.25K/min (pure phase)')
plt.scatter(df1bl['T(K)'],df1bl['cp(J/gK)'],6,linewidths=0.1,marker='.',label='0.25K/min (pure phase)')
plt.xlabel('T (K)')
plt.ylabel('cp (J $g^{-1}$ $K^{-1}$)')
plt.ylim([1,5])
# plt.xlim([200,380])
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'CPM_rateCompare_1.6wt%.png'
plt.savefig(saveFig)
plt.close()

## PLOT 3D SCATTER
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(df1['T(K)'],df1['X'],df1['cp(J/gK)'],marker='o',s=2,label='1.6wt%')
ax.scatter(df2['T(K)'],df2['X'],df2['cp(J/gK)'],marker='o',s=2,label='3.52wt%')
ax.scatter(df3['T(K)'],df3['X'],df3['cp(J/gK)'],marker='o',s=2,label='6wt%')
ax.scatter(df4['T(K)'],df4['X'],df4['cp(J/gK)'],marker='o',s=2,label='8wt%')
ax.scatter(df5['T(K)'],df5['X'],df5['cp(J/gK)'],marker='o',s=2,label='10wt%')
ax.set_xlabel('T (K)')
ax.set_ylabel('X')
ax.set_zlabel('cp (J/gK)')
ax.set_zlim([-5,20])
plt.legend()
# plt.show()
# saveFig = fd + 'cp_3dscatter.png'
# plt.savefig(saveFig)
# plt.close()

ax.view_init(elev=0, azim=270)
# plt.show()
saveFig = fd + 'cp_3dscatter_side.png'
plt.savefig(saveFig)
plt.close()
print('test')

## PLOT OVERLAP OF ALL CP

fig = plt.figure(1)
plt.plot(df1['T(K)'],df1['cp(J/gK)'],label='1.6wt%')
plt.plot(df2['T(K)'],df2['cp(J/gK)'],label='3.2wt%')
plt.plot(df3['T(K)'],df3['cp(J/gK)'],label='6wt%')
plt.plot(df4['T(K)'],df4['cp(J/gK)'],label='8wt%')
plt.plot(df5['T(K)'],df5['cp(J/gK)'],label='10wt%')
plt.xlabel('T (K)')
plt.ylabel('cp (J/gK)')
plt.ylim([-5,20])
plt.legend()
# plt.show()
saveFig = fd + 'cp_all.png'
plt.savefig(saveFig)
plt.close()

## PLOT OVERLAP OF ALL CP

fig = plt.figure(1)
plt.plot(df2['T(K)'],df2['cp(J/gK)'],label='3.2wt%')
plt.plot(df2['T(K)'].iloc[200],df2['cp(J/gK)'].iloc[200],'o',markersize=2)
plt.plot(df2['T(K)'].iloc[1400],df2['cp(J/gK)'].iloc[1400],'o',markersize=2)

plt.xlabel('T (K)')
plt.ylabel('cp (J/gK)')
plt.ylim([-5,20])
plt.legend()
# plt.show()
saveFig = fd + 'cp_test_10Aug.png'
plt.savefig(saveFig)
plt.close()



## debug end line
print('Finished running script')