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
import glob

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['scatter.edgecolors'] = 'k'
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

old_wt_arg = 1

# set working, figure and output directories if required
wd = r"../"
dd = r"i_data_processed/"
ddb = r"i_data_processed/0.25Kmin/"
fd = r"/"
os.chdir(wd)

## FUNCTIONS

## IMPORT DATA FILES

wtpct_ls = [1.42,2.84,5.11,6.88,26.96]
file_ls = [dd+f'{i}wt%_cp_cut.csv' for i in wtpct_ls]
fileL_ls = [dd+f'{i}wt%_cp_cut_liq.csv' for i in wtpct_ls]
files = enumerate(file_ls, 1)
filesL = enumerate(fileL_ls, 1)
dfs= {idx: pd.read_csv(f) for idx, f in files}
dfs_L = {idx: pd.read_csv(f) for idx, f in filesL}

if old_wt_arg ==1:
    dd = r"i_data_processed/old_wt%/"
    fd = r"/t_old_wt%/"
    wtpct_ls = [1.62, 3.24, 5.83, 7.85, 26.96]
    file_ls = [dd + f'{i}wt%_cp_cut.csv' for i in wtpct_ls]
    fileL_ls = [dd + f'{i}wt%_cp_cut_liq.csv' for i in wtpct_ls]
    files = enumerate(file_ls, 1)
    filesL = enumerate(fileL_ls, 1)
    dfs = {idx: pd.read_csv(f) for idx, f in files}
    dfs_L = {idx: pd.read_csv(f) for idx, f in filesL}

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

## PLOT CPM

## PLOT CP

fig = plt.figure()
plt.scatter(dfs[1]['T(K)'],dfs[1]['cp(J/gK)'],1,linewidths=0.15,label='{}wt% (melt)'.format(wtpct_ls[0]))
plt.scatter(dfs[2]['T(K)'],dfs[2]['cp(J/gK)'],1,linewidths=0.15,label='{}wt% (melt)'.format(wtpct_ls[1]))
plt.scatter(dfs[3]['T(K)'],dfs[3]['cp(J/gK)'],1,linewidths=0.15,label='{}wt% (melt)'.format(wtpct_ls[2]))
plt.scatter(dfs[4]['T(K)'],dfs[4]['cp(J/gK)'],1,linewidths=0.15,label='{}wt% (melt)'.format(wtpct_ls[3]))
plt.scatter(dfs[5]['T(K)'],dfs[5]['cp(J/gK)'],1,linewidths=0.15,label='{}wt% (melt)'.format(wtpct_ls[4]))
plt.scatter(dfs_L[1]['T(K)'],dfs_L[1]['cp(J/gK)'],1,marker='s',facecolor='tab:blue',linewidths=0.15)
plt.scatter(dfs_L[2]['T(K)'],dfs_L[2]['cp(J/gK)'],1,marker='s',facecolor='tab:green',linewidths=0.15)
plt.scatter(dfs_L[3]['T(K)'],dfs_L[3]['cp(J/gK)'],1,marker='s',facecolor='tab:orange',linewidths=0.15)
plt.scatter(dfs_L[4]['T(K)'],dfs_L[4]['cp(J/gK)'],1,marker='s',facecolor='tab:red',linewidths=0.15)
plt.scatter(dfs_L[5]['T(K)'],dfs_L[5]['cp(J/gK)'],1,marker='s',facecolor='tab:purple',linewidths=0.15)
plt.scatter(0,0,1,linewidths=0.15,marker='s',facecolor='white',label='liquid phase')
plt.plot(df_33['T(K)'],df_33['cp(J/gK)'],color='k',label='33wt% (literature)',zorder=0)
plt.plot(df_65['T(K)'],df_65['cp(J/gK)'],'--',color='k',label='65wt% (literature)',zorder=0)
# plt.scatter(df_49['T(K)'],df_49['cp(J/gK)'],10,linewidths=0.2,marker='^',facecolor='grey',label='49wt%')
# plt.scatter(df_59['T(K)'],df_59['cp(J/gK)'],10,linewidths=0.2,marker='d',facecolor='darkgrey',label='59wt%')
# plt.scatter(df_61['T(K)'],df_61['cp(J/gK)'],10,linewidths=0.2,marker='s',facecolor='lightgrey',label='61wt%')
# plt.scatter(df_64['T(K)'],df_64['cp(J/gK)'],10,linewidths=0.2,marker='p',facecolor='gainsboro',label='64wt%')
# plt.scatter(df_65['T(K)'],df_65['cp(J/gK)'],10,linewidths=0.2,marker='H',facecolor='white',label='65wt%')
plt.plot(df_SF8['T(K)'],df_SF8['cp(J/gK)'],label='EOS 7.85wt% ',color='cadetblue')
plt.plot(df_SF27['T(K)'],df_SF27['cp(J/gK)'],label='EOS 27.26wt% ',color='peru')
plt.title('cp all SF')
plt.xlabel('T (K)')
plt.ylabel('cp (J $g^{-1}$ $K^{-1}$)')
plt.ylim([2,7])
# plt.ylim([4,4.5])
# plt.xlim([200,380])
plt.legend(prop={'size': 3},markerscale=3)
# plt.show()
saveFig = 'CP_all_SF_old.png'
plt.savefig(saveFig)
plt.close()

## PLOT CP AND FIT

coef = np.polyfit(dfs[1]['T(K)'],dfs[1]['cp(J/gK)'],1); fit1 = np.poly1d(coef)
coef = np.polyfit(dfs[2]['T(K)'],dfs[2]['cp(J/gK)'],1); fit2 = np.poly1d(coef)
coef = np.polyfit(dfs[3]['T(K)'],dfs[3]['cp(J/gK)'],1); fit3 = np.poly1d(coef)
coef = np.polyfit(dfs[4]['T(K)'],dfs[4]['cp(J/gK)'],1); fit4 = np.poly1d(coef)
coef = np.polyfit(dfs[5]['T(K)'],dfs[5]['cp(J/gK)'],1); fit5 = np.poly1d(coef)

fig = plt.figure()
plt.scatter(dfs[1]['T(K)'],dfs[1]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs[2]['T(K)'],dfs[2]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs[3]['T(K)'],dfs[3]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs[4]['T(K)'],dfs[4]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs[5]['T(K)'],dfs[5]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.plot(dfs[1]['T(K)'], fit1(dfs[1]['T(K)']),color='tab:blue',label='{}wt% (melt)'.format(wtpct_ls[0]))
plt.plot(dfs[2]['T(K)'], fit2(dfs[2]['T(K)']),color='tab:green',label='{}wt% (melt)'.format(wtpct_ls[1]))
plt.plot(dfs[3]['T(K)'], fit3(dfs[3]['T(K)']),color='tab:orange',label='{}wt% (melt)'.format(wtpct_ls[2]))
plt.plot(dfs[4]['T(K)'], fit4(dfs[4]['T(K)']),color='tab:red',label='{}wt% (melt)'.format(wtpct_ls[3]))
plt.plot(dfs[5]['T(K)'], fit5(dfs[5]['T(K)']),color='tab:purple',label='{}wt% (melt)'.format(wtpct_ls[4]))
plt.title('cp melt phase linear fit')
plt.xlabel('T (K)')
plt.ylabel('cp (J $g^{-1}$ $K^{-1}$)')
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'CP_melt_fit_old.png'
plt.savefig(saveFig)
plt.close()

## PLOT CP AND FIT OF LIQUID PHASE

coef = np.polyfit(dfs_L[1]['T(K)'],dfs_L[1]['cp(J/gK)'],1); fit1 = np.poly1d(coef)
coef = np.polyfit(dfs_L[2]['T(K)'],dfs_L[2]['cp(J/gK)'],1); fit2 = np.poly1d(coef)
coef = np.polyfit(dfs_L[3]['T(K)'],dfs_L[3]['cp(J/gK)'],1); fit3 = np.poly1d(coef)
coef = np.polyfit(dfs_L[4]['T(K)'],dfs_L[4]['cp(J/gK)'],1); fit4 = np.poly1d(coef)
coef = np.polyfit(dfs_L[5]['T(K)'],dfs_L[5]['cp(J/gK)'],1); fit5 = np.poly1d(coef)

fig = plt.figure()
plt.scatter(dfs_L[1]['T(K)'],dfs_L[1]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs_L[2]['T(K)'],dfs_L[2]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs_L[3]['T(K)'],dfs_L[3]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.scatter(dfs_L[4]['T(K)'],dfs_L[4]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
# plt.scatter(dfs_L[5]['T(K)'],dfs_L[5]['cp(J/gK)'],1,linewidths=0.05,alpha=0.1)
plt.plot(dfs_L[1]['T(K)'], fit1(dfs_L[1]['T(K)']),color='tab:blue',label='{}wt% (melt)'.format(wtpct_ls[0]))
plt.plot(dfs_L[2]['T(K)'], fit2(dfs_L[2]['T(K)']),color='tab:green',label='{}wt% (melt)'.format(wtpct_ls[1]))
plt.plot(dfs_L[3]['T(K)'], fit3(dfs_L[3]['T(K)']),color='tab:orange',label='{}wt% (melt)'.format(wtpct_ls[2]))
plt.plot(dfs_L[4]['T(K)'], fit4(dfs_L[4]['T(K)']),color='tab:red',label='{}wt% (melt)'.format(wtpct_ls[3]))
# plt.plot(dfs_L[5]['T(K)'], fit5(dfs_L[5]['T(K)']),color='tab:purple',label='{}wt% (melt)'.format(wtpct_ls[4]))
plt.title('cp pure phase linear fit')
plt.xlabel('T (K)')
plt.ylabel('cp (J $g^{-1}$ $K^{-1}$)')
plt.legend(prop={'size': 5})
# plt.show()
saveFig = 'CP_liquid_fit.png'
plt.savefig(saveFig)
plt.close()


## PLOT 0.25Kmin vs 0.1Kmin comparisons

plt.scatter(dfs[1]['T(K)'],dfs[1]['cp(J/gK)'],1.5,linewidths=0.1,marker='s',label='0.1K/min')
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
ax.scatter(dfs[1]['T(K)'],dfs[1]['X'],dfs[1]['cp(J/gK)'],marker='o',s=2,label='1.6wt%')
ax.scatter(dfs[2]['T(K)'],dfs[2]['X'],dfs[2]['cp(J/gK)'],marker='o',s=2,label='3.52wt%')
ax.scatter(dfs[3]['T(K)'],dfs[3]['X'],dfs[3]['cp(J/gK)'],marker='o',s=2,label='6wt%')
ax.scatter(dfs[4]['T(K)'],dfs[4]['X'],dfs[4]['cp(J/gK)'],marker='o',s=2,label='8wt%')
ax.scatter(dfs[5]['T(K)'],dfs[5]['X'],dfs[5]['cp(J/gK)'],marker='o',s=2,label='10wt%')
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
plt.plot(dfs[1]['T(K)'],dfs[1]['cp(J/gK)'],label='1.6wt%')
plt.plot(dfs[2]['T(K)'],dfs[2]['cp(J/gK)'],label='3.2wt%')
plt.plot(dfs[3]['T(K)'],dfs[3]['cp(J/gK)'],label='6wt%')
plt.plot(dfs[4]['T(K)'],dfs[4]['cp(J/gK)'],label='8wt%')
plt.plot(dfs[5]['T(K)'],dfs[5]['cp(J/gK)'],label='10wt%')
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
plt.plot(dfs[2]['T(K)'],dfs[2]['cp(J/gK)'],label='3.2wt%')
plt.plot(dfs[2]['T(K)'].iloc[200],dfs[2]['cp(J/gK)'].iloc[200],'o',markersize=2)
plt.plot(dfs[2]['T(K)'].iloc[1400],dfs[2]['cp(J/gK)'].iloc[1400],'o',markersize=2)

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