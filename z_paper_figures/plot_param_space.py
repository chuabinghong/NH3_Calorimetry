"""
created by Bing Hong CHUA 04Aug22

script objective:
plot representations of specific heat
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Rectangle
import pandas as pd
import base
from scipy.interpolate import interp1d

plt.style.use(['science', 'nature', 'no-latex'])
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['scatter.edgecolors'] = 'k'
plt.rcParams.update({'figure.dpi': '600'})
mpl.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

ld = r"i_data_literature/"
fd = r"z_paper_figures/"

## LOAD LITERATURE DATA

# df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
df_TRF = pd.read_csv(ld + 'Tillner-Roth_Friend_liq.csv', header=0)
CG = {'33': pd.read_csv(ld + 'Chan_Giauque_0.33.csv', header=0)}
HG = { '49': pd.read_csv(ld + 'Hildenbrand_Giauque_0.49.csv', header=0),
           '59': pd.read_csv(ld + 'Hildenbrand_Giauque_0.59.csv', header=0),
           '61': pd.read_csv(ld + 'Hildenbrand_Giauque_0.61.csv', header=0),
           '64': pd.read_csv(ld + 'Hildenbrand_Giauque_0.64.csv', header=0),
           '65': pd.read_csv(ld + 'Hildenbrand_Giauque_0.65.csv', header=0)}

df_WK = pd.read_csv(ld + 'Wrewsky_Kaigorodoff.csv', header=0)
##
liq_T = np.linspace(176,273,100)
Xl_calc_v = np.vectorize(base.Xl_calc)
liq_wt = Xl_calc_v(liq_T) *100
liq_T =np.append(liq_T,273.15)
liq_wt =np.append(liq_wt,0)
# liq_T_int = interp1d(liq_wt,liq_T)
##
fig,ax = plt.subplots()
a=ax.add_patch(Rectangle((176,0), 300-176, 15-0,facecolor='#7fcdbb',alpha=0.5))
# b=ax.add_patch(Rectangle((180,2.9), 315-180, 26.9-2.9,facecolor='tab:blue',alpha=0.5,label='This Work'))


wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K
mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    m = mass_ls[idx]
    data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data_processed\{wt}wt%_cp_cut_melt_{m}g.csv',
        header=0)
    plt.scatter(data['T(K)'],np.repeat(wt,len(data)),2,color='#225ea8',edgecolors='none',marker='o')
    data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data_processed\{wt}wt%_cp_cut_pure_{m}g.csv',
        header=0)
    plt.scatter(data['T(K)'],np.repeat(wt,len(data)),2,color='#225ea8',edgecolors='none',marker='o')

plt.scatter(0,0,2,color='#225ea8',edgecolors='none',marker='o',label='This work')
plt.scatter(CG['33']['T(K)'],CG['33']['X']*100,marker='o', facecolors='k', edgecolors='k', linewidths=.1,label='Chan & Giauque')
plt.scatter(HG['49']['T(K)'],HG['49']['X']*100,marker='o', facecolors='w', edgecolors='k',linewidths=.1,label='Hildenbrand & Giauque')
plt.scatter(HG['59']['T(K)'],HG['59']['X']*100,marker='o', facecolors='w', edgecolors='k',linewidths=.1)
plt.scatter(HG['61']['T(K)'],HG['61']['X']*100,marker='o', facecolors='w', edgecolors='k',linewidths=.1)
plt.scatter(HG['64']['T(K)'],HG['64']['X']*100,marker='o', facecolors='w', edgecolors='k',linewidths=.1)
plt.scatter(HG['65']['T(K)'],HG['65']['X']*100,marker='o', facecolors='w', edgecolors='k',linewidths=.1)
plt.scatter(df_WK['T(K)'],df_WK['wt'],marker='x', facecolors='k', linewidths=.3,label='Wrewsky & Kaigorodoff')
plt.plot(liq_T,liq_wt,'k--',label='liquidus')
# ax.fill_between(liq_T,liq_wt,15)

plt.xlim(176,319.55)
plt.ylim(0,100)
plt.xlabel('T (K)')
plt.ylabel('$w_{NH_3}$ (wt%)')
# plt.legend(bbox_to_anchor=(0.5,-0.4),loc='lower center')
plt.legend(fontsize=5)
# plt.show()
saveFig = fd + 'param_space_2.png'
plt.savefig(saveFig)
plt.close()

## debug end line
print('Finished running script')



## UNC STUFF

x = [100,500,1000]
y = [3,1,0.6]

# plt.plot(x,y)
# plt.show()
## PARAM SPACE

T_l = np.arange(100,300,0.01)
Xl_f = np.vectorize(base.Xl_calc)
Xl = Xl_f(T_l)
