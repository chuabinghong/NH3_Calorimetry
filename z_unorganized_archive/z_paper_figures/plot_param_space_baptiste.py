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

ld = r"i_data_literature/Baptiste/"
fd = r"z_paper_figures/"

## LOAD LITERATURE DATA

# df_SF = pd.read_csv(ld + 'SF_cp_liq.csv', header=0)
# df_TRF = pd.read_csv(ld + 'Tillner-Roth_Friend_liq.csv', header=0)
# CG = {'33': pd.read_csv(ld + 'Chan_Giauque_0.33.csv', header=0)}
# HG = { '49': pd.read_csv(ld + 'Hildenbrand_Giauque_0.49.csv', header=0),
#            '59': pd.read_csv(ld + 'Hildenbrand_Giauque_0.59.csv', header=0),
#            '61': pd.read_csv(ld + 'Hildenbrand_Giauque_0.61.csv', header=0),
#            '64': pd.read_csv(ld + 'Hildenbrand_Giauque_0.64.csv', header=0),
#            '65': pd.read_csv(ld + 'Hildenbrand_Giauque_0.65.csv', header=0)}
#
# df_WK = pd.read_csv(ld + 'Wrewsky_Kaigorodoff.csv', header=0)

df_chan = pd.read_csv(ld + 'chan_1964.csv', header=0)
df_allred = pd.read_csv(ld + 'allred_1981.csv', header=0)
df_chernenkaya = pd.read_csv(ld + 'chernenkaya_1971.csv', header=0)
df_fujita = pd.read_csv(ld + 'fujita_2008.csv', header=0)
df_hildenbrand = pd.read_csv(ld + 'hildenbrand_1953.csv', header=0)
df_wrewsky = pd.read_csv(ld + 'wrewsky_1924.csv', header=0)
df_prieto = pd.read_csv(ld + 'prieto_2022_noLiq.csv', header=0)

marker_ls = ['s','d','P','X','^','v','p']

## ADD MELTING CURVES

# define curves
def curve2(T):

    lqd_coeff = [-29.17, 184.7 - T]  # CROFT 1987
    roots = np.roots(lqd_coeff) # get roots of the polynomial
    roots_real = np.real(roots[np.iscomplex(roots) == False]) #remove imaginary roots
    root_actual = roots_real[(0.329 < roots_real) & (roots_real <= 0.353)] #keep roots within peritectic
    if len(root_actual) == 0:
        root_actual = np.nan
    else:
        root_actual = root_actual[0]
    return root_actual

def curve3(T_range):

    out = np.zeros([0,2]); out[:] = np.nan
    for a,T in enumerate(T_range):
        lqd_coeff = [8064, -14479, 8585.6, -1789.77, 248.07 - T]  # CROFT 1987
        roots = np.roots(lqd_coeff) # get roots of the polynomial
        roots_real = np.real(roots[np.iscomplex(roots) == False]) #remove imaginary roots
        root_actual = roots_real[(0.353 < roots_real) & (roots_real <= 0.572)]*100 #keep roots within peritectic
        if len(root_actual) == 0:
            root_actual = np.nan
        else:
            T2 = np.repeat(T,len(root_actual))
            arr = np.array([T2,root_actual])
            out = np.append(out,np.transpose(arr),axis=0)
    return out

def curve4(T_range):
    out = np.zeros([0, 2]);
    out[:] = np.nan
    for a,T in enumerate(T_range):
        lqd_coeff = [-3524, 11357, -14140, 7867.6, -1435.42 - T]  # CROFT 1987
        roots = np.roots(lqd_coeff) # get roots of the polynomial
        roots_real = np.real(roots[np.iscomplex(roots) == False]) #remove imaginary roots
        root_actual = roots_real[(0.572 < roots_real) & (roots_real <= 0.803)]*100 #keep roots within peritectic
        if len(root_actual) == 0:
            root_actual = np.nan
        else:
            T2 = np.repeat(T,len(root_actual))
            arr = np.array([T2,root_actual])
            out = np.append(out,np.transpose(arr),axis=0)
    return out

def curve5(T):

    lqd_coeff = [-241.7, 504.3, -239.48, 172.44 - T]  # CROFT 1987
    roots = np.roots(lqd_coeff) # get roots of the polynomial
    roots_real = np.real(roots[np.iscomplex(roots) == False]) #remove imaginary roots
    root_actual = roots_real[(0.803 < roots_real) & (roots_real <= 1)] #keep roots within peritectic
    if len(root_actual) == 0:
        root_actual = np.nan
    else:
        root_actual = root_actual[0]
    return root_actual

def lineplot(fx,liq_T):
    Xl_calc_v = np.vectorize(fx)
    liq_wt = Xl_calc_v(liq_T) *100
    return(liq_T,liq_wt)


##
fig,ax = plt.subplots()

wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K
mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

# PLOT DATA ABOVE LIQUIDUS
for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    m = mass_ls[idx]
    data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data_processed\{wt}wt%_cp_cut_melt_{m}g.csv',
        header=0)
    # plt.scatter(data['T(K)'],np.repeat(wt,len(data)),2,color='#225ea8',edgecolors='none',marker='o')
    data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data_processed\{wt}wt%_cp_cut_pure_{m}g.csv',
        header=0)
    plt.scatter(data['T(K)'],np.repeat(wt,len(data)),3,color='#225ea8',edgecolors='None', linewidths=.5,marker='o',zorder=5)


plt.scatter(0,0,3,edgecolors='None', linewidths=.5,marker='o',label='This work')
#
# # PLOT DATA ALONG LIQUIDUS
liq_T2 = np.arange(184.5,233.5+.5,.5)
xx,yy = lineplot(base.Xl_calc,liq_T2)
plt.plot(xx,yy,color='#225ea8',linewidth=1.2)
plt.scatter(xx,yy,3,color='#225ea8',edgecolors='None', linewidths=.5,marker='o',zorder=5)
# plt.scatter(xx,yy,3,edgecolors='None', linewidths=.5,marker='o',zorder=5,alpha=0)

# PLOT LITERATURE
plt.scatter(df_chan['T(K)'],df_chan['massfrac']*100,marker=marker_ls[0],label= 'CG64',edgecolors='k', linewidths=.5)
plt.scatter(df_hildenbrand['T(K)'],df_hildenbrand['massfrac']*100,marker=marker_ls[1],label= 'HG53',edgecolors='k', linewidths=.5)
plt.scatter(df_wrewsky['T(K)'],df_wrewsky['massfrac']*100,marker=marker_ls[2],label= 'WK24',edgecolors='k', linewidths=.5)
plt.scatter(df_chernenkaya['T(K)'],df_chernenkaya['massfrac']*100,marker=marker_ls[3],label= 'C71',edgecolors='k', linewidths=.5)
plt.scatter(df_allred['T(K)'],df_allred['massfrac']*100,marker=marker_ls[4],label= 'AW81',edgecolors='k', linewidths=.5)
plt.scatter(df_fujita['T(K)'],df_fujita['massfrac']*100,marker=marker_ls[5],label= 'F08',edgecolors='k', linewidths=.5)
plt.scatter(df_prieto['T(K)'],df_prieto['massfrac']*100,marker=marker_ls[6],label= 'P22',facecolor='#e377c2',edgecolors='k', linewidths=.5)

# PLOT ALL MELTING CURVES
liq_T = np.linspace(160,273.14,1000)
x1,y1 = lineplot(base.Xl_calc,liq_T)
x1 = np.append(x1, 273.15)
y1 = np.append(y1, 0)
x2,y2 = lineplot(curve2,liq_T)
x5,y5 = lineplot(curve5,liq_T)

T_range = np.linspace(172,210,500)
o3 = curve3(T_range)
o4 = curve4(T_range)

o1 = np.column_stack([x1,y1])
o2 = np.column_stack([x2,y2])
o5 = np.column_stack([x5,y5])

liq = np.concatenate([o1,o2,o3,o4,o5])
liq = liq[~np.isnan(liq).any(axis=1)]
liq = liq[liq[:, 1].argsort()]

plt.plot(liq[:,0],liq[:,1],'k',linewidth=1.5,label='Melting Curve')

# PLOT DESIRED PARAMETER RANGE
T_start = base.T_calc(.15)
T_range = np.linspace(T_start,273.14,1000)
xs,ys = lineplot(base.Xl_calc,T_range)
xs = np.append(xs, 273.15)
ys = np.append(ys, 0)

a=ax.fill_between(xs, ys, 15,facecolor='#94e5d2',alpha=1)
a.zorder=-2
a=ax.add_patch(Rectangle((273.15,0), 300-273.15, 15-0,facecolor='#94e5d2',alpha=1))
a.zorder=-2

plt.plot(x1,y1,'#94e5d2',linewidth=.5)
plt.plot(1,1,'#94e5d2',linewidth=1.5,label='Relevant Range')



plt.xlim(172,319.55)
plt.ylim(0,100)
plt.xlabel('Temperature (K)')
plt.ylabel('$NH_3$ Mass Fraction (wt%)')
# plt.legend(bbox_to_anchor=(0.5,-0.4),loc='lower center')
plt.legend(fontsize=4)
# base.show_plot_max()
saveFig = fd + 'param_space_28mar.png'
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
