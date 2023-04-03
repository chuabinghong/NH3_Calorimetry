import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import base
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

ld = r"i_data_literature/"
dd = r"i_data_processed/"
od = r"i_data_processed/"
fd = r"o_specificHeat/"

wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K
mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

colour_ls = ['#d73027','#f46d43','#fdae61','#e8e884','#abd9e9','#74add1','#4575b4']
marker_ls = ['o','s','d','P','X','^','v']
## test deviation with linear ice Lh

dd2 = dd + "ice_lh_test/"

for idx, wt in enumerate(wt_ls):
    # idx = 1 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]
    marker = marker_ls[idx]
    m = mass_ls[idx]

    data = base.DataCP(m, wt)
    data.import_data(dd, mean=0)
    data2 =  base.DataCP(m, wt)
    data2.import_data(dd2, mean=0)

    dev = 100 * (1-data.df_m['cp(J/gK)']/data2.df_m['cp(J/gK)'])

    plt.scatter(data.df_m['T(K)'],dev,1,c= colour,label=f'{wt} wt%')
    plt.axhline(0, linewidth=.3, color='k')
    plt.xlabel("T(K)")
    plt.ylabel("deviation (%)")
    plt.legend(prop={'size':5})

base.show_plot_max()
NH3 = 17.03026
## extract TRF 26.912

H2O = 18.01528
NH3 = 17.03026


data = pd.read_csv(ld + 'Baptiste_TRF.csv',na_values="nan",header=None)

array = np.loadtxt(ld + 'Baptiste_TRF.csv', delimiter=',')

T = np.arange(180,340+1,1)
molfrac = np.arange(0,1.01,.01)
mol_in_g = molfrac*NH3+(1-molfrac)*H2O
wt = molfrac*NH3/mol_in_g
interp = RegularGridInterpolator((T, wt), array,bounds_error=False,fill_value=np.NaN)

T_want = np.arange(196,340,0.5)
for i, wt in enumerate(wt_ls):
    wt_want = np.repeat(wt/100,len(T_want))
    pts = np.dstack((T_want, wt_want))
    out = interp((T_want,wt_want))/1000

    out2 = pd.DataFrame(out,index=T_want)
    out2.to_csv(ld +f'/BapOut/{wt}.csv')

## plot energy decomp

fig, ax = plt.subplots(1,2,sharey=True)
idx = 1
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(
    rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\z_paper_figures\{wt}wt%_hb_{m}g.csv',
    header=0)

ax[0].stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])

ax[0].set_ylim([0, 1])
ax[0].set_xlim(184.5,233.5)

ax[0].text(0.5, 0.08, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.5, 0.38, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.5, 0.78, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes)
ax[0].text(0.72, 0.657, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax[0].transAxes)
ax[0].arrow(220,0.655,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax[0].text(0.05, 0.94, '8.2 wt%', fontsize=4, horizontalalignment='left',
     verticalalignment='center', transform=ax[0].transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))




# handles, labels = ax[0].get_legend_handles_labels()
# legend = ax[0].legend(handles[::-1], labels[::-1],fontsize=4, loc='upper left',frameon=True,fancybox=False,framealpha=1,edgecolor='k')
# legend.get_frame().set_linewidth(.3)

idx = 5
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(
    rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\z_paper_figures\{wt}wt%_hb_{m}g.csv',
    header=0)

ax[1].stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])
ax[1].set_ylim([0, 1])
ax[1].set_xlim(184.5,227.5)

ax[1].text(0.5, 0.11, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.305, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.5, 0.65, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes)
ax[1].text(0.72, 0.453, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax[1].transAxes)
ax[1].arrow(216,0.45,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax[1].text(0.05, 0.94, '20.07 wt%', fontsize=4, horizontalalignment='left',
     verticalalignment='center', transform=ax[1].transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))


base.splt_axis_label(fig,'Temperature (K)','Heat (Normalized)')
plt.tight_layout()

# base.show_plot_max()
saveFig = rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\z_paper_figures\energy_budget.png'
plt.savefig(saveFig)
plt.close()



## plot thermograms

fig = plt.figure()
# plt.subplot(1,1,1)
idx = 1
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Te  chnological University\OFYP\NH3_Calorimetry\i_data_processed\03Jan\{wt}wt%_hf_prepped_{m}g.csv',
        header=0)
plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'{wt} wt%')
idx = 5
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data_processed\03Jan\{wt}wt%_hf_prepped_{m}g.csv',
        header=0)
plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'{wt} wt%')


plt.annotate('Glass',(151.5,-20),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('AMH',(171,0),xytext=(0,10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)
plt.annotate('ADH',(177.5,0),xytext=(0,20),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)
plt.annotate('Liquidus 20.07 wt%',(235,-48),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Liquidus 8.2 wt%',(263,-132),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Artifact',(296,-36),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Endothermic',(145,-240),xytext=(0,20),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=3,headlength=3,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)

plt.annotate('Partial Melting Domain',(185,5),ha='left',va='top',fontsize=4,rotation=-3)
plt.annotate('Liquid Domain',(273,-7),ha='left',va='top',fontsize=4)


plt.xlim([130,313])
plt.legend(prop={'size': 5})
base.splt_axis_label(fig,'Temperature (K)','Heat Flow (mW)')


# plt.tight_layout()
# base.show_plot_max()
plt.savefig(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\z_paper_figures\thermogram.png')
## plot thermogram with liquidus double peak

data = pd.read_csv(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data_processed\03Jan\8.4wt%_hf_prepped_4.5858g.csv',header=0)

fig = plt.figure()
plt.subplot(2,1,1)
plt.plot(data['sampleT(K)'],data['Q_corrected(mW)'],color=colour_ls[2],label='8.4 wt%')
plt.legend(prop={'size': 5})
plt.subplot(2,1,2)
for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    m = mass_ls[idx]
    c = colour_ls[idx]
    data = pd.read_csv(
        rf'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data_processed\03Jan\{wt}wt%_hf_prepped_{m}g.csv',
        header=0)
    plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'{wt} wt%')
    plt.xlim([174,184])
    plt.ylim([-250,0])
plt.tight_layout()
plt.legend(prop={'size': 5})

# base.show_plot_max()
base.splt_axis_label(fig,'T (K)','Heat Flow (mW)')
plt.savefig(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\z_paper_figures\doublepeak.png')
##
m_max = np.zeros([len(wt_ls), 1]); m_max[:] = np.nan
p_max = np.zeros([len(wt_ls), 1]); p_max[:] = np.nan
for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    m = mass_ls[idx]

    # FIND MAX ERROR
    data = base.DataCP(m,wt)
    data.import_data(dd, mean=0)
    m_max[idx] = data.df_m['cp_err(%)'].max()
    p_max[idx] = data.df_p['cp_err(%)'].max()
print('ye')


##

def tanh_eqn(T, a, b, c,d):
    """tanh equation for fitting"""
    return a*np.tanh((T+b)*c)+d

data =pd.read_csv(r"C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data\others\liquidus.csv", header=0)

out, pcov = curve_fit(tanh_eqn,data['Theoretical'],data['Offset'],p0=[.3,-250,.015,.3])

T = np.arange(100,500)
plt.scatter(data['Theoretical'],data['Offset'])
# plt.plot(T,.3*np.tanh((T-250)*.015)+.3)
plt.plot(T,tanh_eqn(T,*out),'k')
plt.title('tanh fit to phase transition offsets')
plt.xlabel('Theoretical T of Transition (K)')
plt.ylabel('Measured T - Theoretical T (K)')
plt.savefig(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\NH3_Calorimetry\i_data\others\liquidus.png')
base.show_plot_max()

## debug end line
print('Finished running script')