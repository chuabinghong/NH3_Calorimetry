import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import base

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'


## FILE DIRECTORY SETTINGS


wd = r"../"
os.chdir(wd)


## PREPARE HEATFLOW FILE (WITHOUT MEAN SMOOTHING)

## FILE DIRECTORY SETTINGS

dd = r"i_data_processed/not_smoothed/"
fd = r"o_supplementaryPlots/"

wt_ls = [8.2, 20.07] # based off liquidus alignment - 0.5K
mass_ls = [4.1943,3.7107]

colour_ls = ['#f46d43','#74add1']


## plot thermograms

w= 2
idx = 0
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
# data = pd.read_csv(dd+
#         rf'{wt}wt%_hf_prepped_{m}g_not_smoothed.csv',
#         header=0)

# fig,ax = plt.subplots(1,figsize=(w,w))
# plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'8.2 wt%')
#
#
# plt.annotate('Eutectic',(179,-153),xytext=(20,-10),textcoords='offset points',
#              arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=6)
# plt.annotate('Liquidus',(263,-125),xytext=(0,-10),textcoords='offset points',
#              arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=6)
#
# # plt.annotate('Partial Melting Domain',(185,5),ha='left',va='top',fontsize=4,rotation=-3)
# # plt.annotate('Liquid Domain',(273,-7),ha='left',va='top',fontsize=4)
#
# # rect = matplotlib.patches.Rectangle((180, -155), 40, 30, linewidth=1, edgecolor='k',facecolor='none')
# # ax.add_patch(rect)
#
# plt.xlim([173,313])
# plt.ylim([-180,0])
# plt.title('Thermogram from calorimeter')
# # plt.legend(prop={'size': 5})
# base.show_plot_max()
# base.splt_axis_label(fig,'Temperature (K)','Heat Flow (mW)')
#
#
# # plt.tight_layout()
#
# plt.savefig(fd+ r'GA_thermogram.png')


##
dd = r"i_data_processed/"
data = pd.read_csv(dd+ rf'{wt}wt%_hb_{m}g.csv', header=0)


fig,ax = plt.subplots(1,figsize=(w,w))

ax.stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])

ax.set_ylim([0, 1])
ax.set_xlim(184.5,233.5)

ax.text(0.5, 0.10, 'melt ice', fontsize=6, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.5, 0.41, 'heat ice fraction', fontsize=6, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.5, 0.81, 'heat liquid fraction', fontsize=6, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.72, 0.687, 'mix molten ice', fontsize=6, horizontalalignment='center',
     verticalalignment='bottom', transform=ax.transAxes)
ax.arrow(220,0.687,0,-.05,head_width=.5,head_length=.005,linewidth=.7)
plt.title('Heat Budget in Partial Melting Domain')
base.splt_axis_label(fig,'Temperature (K)','Heat (Normalized)')
base.show_plot_max()
plt.savefig(fd+ r'GA_heatbudget.png')