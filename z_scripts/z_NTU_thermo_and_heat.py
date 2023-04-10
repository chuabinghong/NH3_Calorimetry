import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import base
import matplotlib.gridspec as gridspec

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'


## FILE DIRECTORY SETTINGS


wd = r"../"
os.chdir(wd)


## PREPARE HEATFLOW FILE (WITHOUT MEAN SMOOTHING)

# peak_arg = 0  # 1 if we are getting peaks (for liquidus point identification) (broken in current implementation (29 Mar 23)
# range_arg = 1  # 1 if using -196C to 46C
# ramp_rate = 0.1
# T_bin = 3  # size of T_averaged bin
#
#
# sd = r"i_data_processed/not_smoothed/"
# fd = r"o_heatFlow/"
#
# if range_arg == 1:
#     dd = r"i_data/46C/"
#     mass_ls = [4.1943, 3.7107]
#     wt_ls = [8.2, 20.07]  # based off liquidus alignment
#     t_end = 13066
#     lb_ls = [289, 289]  # lower and upper bound of kink to be removed
#     ub_ls = [300, 298]
#
# # import blank
# blank = base.DataRaw(0, 'blank', ramp_rate, t_end)
# blank.df = pd.read_csv(dd + 'Blank.csv', skiprows=8)
# blank.correct_HF(blank.df)
#
# for idx, wt in enumerate(wt_ls):
#     # idx = 0 #TODO for single data debugging
#     # wt = wt_ls[idx] #TODO for single data debugging
#
#     m = mass_ls[idx]
#
#     # data processing
#     data = base.DataRaw(m, wt, ramp_rate, t_end)
#     data.import_raw(dd)
#     if range_arg == 1:
#         lb = lb_ls[idx]
#         ub = ub_ls[idx]
#         data.remove_kink(lb, ub)
#     data.correct_HF(blank.df) # calibration with blank
#
#     saveFile = sd + '{}wt%_hf_prepped_{}g_not_smoothed.csv'.format(data.wt, data.m)
#     col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)']
#     data.df2.to_csv(saveFile, index=False, header=col_names)

## FILE DIRECTORY SETTINGS

dd = r"i_data_processed/not_smoothed/"
fd = r"o_finalPlots/"

wt_ls = [8.2, 20.07] # based off liquidus alignment - 0.5K
mass_ls = [4.1943,3.7107]

colour_ls = ['#f46d43','#74add1']


## plot thermograms


# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)

fig = plt.figure()

##
ax = plt.subplot(gs[0, :]) # row 0, col 0

idx = 0
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(dd+
        rf'{wt}wt%_hf_prepped_{m}g_not_smoothed.csv',
        header=0)
plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'{wt} wt%')
idx = 1
wt = wt_ls[idx]
m = mass_ls[idx]
c = colour_ls[idx]
data = pd.read_csv(dd+
        rf'{wt}wt%_hf_prepped_{m}g_not_smoothed.csv',
        header=0)
plt.plot(data['sampleT(K)'], data['Q_corrected(mW)'], color=c,zorder=10-idx,label=f'20.1 wt%')


plt.annotate('Glass',(151.5,-20),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('AMH',(171,0),xytext=(0,10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)
plt.annotate('ADH',(177.5,0),xytext=(0,20),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)
plt.annotate('Liquidus 20.1 wt%',(235,-48),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Liquidus 8.2 wt%',(263,-132),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Artifact',(296,-36),xytext=(0,-10),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=2,headlength=2,edgecolor=None,linewidth=.5),ha='center',va='top',fontsize=4)
plt.annotate('Endothermic',(145,-240),xytext=(0,20),textcoords='offset points',
             arrowprops = dict(facecolor ='k',width=.5,headwidth=3,headlength=3,edgecolor=None,linewidth=.5),ha='center',va='bottom',fontsize=4)

plt.annotate('Partial Melting Domain',(185,11),ha='left',va='top',fontsize=4,rotation=-3)
plt.annotate('Liquid Domain',(273,-7),ha='left',va='top',fontsize=4)


plt.xlim([130,313])
plt.legend(prop={'size': 5})
plt.ylabel('Heat Flow (mW)')

plt.text(0.025, 0.98, 'a)', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

##
ld = r"i_data_literature/"
dd = r"i_data_processed/"
fd = r"o_finalPlots/"

ax2 = plt.subplot(gs[1, 0]) # row 0, col 0


idx = 0
wt = wt_ls[idx]
m = mass_ls[idx]
data = pd.read_csv(dd+ rf'{wt}wt%_hb_{m}g.csv', header=0)

ax2.stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])

ax2.set_ylim([0, 1])
ax2.set_xlim(184.5,233.5)

ax2.text(0.5, 0.10, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax2.transAxes)
ax2.text(0.5, 0.41, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax2.transAxes)
ax2.text(0.5, 0.81, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax2.transAxes)
ax2.text(0.72, 0.687, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax2.transAxes)
ax2.arrow(220,0.685,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax2.text(0.95, 0.94, '8.2 wt%', fontsize=4, horizontalalignment='right',
     verticalalignment='center', transform=ax2.transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))
ax2.set_ylabel('Heat (Normalized)')
ax2.text(0.025, 0.98, 'b)', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)



##

ax = plt.subplot(gs[1, 1], sharex = ax2) # row 0, col 0

idx = 1
wt = wt_ls[idx]
m = mass_ls[idx]
data = pd.read_csv(dd+ rf'{wt}wt%_hb_{m}g.csv', header=0)

ax.stackplot(data['T(K)'], data['c_melt'] / data['j_total'],
              data['c_heat'] / data['j_total'],
              data['l_mix'] / data['j_total'],
              data['l_heat'] / data['j_total'],
              labels=['melt ice', 'heat ice', 'heat liquid'], colors = ['#abd9e9','#4575b4','#e8e884','#f46d43'])
ax.set_ylim([0, 1])
ax.set_xlim(184.5,227.5)

ax.text(0.5, 0.145, 'melt ice', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.5, 0.352, 'heat ice fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.5, 0.68, 'heat liquid fraction', fontsize=4, horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes)
ax.text(0.72, 0.487, 'mix molten ice', fontsize=4, horizontalalignment='center',
     verticalalignment='bottom', transform=ax.transAxes)
ax.arrow(216,0.485,0,-.05,head_width=.5,head_length=.005,linewidth=.5)

ax.text(0.95, 0.94, '20.1 wt%', fontsize=4, horizontalalignment='right',
     verticalalignment='center', transform=ax.transAxes,bbox=dict(boxstyle='square', facecolor='white', alpha=1,linewidth=0.5))

ax.set_yticklabels([])
ax.text(0.025, 0.98, 'c)', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)


base.splt_axis_label(fig,'Temperature (K)','')
##
w = 3.3
fig.set_size_inches(w, w * 1)

plt.tight_layout()
# base.show_plot_max()
plt.savefig(fd+ r'NTU_thermograms.png')
