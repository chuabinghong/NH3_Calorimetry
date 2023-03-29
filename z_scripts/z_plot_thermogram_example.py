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

fig = plt.figure()
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
plt.savefig(fd+ r'thermograms.png')
