"""
created by Bing Hong CHUA 08Dec22

script objective:
identify oscillation of cp from shomate fit
"""

import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import base
from scipy.interpolate import UnivariateSpline

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
fd = r"t_all_cp/"

## SCRIPT ARGUMENTS
range_arg = 1
melt_arg = 0
pure_arg = 1
ramp_rate = 0.1
t_range= [9000,10300]
# wt_ls = [5.2, 5.2, 8.2, 8.4, 10.0, 14.4, 20.1, 26.93] # based off liquidus alignment - 0.5K
# mass_ls = [4.6293,4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]
colour_ls = ['tab:blue','tab:green','tab:orange','tab:red','tab:purple','tab:pink','tab:olive','dimgray']

mass_ls = [4.6293,4.5858,4.5202,3.7778]
wt_ls= [ 5.2, 8.4, 10.0, 26.93]
t_end = 12981
dd2 = r"i_data/0.10Kmin/"


# adjust script arguments according to range_arg
if range_arg == 1:
    dd2 = r"i_data/46C/"
    sd = r"i_data_processed/"
    fd = r"o_heatFlow/"
    mass_ls = [4.5386,4.1943,3.8153,3.7107]
    wt_ls = [5.2, 8.2, 14.4, 20.1]  # based off liquidus alignment
    t_end = 13066
if range_arg == 2:
    dd2 = r"i_data/0.25Kmin/"
    mass_ls = [4.9859,4.8740,4.7648,4.7459,4.5983,4.5500,4.3987]
    wt_ls = [0,2,4,4,6,8,10]
    t_end = 13066
    ramp_rate = 0.25
    t_range = [10000,11000]

## CHECK T RAMP SIMILARITY
# ---- T RAMP ALIGNS WELL ----

# BLANK
blank = base.DataRaw(0,'blank',ramp_rate,t_end)
blank.df = pd.read_csv(dd2 + 'Blank.csv', skiprows=8)
blank.correct_HF(blank.df)
plt.plot(blank.df2['Time(s)'],blank.df2['Sample Temperature(K)'])
base.show_plot_max()

figall, axsall = plt.subplots(len(wt_ls),1,sharex=True)
for idx, wt in enumerate(wt_ls):
    # idx = 5 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging
    colour = colour_ls[idx]
    m = mass_ls[idx]

    data = base.DataRaw(m,wt,ramp_rate,t_end)
    data.import_raw(dd2)
    data.dfcut = data.df.iloc[t_range[0]:t_range[1]]
    # data.mean_sampling(data.dfcut,3,od)
    s = UnivariateSpline(data.dfcut['Time(s)'],data.dfcut['HeatFlow(mW)'],s=10)

    dev = data.dfcut['HeatFlow(mW)'] - s(data.dfcut['Time(s)'])

    plt.figure()
    plt.plot(data.dfcut['Time(s)'],data.dfcut['HeatFlow(mW)'],label='data')
    plt.plot(data.dfcut['Time(s)'],s(data.dfcut['Time(s)']),linewidth=0.5,label='spline')
    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('heatflow (mW)')
    plt.savefig(r't_fft/' + f'{wt}wt%_dev_{m}g.png')

    cutoff = .05
    dt= data.df['Time(s)'][1] - data.df['Time(s)'][0]
    # f = data.dfcut['HeatFlow(mW)']
    f = dev
    t = data.dfcut['Time(s)']

    n = len(f)
    fhat = np.fft.fft(f,len(f))
    PSD = fhat * np.conj(fhat) / n
    freq = (1/(dt*n)) * np.arange(n)
    L = np.arange(1,np.floor(n/2),dtype='int')

    indices = PSD < cutoff
    PSDclean = PSD * indices
    fhatclean = indices * fhat
    ffilt = np.fft.ifft(fhatclean)

    fig,axs = plt.subplots(2,1)
    # plt.title('FFT of deviations')
    plt.sca(axs[0])
    plt.plot(t,f)
    plt.plot(t, ffilt,label='inverse fft')
    axs[0].set(xlabel='time (s)',ylabel='data - spline')

    plt.sca(axs[1])
    plt.plot(freq[L],PSD[L])
    plt.axhline(y=cutoff, color='k')
    plt.xlim(freq[L[0]],.005)
    axs[1].set(xlabel='frequency (Hz)', ylabel='power spectrum density')
    plt.savefig(r't_fft/' + f'{wt}wt%_fft_{m}g_fit.png')

    plt.sca(axsall[idx])
    plt.plot(freq[L],PSD[L],label=f'{wt}wt%')
    plt.legend()
    plt.xlim(freq[L[0]],0.005) #TODO hardcoded upper lim

# plt.savefig(r't_fft/' + f'fft_all_rangearg={range_arg}.png')
    # plt.sca(axs[1])
    # plt.plot(t,f)
    # plt.plot(t,ffilt)
    #
    # plt.sca(axs[2])
    # plt.plot(freq[L],PSD[L])
    # plt.axhline(y=cutoff,color='k')
    # plt.xlim(freq[L[0]],freq[L[-1]])

    # base.show_plot_max()
print('ye')

## CHECK BLANKS
#

T_bin = 10
# BLANK
blank = base.DataRaw(0,'blank',ramp_rate,t_end)
blank.df = pd.read_csv(dd2 + 'Blank.csv', skiprows=8)
blank.correct_HF(blank.df)
blank.mean_sampling(blank.df2,T_bin,od)
#
# # s = UnivariateSpline(blank.df_mean['Sample Temperature(K)'],blank.df_mean['Q_Corrected'], s=1)
#
plt.plot(blank.df2['Sample Temperature(K)'],blank.df2['Q_Corrected'],label='data')
plt.plot(blank.df_mean['Sample Temperature(K)'],blank.df_mean['Q_Corrected'],label=f'{T_bin} K bins')
plt.legend()
plt.savefig(r't_fft/' + f'blank_mean_{T_bin}K.png')
# plt.plot(blank.df_mean['Sample Temperature(K)'],s(blank.df_mean['Sample Temperature(K)']))
#
base.show_plot_max()
# ## CHECK CP FIT DEVIATIONS
# for idx, wt in enumerate(wt_ls):
#     # idx = 5 #TODO for single data debugging
#     # wt = wt_ls[idx] #TODO for single data debugging
#     colour = colour_ls[idx]
#     m = mass_ls[idx]
#
#     data = base.DataCP(m,wt)
#     data.import_data(dd,mean=1)
#
#     data.calc_shomate('mean')
#     dev = 20*(data.df_m_mean['cp(J/gK)'] - base.shomate_eqn(data.df_m_mean['T(K)'], *data.shomate_m))
#
#
#
#     plt.plot(data.df_m_mean['T(K)'],dev)
# plt.axhline(0,color='k')
# base.show_plot_max()
