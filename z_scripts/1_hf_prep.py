"""
created by Bing Hong CHUA 27Sep22

script objective:
from raw data obtained from the calorimeter
output corrected heatflow data to be used in further analysis
plot thermograms, crystal fraction curves and identify peaks
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import base
from scipy.signal import argrelextrema
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

sd = r"i_data_processed/"
fd = r"o_supplementaryPlots/"

## SCRIPT ARGUMENTS
peak_arg = 0 #1 if we are getting peaks (for liquidus point identification) (broken in current implementation (29 Mar 23)
range_arg = 0 #1 if using -196C to 46C
ramp_rate = 0.1
T_bin = 3 #size of T_averaged bin

# adjust script arguments according to range_arg
if range_arg == 0:
    dd = r"i_data/26C/"
    mass_ls = [4.5858, 4.5202, 3.7778]
    wt_ls = [8.4, 10.0, 26.912]  # based off liquidus alignment offset
    t_end = 12981

elif range_arg == 1:
    dd = r"i_data/46C/"
    mass_ls = [4.5386,4.1943,3.8153,3.7107]
    wt_ls = [5.2, 8.2, 14.3, 20.07]  # based off liquidus alignment
    t_end = 13066
    lb_ls = [289, 289, 293, 289] # lower and upper bound of kink to be removed
    ub_ls = [301.5, 300, 305, 298]

## DATA PROCESSING

# import blank
blank = base.DataRaw(0,'blank',ramp_rate,t_end)
blank.df = pd.read_csv(dd + 'Blank.csv', skiprows=8)
blank.t_start = blank.df['Sample Temperature(K)'].idxmin()
blank.df2 = blank.df[blank.t_start:blank.t_end].copy()
blank.mean_smoothing(blank.df2,T_bin,sd)
blank.df2 = blank.df_mean.copy()
blank.df2['Q_Corrected'] = blank.df2['HeatFlow(mW)']
blank.plot_Q_T(fd,plt_err=1,save=1)

for idx, wt in enumerate(wt_ls):
    # idx = 0 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging

    m = mass_ls[idx]

    # data processing
    data = base.DataRaw(m, wt,ramp_rate,t_end)
    data.import_raw(dd)
    if range_arg ==1:
        lb = lb_ls[idx]
        ub = ub_ls[idx]
        data.remove_kink(lb,ub)
    data.mean_smoothing(data.df2,T_bin,sd)

    data.correct_HF_smooth(blank.df_mean) # calibration with blank

    data.plot_Q_T(fd,plt_err=1,save=1)
    data.calc_Fc() # calculating crystal fraction
    # data.calc_liquidus()
    data.plot_Fc_T(fd,1)
    
    data.save_data(sd)
    # data.save_liquidus(sd)


    ## FIND PEAKS (MANUALLY SELECTED FOR EACH RUN)

    if peak_arg == 1:

        min_id_r = np.ndarray.flatten(np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.less,order=90)])) #find minima peak indexes, plus index offset
        max_id_r = np.ndarray.flatten(np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.greater,order=90)])) #find minima peak indexes, plus index offset

        if range_arg == 1:
            if idx == 0:
                min_id = [min_id_r[i] for i in [3,5]]
                max_id = [max_id_r[i] for i in [1]]
            elif idx == 1:
                min_id = [min_id_r[i] for i in [3,4,5]]
                max_id = [max_id_r[i] for i in [2]]
            elif idx == 2:
                min_id = [min_id_r[i] for i in [3, 4]]
                max_id = [max_id_r[i] for i in [1]]
            elif idx == 3:
                min_id = [min_id_r[i] for i in [3,5]]
                max_id = [max_id_r[i] for i in [1]]
            elif idx == 4:
                min_id = [min_id_r[i] for i in [3, 4]]
                max_id = [max_id_r[i] for i in [2]]
            elif idx == 5:
                min_id_r = np.ndarray.flatten(
                    np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.less, order=10)]))
                min_id = [min_id_r[i] for i in [41]]
                max_id = [max_id_r[i] for i in [1]]

        elif range_arg != 1:
            if idx==0:
                min_id = np.ndarray.flatten(
                    np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.less, order=1000)]))
                max_id = min_id
            elif idx == 1:
                min_id = [min_id_r[i] for i in [2, 3, 6]]
                max_id = [max_id_r[i] for i in [2,4]]
            elif idx == 2:
                min_id = [min_id_r[i] for i in [2, 3, 5]]
                max_id = [max_id_r[i] for i in [1, 2]]
            elif idx == 3:
                min_id = [min_id_r[i] for i in [3, 4, 5]]
                max_id = [max_id_r[i] for i in [3, 4]]
                spl = UnivariateSpline(data.df2['Sample Temperature(K)'][300:-1],data.df2['Q_Corrected'][300:-1],k=4)
                min_2 = np.ndarray.flatten(np.asarray([x for x in argrelextrema(np.array(spl(data.df2['Sample Temperature(K)'])), np.less,order=200)]))
                dfpeak2 = data.df2.iloc[min_2[2]]
            elif idx==4:
                min_id = [min_id_r[i] for i in [3,4,7]]
                max_id = [max_id_r[i] for i in [2]]
                spl = UnivariateSpline(data.df2['Sample Temperature(K)'][300:-1], data.df2['Q_Corrected'][300:-1], k=4)
                min_2 = np.ndarray.flatten(np.asarray(
                    [x for x in argrelextrema(np.array(spl(data.df2['Sample Temperature(K)'])), np.less, order=200)]))
                dfpeak2 = data.df2.iloc[min_2[3]]
            elif idx == 5:
                min_id = [min_id_r[i] for i in [2, 3, 4]]
                max_id = [max_id_r[i] for i in [3]]
            elif idx==6:
                max_id_r = np.ndarray.flatten(
                    np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.greater, order=25)]))
                min_id = [min_id_r[i] for i in [2,3,4]]
                max_id = [max_id_r[i] for i in [15]]
            else:
                min_id = min_id_r
                max_id = max_id_r

        peak_id = np.concatenate((min_id, max_id), axis=None)
        dfpeak = data.df2.iloc[peak_id]
        dfpeak = dfpeak.sort_index()

        if range_arg != 1 and (idx==3 or idx==4):
            dfpeak = pd.concat([dfpeak,dfpeak2.to_frame().T],ignore_index=True)


        # --- plot Peaks --- #
        plt.plot(data.df2['Sample Temperature(K)'], data.df2['Q_Corrected'], 'g-', label=r'Constant $\Delta$T')
        plt.title('{} wt% Thermogram Peaks'.format(data.wt))
        plt.xlabel("Temperature (K)")
        plt.ylabel("Heat Flow (mW)")
        plt.plot(dfpeak["Sample Temperature(K)"], dfpeak["HeatFlow(mW)"], 'kx',markersize=3)
        if range_arg != 1 and (idx==3 or idx==4):
            plt.plot(dfpeak2["Sample Temperature(K)"], dfpeak2["HeatFlow(mW)"], 'rx',markersize=3)
        # plt.xlim([230,240])
        # plt.ylim([-50,-20])
        # base.show_plot_max()
        saveFig = fd + '{}wt%_Q-T_peaks.png'.format(data.wt)
        plt.savefig(saveFig)
        plt.close()

        saveFile = sd + '{}wt%_peaks.txt'.format(data.wt)
        with open(saveFile, 'w') as fileOut:
            lines = dfpeak.to_string(index=False, header=True)
            fileOut.writelines(lines)

    print('Finished loop {}/{}'.format(idx+1,len(wt_ls)))

## debug end line

print('Finished running script')