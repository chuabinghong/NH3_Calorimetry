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

wd = r"../../"
os.chdir(wd)

dd = r"i_data/water_tests/"
sd = r"i_data/water_tests/"
fd = r"o_heatFlow/"

## SCRIPT ARGUMENTS
peak_arg = 0 #1 if we are getting peaks
range_arg = 0 #1 if using -196C to 46C
ramp_rate = 0.1
T_bin = 3 #size of T_averaged bin

mass_ls = [0.190,0.340,6.000]
wt_ls = [0, 0, 0]  # based off liquidus alignment offset
t_end = 13066

mass_ls = [6.00000001]
wt_ls = [0]  # based off liquidus alignment offset
t_end = 13066

## DATA PROCESSING

# import blank
blank = base.DataRaw(0,'blank',ramp_rate,t_end)
blank.df = pd.read_csv(dd + '2014_blank.csv', skiprows=8)
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
    data.mean_smoothing(data.df2,T_bin,sd)

    data.correct_HF_smooth(blank.df_mean) # calibration with blank

    data.plot_Q_T(fd,plt_err=1,save=1)
    data.calc_Fc() # calculating crystal fraction
    # data.calc_liquidus()
    data.plot_Fc_T(fd,1)
    
    data.save_data(sd)
    # data.save_liquidus(sd)


    ## FIND PEAKS (MANUAL FOR EACH RUN)

    if peak_arg == 1:

        # plt.close()
        # spl = UnivariateSpline(data.df2['Sample Temperature(K)'][300:-1],data.df2['Q_Corrected'][300:-1],k=5)
        # # spl = UnivariateSpline(data.df2['Sample Temperature(K)'][3700:4200],data.df2['Q_Corrected'][3700:4200],k=2)
        # plt.plot(data.df2['Sample Temperature(K)'][3700:4200], data.df2['Q_Corrected'][3700:4200], label=r'Constant $\Delta$T')
        # plt.plot(data.df2['Sample Temperature(K)'],spl(data.df2['Sample Temperature(K)']))
        # # plt.xlim([200,215])
        # # plt.ylim([-30,-20])
        # plt.show()


        min_id_r = np.ndarray.flatten(np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.less,order=90)])) #find minima peak indexes, plus index offset
        max_id_r = np.ndarray.flatten(np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.greater,order=90)])) #find minima peak indexes, plus index offset


        if smallHF_arg == 1:
            if idx == 0:
                min_id_r = np.ndarray.flatten(
                    np.asarray([x for x in argrelextrema(np.array(data.df2['Q_Corrected']), np.less, order=10)]))
                min_id = [min_id_r[i] for i in [118]]
                max_id = [max_id_r[i] for i in [1]]

        elif range_arg == 1 and smallHF_arg!=1:
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

        elif range_arg != 1 and smallHF_arg!=1:
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

        if range_arg != 1 and smallHF_arg!=1 and (idx==3 or idx==4):
            dfpeak = pd.concat([dfpeak,dfpeak2.to_frame().T],ignore_index=True)


        # --- plot Peaks --- #
        plt.plot(data.df2['Sample Temperature(K)'], data.df2['Q_Corrected'], 'g-', label=r'Constant $\Delta$T')
        plt.title('{} wt% Thermogram Peaks'.format(data.wt))
        plt.xlabel("Temperature (K)")
        plt.ylabel("Heat Flow (mW)")
        plt.plot(dfpeak["Sample Temperature(K)"], dfpeak["HeatFlow(mW)"], 'kx',markersize=3)
        if range_arg != 1 and smallHF_arg!=1 and (idx==3 or idx==4):
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

## EXTRAS

# PLOT BLANK HF COMPARISON
# blank1 = base.DataRaw(0,'blank',ramp_rate,12981)
# blank1.df = pd.read_csv(r"i_data/0.10Kmin/" + 'Blank.csv', skiprows=8)
# blank1.correct_HF(blank.df)
#
# blank2 = base.DataRaw(0,'blank',ramp_rate,13066)
# blank2.df = pd.read_csv(r"i_data/46C/" + 'Blank.csv', skiprows=8)
# blank2.correct_HF(blank.df)
#
# blank3 = base.DataRaw(0,'blank',ramp_rate,13066)
# blank3.df = pd.read_csv(r"i_data/smallHF/" + 'Blank.csv', skiprows=8)
# blank3.correct_HF(blank.df)
#
# plt.figure()
# plt.plot(blank1.df2['Sample Temperature(K)'], blank1.df2['Q_Corrected'],label='26C blank')
# plt.plot(blank2.df2['Sample Temperature(K)'], blank2.df2['Q_Corrected'],label='46C blank')
# plt.plot(blank3.df2['Sample Temperature(K)'], blank3.df2['Q_Corrected'],label='46C blank small HF')
# plt.xlabel("Temperature (K)");
# plt.ylabel('Heat Flow (mW)')
# plt.title('Blank Heatflow')
# plt.legend()
# plt.ylim([-2, 2])
# plt.savefig(fd + 'blank_Q-T.png')
# plt.close()
# ##
#
# a = interp1d(blank2.df2['Sample Temperature(K)'],blank2.df2['Q_Corrected'],bounds_error=False,fill_value=0)
# b = interp1d(blank1.df2['Sample Temperature(K)'],blank1.df2['Q_Corrected'],bounds_error=False,fill_value=0)
#
# dev = a(blank1.df2['Sample Temperature(K)']) - b(blank1.df2['Sample Temperature(K)'])
# plt.plot(blank1.df2['Sample Temperature(K)'],dev)
# plt.title("46C blank - 26C blank")
# plt.savefig(fd + 'blank_dev.png')

# # plot heatflow of blank wrt T
# plt.figure()
# plt.plot(blank.df["Furnace Temperature(K)"][blank.t_start:blank.t_end], blank.df["HeatFlow(mW)"][blank.t_start:blank.t_end],
#          'r--', label='Furnance')
# plt.plot(blank.df["Sample Temperature(K)"][blank.t_start:blank.t_end], blank.df["HeatFlow(mW)"][blank.t_start:blank.t_end],
#          'c-', label='Sample')
# plt.xlabel('T (K)')
# plt.ylabel('Heat Flow (mW)')
# plt.title('Heat Flows of Blank')
# plt.legend()
# # plt.show()
# saveFig = fd + 'blank_heatflow.png'
# plt.savefig(saveFig)
# plt.close()
#
# # plot T wrt t of sample and blank
# plt.figure()
# plt.plot(data.df["Time(s)"], data.df["Sample Temperature(K)"], 'c-', label='Sample')
# plt.plot(blank.df["Time(s)"], blank.df["Sample Temperature(K)"], 'r', label='Blank')
# plt.xlabel('t (s)')
# plt.ylabel('T (K)')
# plt.title('{}wt%'.format(data.wt))
# plt.legend()
# # plt.show()
# saveFig = fd + '{}wt%_T-t.png'.format(data.wt)
# plt.savefig(saveFig)
# plt.close()
