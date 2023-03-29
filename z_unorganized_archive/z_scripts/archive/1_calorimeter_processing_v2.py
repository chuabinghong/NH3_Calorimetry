"""
created by Jack DIAB 26Jun22
modified by Bing Hong CHUA

script objective:
exploratory analysis of calorimeter output data
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
plt.rcParams.update({'figure.dpi': '600'})

##

old_arg = 0 #1 if we are testing 0.25Kmin
peak_arg = 1 #1 if we are getting peaks
old_wt_arg = 0 #1 if we are testing old wt% measurements

# set working and save directories if required
dd = r"../i_data/0.10Kmin/"
sd = r"../i_data_processed/"
fd = r"../o_heatFlow/"

# --- hardcode start and end t --- #
start_offset = 6000  # FOR 0.1KMIN
end_offset = 12900

# data files
wtpct_ls = [0,1.42,2.84,5.11,6.88,8.73,'newBottle']
mass_ls_in = [2510.8, 4887.3,4840.2,4629.3,4585.8,4520.2,3777.8]
wtpct_ls = [2.84]
mass_ls = [4520.2]
file_ls = [f'H2O-NH3_{i}wt%.csv' for i in wtpct_ls]
mass_ls = [x*0.001 for x in mass_ls_in]

if old_wt_arg ==1:
    dd = r"../i_data/0.10Kmin/old_wt%/"
    sd = r"../i_data_processed/old_wt%/"
    fd = r"../o_heatFlow/old_wt%/"
    wtpct_ls = [0, 1.62, 3.24, 5.83, 7.85, 'newBottle']
    file_ls = [f'H2O-NH3_{i}wt%.csv' for i in wtpct_ls]
    mass_ls = [x * 0.001 for x in mass_ls_in]

# set working and save directories for 0.25Kmin
if old_arg==1:
    dd = r"../i_data/0.25Kmin/"
    sd = r"../i_data_processed/0.25Kmin/"
    fd = r"../o_heatFlow/0.25Kmin/"
    start_offset = 8500 #FOR 0.25KMIN
    end_offset = 12100
    wtpct_ls = [0,1.42,2.84,5.11,6.88,10]
    file_ls = ["H2O-NH3_0wt%.csv","H2O-NH3_2wt%.csv", "H2O-NH3_4wt%.csv", "H2O-NH3_6wt%.csv", "H2O-NH3_8wt%.csv", "H2O-NH3_10wt%.csv"]
    mass_ls = [4.9859,4.8389, 4.7459, 4.5983, 4.5500, 4.3987]


## FUNCTIONS
def round_dec(x, factor=1000000):
    """Takes in number to be rounded and factor to determine decimal place (default = 1000000).
    e.g. A factor of 1000 would round to the thousanths place.
    Returns rounded number"""
    x1 = x * factor
    x2 = round(x1)
    x3 = x2 / factor
    return x3

def cp_calc(mW, mass, dTdt):
    """Calculates specific heat capacity
    inputs: heat flow (mW), mass (g), delta_temperature over sampling rate (K/s)
    outputs: specific heat"""
    # Cp(T) =
    # (dq/dt)/(m × (dT/dt))
    # where m is mass
    # mW to W, then W/(g * K), 1 W = 1 J/s, ==> J/(s * g * K)
    Cp = np.abs((mW * 0.001) / (mass * dTdt))
    return Cp

def Cp_T(T, a, b, c, d, e):
    """Shomate Formulation equation for fitting"""
    return a + (b * T) + (c * T ** 2) + (d * T ** 3) + (e / (T ** 2))

def Xl_calc(T):
    """Croft equation for liquidus of water-rich ammonia-water mixutre, up till eutectic (X<=0.329)
    input: current tempertaure
    output: ammonia weight%"""
    lqd_coeff = [56695,-46269,11842,-1651.4,-53.07,273.15-T] #CROFT 1987
    roots = np.roots(lqd_coeff)
    roots_real = np.real(roots[np.iscomplex(roots)==False])
    root_actual = roots_real[(0<roots_real) & (roots_real<=0.329)]
    if len(root_actual)==0:
        root_actual = np.nan
    else:
        root_actual = root_actual[0]
    return root_actual

## BLANK RUN CHECKS


df_blank = pd.read_csv(dd+"Blank.csv", skiprows=8)

fig = plt.figure(1)
# ax = fig.subplot
plt.plot(df_blank["Furnace Temperature(K)"][start_offset:end_offset], df_blank["HeatFlow(mW)"][start_offset:end_offset], 'r--', label='Furnance')
plt.plot(df_blank["Sample Temperature(K)"][start_offset:end_offset], df_blank["HeatFlow(mW)"][start_offset:end_offset], 'c-', label='Sample')
plt.xlabel('T (K)')
plt.ylabel('Q (mW)')
# plt.title('testing')
plt.legend()
# plt.show()
saveFig = fd + 'blank_Q-T.png'
plt.savefig(saveFig)
plt.close()

## IMPORT DATA FILES

# --- run following for all samples --- #
for samp_id, samp in enumerate(wtpct_ls):
    m = mass_ls[samp_id]
    df = pd.read_csv(dd + file_ls[samp_id], skiprows=8)

    samp = 2.842 #########################################
    m = 0.9642
    df = pd.read_csv(r'../i_data/H2O-NH3_2.84wt%_1ml.csv', skiprows=8)

    if samp == 'newBottle':
        samp = 26.96

    ## T OVER t

    fig = plt.figure(2)
    plt.plot(df["Time(s)"], df["Sample Temperature(K)"], 'c-', label='Sample')
    plt.plot(df_blank["Time(s)"], df_blank["Sample Temperature(K)"], 'r', label='Blank')
    plt.xlabel('t (s)')
    plt.ylabel('T (K)')
    plt.title('{}wt%'.format(samp))
    plt.legend()
    # plt.show()
    saveFig = fd + '{}wt%_T-t.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()


    # --- Rate of T increase --- #
    T_slope = ((df["Sample Temperature(K)"][end_offset]
                - df["Sample Temperature(K)"][start_offset])
               / (df["Time(s)"][end_offset] -
                 df["Time(s)"][start_offset]))
    T_slope = round_dec(T_slope)

    # --- find corrected Q within T climb period --- #
    df['Q_Corrected'] = df["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"]
    dfclimb = df[start_offset:end_offset]

    hf_corrected = df["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"]
    climb_temp = list(df["Sample Temperature(K)"][start_offset:end_offset])
    climb_hf = hf_corrected[start_offset:end_offset]
    # fig = plt.figure(3)
    # plt.plot(df["Sample Temperature(K)"], hf_corrected, 'm-', label='Corrected Sample')
    # plt.axhline(0)

    # --- plot Q-T --- #
    plt.plot(climb_temp, climb_hf, label=r'Corrected')
    # plt.plot(df["Sample Temperature(K)"][start_offset:end_offset],
    #          df["HeatFlow(mW)"][start_offset:end_offset], '--', label='Uncorrected')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat Flow (mW)")
    plt.title(r'{}wt% $NH_3$'.format(samp))
    # plt.legend()
    # plt.show()
    saveFig = fd + '{}wt%_Q-T_2.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    ## FIND PEAKS (MANUAL FOR EACH RUN)

    if peak_arg == 1:

        min_id_r = np.ndarray.flatten(np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.less,order=90)])) #find minima peak indexes, plus index offset
        max_id_r = np.ndarray.flatten(np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.greater,order=90)])) #find minima peak indexes, plus index offset

        if samp_id==0:
            min_id = np.ndarray.flatten(
                np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.less, order=1000)]))
            max_id = min_id
        elif samp_id==4:
            min_id_r = np.ndarray.flatten(
                np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.less, order=31)]))
            min_id = [min_id_r[i] for i in [6,7,15,16,17]]
            max_id_r = np.ndarray.flatten(
                np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.greater, order=20)]))
            max_id = [max_id_r[i] for i in [10,12,32,33]]
        elif samp_id==6:
            min_id = [min_id_r[i] for i in [2,3,4]]
            max_id = max_id_r[2]
        elif samp_id==1:
            min_id = [min_id_r[i] for i in [2,3,6]]
            max_id = [max_id_r[i] for i in [2,4]]
        elif samp_id==2:
              min_id = [min_id_r[i] for i in [2, 3,5]]
              max_id = [max_id_r[i] for i in [1, 2]]
        elif samp_id == 3:
              min_id = [min_id_r[i] for i in [3,4,5]]
              max_id = [max_id_r[i] for i in [3,4]]
        elif samp_id == 5:
            min_id = [min_id_r[i] for i in [2,3,4]]
            max_id = [max_id_r[i] for i in [3]]
        else:
            min_id = min_id_r
            max_id = max_id_r

        peak_id = np.concatenate((min_id, max_id), axis=None)
        dfpeak = df.iloc[peak_id]
        dfpeak = dfpeak.sort_index()

        # --- plot Peaks --- #
        plt.plot(climb_temp, climb_hf, 'g-', label=r'Constant $\Delta$T')
        plt.title('{}wt%'.format(samp))
        plt.xlabel("T (K)")
        plt.ylabel("Q (mW)")
        plt.plot(dfpeak["Sample Temperature(K)"], dfpeak["HeatFlow(mW)"], 'kx',markersize=3)
        # plt.show()
        saveFig = fd + '{}wt%_Q-T_peaks.png'.format(samp)
        plt.savefig(saveFig)
        plt.close()

        saveFile = sd + '{}wt%_peaks.txt'.format(samp)
        with open(saveFile, 'w') as fileOut:
            lines = dfpeak.to_string(index=False, header=True, col_space=0)
            fileOut.writelines(lines)

    ## SPECIFIC HEAT CALCULATIONS

    # Cp = [cp_calc(i, m, T_slope) for i in climb_hf]
    #
    # fig = plt.figure(4, dpi=500)
    # plt.plot(climb_temp, Cp, 'k')
    # points = [50, 700, 1500, 2500, 3597, 3599] #pick out regions with linear Cp
    # plt.plot([climb_temp[i] for i in points], [Cp[i] for i in points], 'ro')
    # plt.title('{}wt%'.format(samp))
    # plt.xlabel('T (K)')
    # plt.ylabel(r'Cp $(\frac{J}{g ⋅ K})$')
    # # plt.show()
    # saveFig = fd + '{}wt%_Cp-T.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close()
    #
    # # --- summary prints --- #
    # print(f'Steps of slope (# steps between yellow and cyan dots): {len(hf_corrected[start_offset:end_offset])}')
    # print(
    #     f'Time of slope (i.e., time in seconds from yellow to cyan dot): {len(hf_corrected[start_offset:end_offset]) * 11.6}')
    # print(f'Number of Heat Capacity data points: {len(Cp)}')
    # print(f'Average Cp ~200K: {np.average(Cp[1600:1800])}')


    ## DETERMINE CRYSTAL FRACTION AND LIQUIDUS TEMPERATURE

    Xl = pd.Series([Xl_calc(T) for z, T in dfclimb['Sample Temperature(K)'].iteritems()])

    Xsamp = samp/100
    Xc = 0
    Fc = (Xl - Xsamp) / (Xl-Xc)
    Fc[Fc<0] = np.nan

    Fc.index = np.arange(start_offset, end_offset,1)
    dfclimb['Fc'] = Fc


    fig = plt.figure()
    plt.plot(dfclimb['Fc'],dfclimb['Sample Temperature(K)'])
    plt.title('{}wt%'.format(samp))
    plt.xlim([0,1])
    plt.xlabel('Fraction of Crystal')
    plt.ylabel('Temperature (K)')
    # plt.show()
    saveFig = fd + '{}wt%_Fc.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    dfliquidus = dfclimb.loc[dfclimb['Fc'].idxmin()]




    ##

    saveFile = sd + '{}wt%_liquidus.txt'.format(samp)
    with open(saveFile, 'w') as fileOut:
        lines = dfliquidus.to_string(index=True, header=True)
        fileOut.writelines(lines)

    saveFile = sd + '{}wt%_climb_proc.txt'.format(samp)
    with open(saveFile, 'w') as fileOut:
        lines = dfclimb.to_string(index=False, header=True)
        fileOut.writelines(lines)

# # debug end line
print('Finished running script')