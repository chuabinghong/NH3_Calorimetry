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

# set working and save directories if required
wd = r"../i_data/"
sd = r"../"
os.chdir(wd)

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

df_blank = pd.read_csv("H2O-NH3_blank.csv", skiprows=8)

fig = plt.figure(1)
# ax = fig.subplot
plt.plot(df_blank["Furnace Temperature(K)"], df_blank["HeatFlow(mW)"], 'r--', label='Furnance')
plt.plot(df_blank["Sample Temperature(K)"], df_blank["HeatFlow(mW)"], 'c-', label='Sample')
plt.xlabel('T (K)')
plt.ylabel('Q (mW)')
# plt.title('testing')
plt.legend()
# plt.show()
saveFig = sd + 'blank_Q-T.png'
plt.savefig(saveFig)
plt.close()

## IMPORT DATA FILES

wtpct_ls = [0,1.76,3.52,6,8,10]
file_ls = ["H2O-NH3_0wt%.csv","H2O-NH3_2wt%.csv", "H2O-NH3_4wt%.csv", "H2O-NH3_6wt%.csv", "H2O-NH3_8wt%.csv", "H2O-NH3_10wt%.csv"]
mass_ls = [4.9859,4.8389, 4.7459, 4.5983, 4.5500, 4.3987]

# wtpct_ls = [1.76,3.52,6,8,10]
# file_ls = ["H2O-NH3_2wt%.csv", "H2O-NH3_4wt%.csv", "H2O-NH3_6wt%.csv", "H2O-NH3_8wt%.csv", "H2O-NH3_10wt%.csv"]
# mass_ls = [4.8389, 4.7459, 4.5983, 4.5500, 4.3987]

wtpct_ls = [3.52]
file_ls = ["H2O-NH3_4wt%_2.csv"]
mass_ls = [4.7648]

# --- run following for all samples --- #
for samp_id, samp in enumerate(wtpct_ls):
    m = mass_ls[samp_id]
    df = pd.read_csv(file_ls[samp_id], skiprows=8)


    ## T OVER t

    fig = plt.figure(2)
    plt.plot(df["Time(s)"], df["Sample Temperature(K)"], 'c-', label='Sample')
    plt.plot(df_blank["Time(s)"], df_blank["Sample Temperature(K)"], 'r', label='Blank')
    plt.xlabel('t (s)')
    plt.ylabel('T (K)')
    plt.title('{}wt%'.format(samp))
    plt.legend()
    # plt.show()
    saveFig = sd + '{}wt%_T-t.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()


    # --- hardcode start and end t --- #
    start_offset = 8500
    end_offset = 12100

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
    plt.plot(climb_temp, climb_hf, 'g', label=r'Constant $\Delta$T ')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat Flow (mW)")
    plt.title(r'{}wt% $NH_3$'.format(samp))
    # plt.legend()
    # plt.show()
    saveFig = sd + '{}wt%_Q-T_2.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    ## FIND PEAKS (MANUAL FOR EACH RUN)

    # min_id_r = np.ndarray.flatten(np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.less,order=90)])) #find minima peak indexes, plus index offset
    # max_id_r = np.ndarray.flatten(np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.greater,order=90)])) #find minima peak indexes, plus index offset
    #
    #
    # if samp ==1.76:
    #     min_id = min_id_r[1:4]
    #     max_id = max_id_r[0:2]
    # elif samp ==0:
    #     min_id = min_id_r[0:4]
    #     max_id = max_id_r[0:2]
    # elif samp ==3.52:
    #     min_id_r = np.ndarray.flatten(
    #         np.asarray([start_offset + x for x in argrelextrema(np.array(climb_hf), np.less, order=50)]))
    #     min_id = min_id_r[1:5]
    #     max_id = max_id_r[0:2]
    # elif samp ==6:
    #     min_id = min_id_r[2:6]
    #     max_id = max_id_r[0:2]
    # elif samp ==8:
    #     min_id = min_id_r[2:6]
    #     max_id = max_id_r[0:2]
    # elif samp ==10:
    #     min_id = min_id_r[2:5]
    #     max_id = max_id_r[1:3]
    #
    # peak_id = np.concatenate((min_id, max_id), axis=None)
    # dfpeak = df.iloc[peak_id]
    # dfpeak = dfpeak.sort_index()
    #
    # # --- plot Q-T --- #
    # plt.plot(climb_temp, climb_hf, 'g-', label=r'Constant $\Delta$T')
    # plt.title('{}wt%'.format(samp))
    # plt.xlabel("T (K)")
    # plt.ylabel("Q (mW)")
    # plt.plot(dfpeak["Sample Temperature(K)"], dfpeak["HeatFlow(mW)"], 'kx',markersize=3)
    # # plt.show()
    # saveFig = sd + '{}wt%_Q-T.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close()

    # saveFile = sd + 'data_proc/{}wt%_peaks.txt'.format(samp)
    # with open(saveFile, 'w') as fileOut:
    #     lines = dfpeak.to_string(index=False, header=True, col_space=0)
    #     fileOut.writelines(lines)

## SPLINE INTERPOLATION

    # from scipy.interpolate import interp1d
    # f = interp1d(dfclimb["Sample Temperature(K)"], dfclimb['HeatFlow(mW)'])
    #
    # from scipy.interpolate import UnivariateSpline
    # y_spl = UnivariateSpline(dfclimb["Sample Temperature(K)"], dfclimb['HeatFlow(mW)'], s=0, k=4)
    #
    # dev = y_spl.derivative(n=2)
    #
    #
    # fig = plt.figure(9)
    # plt.plot(dfclimb["Sample Temperature(K)"][1200:1400],y_spl(dfclimb["Sample Temperature(K)"][1200:1400]))
    # # plt.plot(dfclimb["Sample Temperature(K)"],dfclimb["HeatFlow(mW)"])
    #
    # fig = plt.figure(8)
    # plt.plot(dfclimb["Sample Temperature(K)"][1200:1400],dev(dfclimb["Sample Temperature(K)"][1200:1400]))



    ## SPECIFIC HEAT CALCULATIONS

    Cp = [cp_calc(i, m, T_slope) for i in climb_hf]

    fig = plt.figure(4, dpi=500)
    plt.plot(climb_temp, Cp, 'k')
    points = [50, 700, 1500, 2500, 3597, 3599] #pick out regions with linear Cp
    plt.plot([climb_temp[i] for i in points], [Cp[i] for i in points], 'ro')
    plt.title('{}wt%'.format(samp))
    plt.xlabel('T (K)')
    plt.ylabel(r'Cp $(\frac{J}{g ⋅ K})$')
    # plt.show()
    saveFig = sd + '{}wt%_Cp-T.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # --- summary prints --- #
    print(f'Steps of slope (# steps between yellow and cyan dots): {len(hf_corrected[start_offset:end_offset])}')
    print(
        f'Time of slope (i.e., time in seconds from yellow to cyan dot): {len(hf_corrected[start_offset:end_offset]) * 11.6}')
    print(f'Number of Heat Capacity data points: {len(Cp)}')
    print(f'Average Cp ~200K: {np.average(Cp[1600:1800])}')

## FIT SPECIFIC HEAT CURVE IN REGIONS WITH LINEAR CP OVER T

    # Cp(T) = a + bT + cT2 + dT3 + e/T2

    Cp_cuts = points
    climb_temp_np = climb_temp[Cp_cuts[0]:Cp_cuts[1]] + climb_temp[Cp_cuts[2]:Cp_cuts[3]]\
                    + climb_temp[Cp_cuts[4]:Cp_cuts[5]]
    Cp_np = Cp[Cp_cuts[0]:Cp_cuts[1]] + Cp[Cp_cuts[2]:Cp_cuts[3]] + Cp[Cp_cuts[4]:Cp_cuts[5]]

    pCp, C_Cp = curve_fit(Cp_T, climb_temp_np, Cp_np, p0=[0, 0, 0, 0, 0]) #TODO figure how curve_fit works
    Cp_pointfit = [Cp_T(i, *pCp) for i in climb_temp_np]


    print(*pCp)
    temp_array = np.linspace(100, 300, 500)
    Cp_curve = Cp_T(temp_array, *pCp)

    plt.figure(5, dpi=500)
    plt.plot(climb_temp_np, Cp_np, 'k.', label="Sample Heat Capacity (peaks removed)")
    plt.plot(temp_array, Cp_curve, 'c', label="Fit")
    plt.title(r'${Cp\: Curve\: \: (Cp(T) = a + bT + cT^2 + dT^3 + \frac{e}{T^2})}$')
    plt.xlabel('Sample Temperature (K)')
    plt.ylabel(r'Heat Capacity $(\frac{J}{g ⋅ K})$')
    plt.legend()
    # plt.show()
    saveFig = sd + '{}wt%_Cp_fit.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()


    print(f"Fit (Params (a, b, c, d): {pCp}")


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
    saveFig = sd + '{}wt%_Fc.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    dfliquidus = dfclimb.loc[dfclimb['Fc'].idxmin()]




    ##

    saveFile = sd + 'data_proc/{}wt%_liquidus.txt'.format(samp)
    with open(saveFile, 'w') as fileOut:
        lines = dfliquidus.to_string(index=True, header=True)
        fileOut.writelines(lines)

    saveFile = sd + 'data_proc/{}wt%_climb_proc.txt'.format(samp)
    with open(saveFile, 'w') as fileOut:
        lines = dfclimb.to_string(index=False, header=True)
        fileOut.writelines(lines)

# # debug end line
print('Finished running script')