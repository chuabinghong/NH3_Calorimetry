"""
created by Bing Hong CHUA 27Jul22

script objective:
calculate specific heat
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline


plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FUNCTIONS

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

## INPUTS

old_arg = 0
old_wt_arg = 0
shift_arg = 0
ice_arg = 1

# set working, figure and output directories if required
wd = r".."
dd = r"i_data_processed/"
ddb = r"i_data_processed/0.25Kmin/"
fd = r"o_specificHeat/"
# fd = r"t_seafreeze_lh/"
od = r"i_data_processed/"
# od = r"t_seafreeze_lh/"

wtpct_ls = [1.42,2.84,5.11,6.88,8.73,26.96]
mass_ls_in = [4887.3,4840.2,4629.3,4585.8,4520.2,3777.8]
# wtpct_ls = [1.76]
# mass_ls = [4887.3]

file_ls = [dd+f'{i}wt%_climb_proc.txt' for i in wtpct_ls]
mass_ls = [x*0.001 for x in mass_ls_in]
peak_ls = [11851,20]

if old_wt_arg ==1:
    dd = r"i_data_processed/old_wt%/"
    od = r"i_data_processed/old_wt%/"
    fd = r"o_specificHeat/old_wt%/"
    wtpct_ls = [1.62, 3.24, 5.83, 7.85, 26.96]
    file_ls = [dd+f'{i}wt%_climb_proc.txt' for i in wtpct_ls]
    mass_ls = [x * 0.001 for x in mass_ls_in]

# set working and save directories for 0.25Kmin
if old_arg==1:
    dd = r"i_data_processed/0.25Kmin/"
    od = r"i_data_processed/0.25Kmin/"
    fd = r"o_specificHeat/0.25Kmin/"
    # 0.25Kmin
    wtpct_ls = [1.6,3.2,6,8,10]
    file_ls = ["1.6wt%_climb_proc.txt", "3.2wt%_climb_proc.txt", "6wt%_climb_proc.txt", "8wt%_climb_proc.txt", "10wt%_climb_proc.txt"]
    file_ls = [dd + x for x in file_ls]
    mass_ls = [4.8389, 4.7648, 4.5983, 4.5500, 4.3987]
    peak_ls = [11858,11780,11682,11484,11476]

os.chdir(wd)

## IMPORT DATA FILES

latent_heat = 334 #water-ice, joules/gram

# if arg = 1, get ice cp from our own 0wt% experiments
if ice_arg == 1:
    col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
    df_0 = pd.read_csv(dd + '0wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)
    m = 2.511

    Cp = np.zeros([len(df_0), 1]); Cp[:] = np.nan
    Cp_T = np.zeros([len(df_0), 1])

    for i in range(len(df_0) - 1):
        Q_integral = -(df_0['Q_corrected(mW)'].iloc[i + 1] / 1000 + df_0['Q_corrected(mW)'].iloc[i] / 1000) * (
                df_0['time(s)'].iloc[i + 1] - df_0['time(s)'].iloc[i]) / 2
        dT = df_0['sampleT(K)'].iloc[i + 1] - df_0['sampleT(K)'].iloc[i]

        Cp[i] = Q_integral / (m * dT)
        Cp_T[i] = (df_0['sampleT(K)'].iloc[i + 1] + df_0['sampleT(K)'].iloc[i]) / 2

    df_cp = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'cp(J/gK)': np.ravel(Cp)})
    df_cp = df_cp[700:-750].copy()

    df_iceCp = pd.read_csv(dd + 'ice_Cp.csv', names=['T', 'Cp'])
    iceCp_f = UnivariateSpline(df_cp['T(K)'], df_cp['cp(J/gK)'], k=3,ext=0)

else:
    # ice specific heat based on Feistel and Wagner 2006
    df_iceCp = pd.read_csv(dd+'ice_Cp.csv', names=['T','Cp'])
    iceCp_f = interp1d(df_iceCp['T'],df_iceCp['Cp']/1000)

# ice sublimation latent heat based on Feistel and Wagner 2007
# df_iceLh = pd.read_csv(dd+'ice_Lh.csv', names=['T','Lh'])
# iceLh_f = interp1d(df_iceLh['T'],df_iceLh['Lh']/1000)

# ice latent heat of melting based on SeaFreeze derivation (>=240K) already in J/g
df_iceLh = pd.read_csv(dd+'SF_Lh.csv',skiprows=1, names=['T','Lh'])
iceLh_f = interp1d(df_iceLh['T'],df_iceLh['Lh'],fill_value='extrapolate')

# # plot Lh
# fig = plt.figure()
# plt.plot(np.arange(180,300),iceLh_f(np.arange(180,300)),label='Extrapolation')
# plt.plot(df_iceLh['T'],df_iceLh['Lh'],linewidth=1.5,label='SeaFreeze')
# plt.plot(273.15,333.55,'x',color='k')
# plt.xlabel('T (K)')
# plt.ylabel(r'$\Delta$H (J $g^{-1}$)')
# plt.title('Latent Heat of Melting (Pure Ice)')
# plt.legend()
# # plt.show()
# saveFig = fd + 'Lh_SF.png'
# plt.savefig(saveFig)
# plt.close()

# --- run following for all samples --- #
for samp_id, samp in enumerate(wtpct_ls):
    # samp_id = 4
    # samp = 8.73
    m = mass_ls[samp_id]
    col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Fc']
    df = pd.read_csv(file_ls[samp_id], skiprows=1, names=col_names, delim_whitespace=True)

    # samp = 2.842
    # m = 0.9642
    # df = pd.read_csv(r'i_data_processed/2.842wt%_climb_proc.txt', skiprows=1, names=col_names, delim_whitespace=True)

## 1) CALCULATE CP VIA LINEAR COMBINATION OF HEAT USED IN MELTING ICE, HEATING ICE AND HEATING LIQUID
# in this section, heatflow peak has not been aligned with liquidus. we expect an overestimate of the contributions of
# ice melt.

    if shift_arg ==1:
        peak_id = peak_ls[samp_id]
        FcLast_id = df['index'][df["Fc"].last_valid_index()]
        df['Q_corrected(mW)2'] = df['Q_corrected(mW)'].shift(periods=FcLast_id-peak_id)

    dfB = df[df["Fc"].first_valid_index():df["Fc"].last_valid_index()].copy()

    # arrays to plot heat budget
    c_melt = np.zeros([len(dfB),1]); c_melt[:] = np.nan
    c_heat = np.zeros([len(dfB),1]); c_heat[:] = np.nan
    l_heat = np.zeros([len(dfB),1]); l_heat[:] = np.nan
    j_total= np.zeros([len(dfB),1]); j_total[:] = np.nan
    j_sum= np.zeros([len(dfB),1]); j_sum[:] = np.nan

    Cp = np.zeros([len(dfB),1]); Cp[:] = np.nan
    Cp_T = np.zeros([len(dfB),1])
    for i in range(len(dfB)-1):
        Q_integral = -(dfB['Q_corrected(mW)'].iloc[i+1]/1000 + dfB['Q_corrected(mW)'].iloc[i]/1000) * (
                    dfB['time(s)'].iloc[i+1] - dfB['time(s)'].iloc[i]) / 2
        dfFc = dfB['Fc'].iloc[i] - dfB['Fc'].iloc[i+1]
        dT = dfB['sampleT(K)'].iloc[i+1] - dfB['sampleT(K)'].iloc[i]

        # Cp[i] = (Q_integral - dfFc * m * latent_heat - dfB['Fc'].iloc[i] * m * iceCp_f(dfB['sampleT(K)'].iloc[i+1]) * dT) / (m * (1 - dfB['Fc'].iloc[i+1]) * dT)
        Cp[i] = (Q_integral - dfFc * m * iceLh_f(dfB['sampleT(K)'].iloc[i+1]) - dfB['Fc'].iloc[i] * m * iceCp_f(dfB['sampleT(K)'].iloc[i+1]) * dT) / (m * (1 - dfB['Fc'].iloc[i+1]) * dT)
        Cp_T[i] = (dfB['sampleT(K)'].iloc[i+1] + dfB['sampleT(K)'].iloc[i]) / 2

        j_total[i] = Q_integral
        # c_melt[i] = dfFc * m * latent_heat
        c_melt[i] = dfFc * m * iceLh_f(dfB['sampleT(K)'].iloc[i+1])
        c_heat[i] = dfB['Fc'].iloc[i] * m * iceCp_f(dfB['sampleT(K)'].iloc[i+1]) * dT
        l_heat[i] = m * (1 - dfB['Fc'].iloc[i+1]) * dT * Cp[i]
        j_sum[i] = c_melt[i] + c_heat[i] + l_heat[i]

    ## SAVE CP VALUES
    Xl_calc_v = np.vectorize(Xl_calc)
    Cp_X = Xl_calc_v(Cp_T)
    Cp_m = Cp_X/17.031/((1-Cp_X)/1000)

    dfCp = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'X': np.ravel(Cp_X), 'm': np.ravel(Cp_m),'cp(J/gK)': np.ravel(Cp)})
    dfCp['cpm(J/molK)'] = dfCp['cp(J/gK)'] * (17.031 * dfCp['X'] + 18.0528 * (1 - dfCp['X']))
    dfCp = dfCp.dropna()
    saveFile = od + '{}wt%_cp.csv'.format(samp)
    dfCp.to_csv(saveFile, index=False)

## EXTRACT CP IN LIQUID ONLY PHASE
    dfL = df[df["Fc"].last_valid_index()+1:df.index[-1]].copy()

    Cp2 = np.zeros([len(dfL),1]); Cp2[:] = np.nan
    Cp_T2 = np.zeros([len(dfL),1])
    for i in range(len(dfL)-1):
        Q_integral = -(dfL['Q_corrected(mW)'].iloc[i+1]/1000 + dfL['Q_corrected(mW)'].iloc[i]/1000) * (
                    dfL['time(s)'].iloc[i+1] - dfL['time(s)'].iloc[i]) / 2
        dT = dfL['sampleT(K)'].iloc[i+1] - dfL['sampleT(K)'].iloc[i]

        Cp2[i] = Q_integral  / (m * dT)
        Cp_T2[i] = (dfL['sampleT(K)'].iloc[i+1] + dfL['sampleT(K)'].iloc[i]) / 2

    ## SAVE CP VALUES
    Cp_m2 = samp/100 / 17.031 / ((1 - samp/100) / 1000)
    dfCp2 = pd.DataFrame({'T(K)': np.ravel(Cp_T2), 'X': samp/100, 'm': Cp_m2, 'cp(J/gK)': np.ravel(Cp2)})
    dfCp2['cpm(J/molK)'] = dfCp2['cp(J/gK)'] * (17.031 * dfCp2['X'] + 18.0528 * (1 - dfCp2['X']))
    dfCp2 = dfCp2.dropna()
    saveFile = od + '{}wt%_cp_liquidPhase.csv'.format(samp)
    dfCp2.to_csv(saveFile, index=False)

    ## PLOT FIGURES

    # CP NO ZOOM
    fig = plt.figure()
    plt.plot(Cp_T,Cp,'o',markersize=0.2, label='Partial Melting')
    plt.plot(Cp_T2,Cp2,'o',markersize=0.2, label='Liquid Phase')
    plt.xlabel('T (K)')
    plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    plt.title('{}wt% no offset'.format(samp))
    plt.legend(markerscale=5)
    # plt.show()
    saveFig = fd + '{}wt%_cp.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # CP ZOOM
    plt.figure()
    plt.plot(Cp_T, Cp, 'o', markersize=0.2, label='Partial Melting')
    plt.plot(Cp_T2, Cp2, 'o', markersize=0.2, label='Liquid Phase')
    plt.xlabel('T (K)')
    plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    plt.title('{}wt% no offset'.format(samp))
    plt.ylim([2,7])
    plt.legend(markerscale=5)
    # plt.show()
    saveFig = fd + '{}wt%_cp_zoom.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # CPM NO ZOOM
    plt.figure()
    plt.plot(dfCp['T(K)'],dfCp['cpm(J/molK)'],'o',markersize=0.2)
    plt.xlabel('T (K)')
    plt.ylabel(r'cpm $(\frac{J}{mol ⋅ K})$')
    plt.title('{}wt% no offset'.format(samp))
    # plt.show()
    saveFig = fd + '{}wt%_cpm.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # CPM ZOOM
    plt.figure()
    plt.plot(dfCp['T(K)'],dfCp['cpm(J/molK)'],'o',markersize=0.2)
    plt.xlabel('T (K)')
    plt.ylabel(r'cpm $(\frac{J}{mol ⋅ K})$')
    plt.title('{}wt% no offset'.format(samp))
    plt.ylim([0, 120])
    # plt.show()
    saveFig = fd + '{}wt%_cpm_zoom.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # HEAT BUDGET
    plt.figure()
    plt.axhline(y=0,color='k')
    plt.plot(Cp_T,j_total,'k',label='total')
    plt.plot(Cp_T,c_melt,'r',label='phase change')
    plt.plot(Cp_T,c_heat,'m',label='heat solid')
    plt.plot(Cp_T,l_heat,'b',label='heat liquid')
    plt.plot(Cp_T,j_sum,'y--',label='sum')
    plt.ylim([-1,2])
    plt.xlabel('T (K)')
    plt.ylabel(r'Heat (J)')
    plt.title('{}wt% no offset'.format(samp))
    plt.legend(fontsize=5)
    # plt.show()
    saveFig = fd + '{}wt%_cp_0_heatBudget.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

    # STACKED HEAT BUDGET
    plt.figure()
    plt.title('{}wt% normalized heat budget'.format(samp))
    plt.stackplot(Cp_T.flatten(),c_melt.flatten()/j_total.flatten(),c_heat.flatten()/j_total.flatten(),  l_heat.flatten()/j_total.flatten(),labels=['melt ice', 'heat ice', 'heat liquid'])
    # plt.stackplot(Cp_T.flatten(), c_heat.flatten()/j_total.flatten(), labels=['melt ice', 'heat ice', 'heat liquid'])
    plt.xlabel('T (K)')
    plt.ylim([0,1])
    plt.legend(fontsize=5,frameon = 1,loc='lower right')
    # plt.show()
    saveFig = fd + '{}wt%_cp_0_heatBudget_normStacked.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

## PLOT ALL WRT TIME

    # plt.figure()
    # plt.plot(dfB['time(s)'],dfB['Q(mW)'],label='Raw Heat Flow')
    # plt.plot(dfB['time(s)'],dfB['Q_corrected(mW)'],label='Corrected Heat Flow')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flow (mW)')
    # plt.legend()
    # # plt.plot(dfB['time(s)'],dfB['sampleT(K)'])
    # # plt.xlabel('Time (s)')
    # # plt.ylabel('T (K_')
    # # plt.ylim([-25,-12])
    # saveFig = fd + 't_heatflow_periodicity.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close()



    ## 1) CALCULATE CP VIA LINEAR COMBINATION OF HEAT USED IN MELTING ICE, HEATING ICE AND HEATING LIQUID
    # in this section, heatflow peak is arbitrarily aligned with liquidus. This may not be the true temperature offset

#     peak_id = peak_ls[samp_id]
#     FcLast_id = df['index'][df["Fc"].last_valid_index()]
#     df['Q_shift'] = df['Q_corrected(mW)'].shift(periods=FcLast_id-peak_id)
#
#     dfB = df[df["Fc"].first_valid_index():df["Fc"].last_valid_index()].copy()


## debug end line
print('Finished running script')