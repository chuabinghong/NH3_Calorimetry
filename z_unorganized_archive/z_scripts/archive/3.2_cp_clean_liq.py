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
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode

## INPUTS

old_arg = 0
old_wt_arg = 0
# set working, figure and output directories if required
wd = r".."
dd = r"i_data_processed/"
ddb = r"i_data_processed/0.25Kmin/"
fd = r"o_specificHeat/cut/"
# fd = r"t_seafreeze_lh/"
od = r"i_data_processed/"
# od = r"t_seafreeze_lh/"


wtpct_ls = [1.42,2.84,5.11,6.88,8.73,26.96]
# wtpct_ls = [1.42,2.84,6.88,26.96]
mass_ls = [4887.3,4840.2,4629.3,4585.8,4.5500,3777.8]
# mass_ls = [4887.3,4840.2,4585.8,3777.8]
lowC_ls = [250,190,190,110,120,129]
highC_ls = [-1,-1,-1,-1,-1,-1]

file_ls =[dd+f'{i}wt%_cp_liquidPhase.csv' for i in wtpct_ls]


if old_wt_arg ==1:
    dd = r"i_data_processed/old_wt%/"
    od = r"i_data_processed/old_wt%/"
    fd = r"o_specificHeat/old_wt%/cut/"
    wtpct_ls = [1.62, 3.24, 5.83, 7.85, 26.96]
    file_ls =[dd+f'{i}wt%_cp_liquidPhase.csv' for i in wtpct_ls]
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

# --- run following for all samples --- #
for samp_id, samp in enumerate(wtpct_ls):
    # samp = 8.73
    # samp_id = 4
    m = mass_ls[samp_id]
    lowC = lowC_ls[samp_id]
    highC = highC_ls[samp_id]
    df = pd.read_csv(file_ls[samp_id])

## TEST LOW T CUTOFF

    df_c = df[120:-1].copy()
    # plt.close()
    fig = plt.figure()
    plt.plot(df_c['T(K)'], df_c['cp(J/gK)'], 'o', markersize=0.2)
    plt.xlabel('T (K)')
    plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    plt.title('{}wt% low T cutoff test'.format(samp))
    plt.ylim([2, 6])
    # plt.xlim([260, 290])
    plt.legend(markerscale=5)
    # plt.show()
    # saveFig = fd + '{}wt%_cp_zoom.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close()

    ## TEST HIGH T CUTOFF

    df_c = df[1:-1].copy()

    plt.close()
    fig = plt.figure(figsize=(3,2))
    plt.plot(df_c['T(K)'], df_c['cp(J/gK)'], 'o', markersize=0.2)
    plt.xlabel('T (K)')
    plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    plt.title('{}wt% high T cutoff test'.format(samp))
    plt.ylim([2, 6])
    # plt.xlim([230, 250])
    # plt.legend(markerscale=5)
    # plt.show()
    # saveFig = fd + '{}wt%_cp_zoom.png'.format(samp)
    # plt.savefig(saveFig)
    # plt.close()

    ## CHECK BOTH CUTOFF


    df_c = df[lowC:highC].copy()
    plt.close()
    fig = plt.figure(figsize=(3,2))
    plt.plot(df['T(K)'], df['cp(J/gK)'], 'o', markersize=0.2, label='original')
    plt.plot(df_c['T(K)'], df_c['cp(J/gK)'], 'o', markersize=0.2, label='cut')
    plt.xlabel('T (K)')
    plt.ylabel(r'cp $(\frac{J}{g ⋅ K})$')
    plt.title('{}wt% cutoff check'.format(samp))
    plt.ylim([0, 10])
    # plt.xlim([230, 250])
    plt.legend(markerscale=5)
    # plt.show()
    saveFig = fd + '{}wt%_cp_cut_liq.png'.format(samp)
    plt.savefig(saveFig)
    plt.close()

## SAVE

    saveFile = od + '{}wt%_cp_cut_liq.csv'.format(samp)
    df_c.to_csv(saveFile, index=False)

## debug end line
print('Finished running script')