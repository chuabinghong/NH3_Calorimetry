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

dd = r"i_data/0.10Kmin/"
sd = r"i_data_processed/"
fd = r"o_heatFlow/"

## SCRIPT ARGUMENTS
peak_arg = 1 #1 if we are getting peaks
range_arg = 0 #1 if using -196C to 46C
wt_arg = 1 # 0 for lab measurements, 1 for liquidus alignment, 2 for liquidus alignment -0.5K
smallHF_arg = 0 #1 if using small heat flow experiments
ramp_rate = 0.1


mass_ls = [2.5108,4.8873,4.8402,4.6293,4.5858,4.5202,3.7778]
if wt_arg == 0:
    wt_ls = [0, 1.42, 2.84, 5.11, 6.88, 8.73, 25.26]  # based off lab measurements
elif wt_arg ==1:
    wt_ls = [0, 1.42, 2.5, 4.9, 8.1, 9.7, 26.9]  # based off liquidus alignment
elif wt_arg ==2:
    wt_ls = [0, 1.42, 2.9, 5.2, 8.4, 10.0, 27.0]  # based off liquidus alignment
t_end = 12981

# adjust script arguments according to range_arg
if range_arg == 1:
    dd = r"i_data/46C/"
    sd = r"i_data_processed/46C/"
    fd = r"o_heatFlow/46C/"
    mass_ls = [4.8505,4.6844,4.5386,4.1943,3.8153,3.7107]
    if wt_arg == 0:
        wt_ls = [2.00,4.03,9.45,16.00,21.01] # based off lab measurements
    elif wt_arg == 1:
        wt_ls = [2.0, 2.9, 4.91, 7.9, 14.2, 20]  # based off liquidus alignment
    elif wt_arg == 2:
        wt_ls = [2.0, 3.3, 5.21, 8.2, 14.4, 20.2]  # based off liquidus alignment
    t_end = 13066

## DATA PROCESSING

# import blank
blank = base.DataRaw(0,'blank',ramp_rate,t_end)
blank.df = pd.read_csv(dd + 'Blank.csv', skiprows=8)
blank.correct_HF(blank.df)
blank.mean_sampling(dd)

blankdf = pd.read_csv(dd + 'mean_H2O-NH3_0g.csv', header=0)
for idx, wt in enumerate(wt_ls):
    # idx = 3 #TODO for single data debugging
    # wt = wt_ls[idx] #TODO for single data debugging

    m = mass_ls[idx]

    # data processing
    data = base.DataRaw(m, wt,ramp_rate,t_end)
    data.import_raw(dd)
    data.mean_sampling(dd)

