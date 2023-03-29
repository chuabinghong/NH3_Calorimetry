import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import base
from scipy.optimize import curve_fit


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
fd = r"o_specificHeat/"

# data =pd.read_csv(r"C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data_processed\20.1wt%_hf_prepped_3.7107g.csv", header=0)
##
# plt.scatter(data['sampleT(K)'],data['Q_corrected(mW)'], 1, marker='o')
# base.show_plot_max()

def tanh_eqn(T, a, b, c,d):
    """Shomate Formulation equation for fitting"""
    return a*np.tanh((T+b)*c)+d

data =pd.read_csv(r"C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data\others\liquidus.csv", header=0)

out, pcov = curve_fit(tanh_eqn,data['Theoretical'],data['Offset'],p0=[.3,-250,.015,.3])

T = np.arange(100,500)
plt.scatter(data['Theoretical'],data['Offset'])
# plt.plot(T,.3*np.tanh((T-250)*.015)+.3)
plt.plot(T,tanh_eqn(T,*out),'k')
plt.title('tanh fit to phase transition offsets')
plt.xlabel('Theoretical T of Transition (K)')
plt.ylabel('Measured T - Theoretical T (K)')
plt.savefig(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data\others\liquidus.png')
base.show_plot_max()

## debug end line
print('Finished running script')