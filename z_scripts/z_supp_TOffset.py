import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.optimize import curve_fit
import base

plt.style.use(['science', 'nature', 'no-latex'])
matplotlib.rcParams['axes.unicode_minus'] = False
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.use('Qt5Agg')  # fixes matplotlib unresponsive issues while debugging in cell mode
matplotlib.rcParams['mathtext.default'] = 'regular'

## FILE DIRECTORY SETTINGS

wd = r"../"
os.chdir(wd)

##
def tanh_eqn(T, a, b, c,d):
    """tanh equation for fitting"""
    return a*np.tanh((T+b)*c)+d

data =pd.read_csv(r"i_data\T_offset\phase_transition_T.csv", header=0)

out, pcov = curve_fit(tanh_eqn,data['Theoretical'],data['Offset'],p0=[.3,-250,.015,.3])

T = np.arange(100,500)
plt.scatter(data['Theoretical'],data['Offset'],label='measured T of transition')
# plt.plot(T,.3*np.tanh((T-250)*.015)+.3)
plt.plot(T,tanh_eqn(T,*out),'k',label='tanh fit')
eqn = '$y = %.3ftanh[%.3f(T%.3f)]+%.3f$' % (out[0],out[2],out[1],out[3])
plt.annotate(eqn,(300,0),ha='left',va='center',fontsize=4)


plt.title('tanh fit to phase transition offsets')
plt.xlabel('Theoretical T of Transition (K)')
plt.ylabel('Measured T - Theoretical T (K)')
plt.savefig(r'o_supplementaryPlots\z_T_offset.png')
# base.show_plot_max()

## debug end line
print('Finished running script')