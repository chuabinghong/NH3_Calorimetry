import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import base
import matplotlib
import pandas as pd

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'figure.dpi': '600'})
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('Qt5Agg')

def fitfunction(v, a, b, c, d):
    return (a * np.tanh((b * v) + c) + d)


# x data = x
# y data = np.log10(y)

# read data from file
readdata = True

if (readdata == True):
    # arr = np.loadtxt("anode_90ma_test1.csv", delimiter=",")
    df = pd.read_csv(r'C:\1_data\OneDrive - Nanyang Technological University\OFYP\CalorimetryAnalysis\i_data\others\toffset.csv',header=0)
    x = df['T']
    # option of using raw data or log of raw data
    y = df['diff']
    # y=np.log10(arr[:,1])

# else:
#     a_coeff = 100
#     b_coeff = 1
#     c_coeff = 0.2
#     d_coeff = -99
#     x = np.arange(-4, 4, 0.1)
#     createvalues = np.vectorize(fitfunction)
#     y = createvalues(x, a_coeff, b_coeff, c_coeff, d_coeff)
#     y = y + np.random.normal(-100,100,len(y))
#     # optionally take log10
#     # y=np.log10(y)

# pars, cov = curve_fit(fitfunction,x,np.log10(y))
pars, cov = curve_fit(fitfunction, x, y)
print("fit pars = ", pars)

plt.plot(x, fitfunction(x, *pars), 'r-', linewidth='1', label='Line of Best Fit')
plt.scatter(x, y,20, marker='.', label='Data')

plt.title('Graph of Line of Best Fit of Anode Current against Anode Potential')
plt.grid(True)
plt.xlabel('Voltage (V)')
plt.ylabel('Current (A)')
plt.legend()
base.show_plot_max()
plt.show()