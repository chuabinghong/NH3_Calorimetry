#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 11:42:35 2022

@author: jackdiab
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use(['science','ieee','no-latex'])


df_blank = pd.read_csv("H2O-NH3_blank.csv", skiprows=8)

df_2 = pd.read_csv("H2O-NH3_2wt%.csv", skiprows=8)
df_4 = pd.read_csv("H2O-NH3_4wt%.csv", skiprows=8)
df_6 = pd.read_csv("6wt%_NH4-H2O.csv", skiprows=8)
df_8 = pd.read_csv("H2O-NH3_8wt%.csv", skiprows=8)
df_10 = pd.read_csv("10wt%_NH4-H2O.csv", skiprows=8)
    
sample_data = {2:[df_2, 4.8389], 4:[df_4, 4.7459], 6:[df_6, 4.5983], 8:[df_8, 4.5500], 10:[df_10, 4.3987]}
    
### Heat Capacity
# Cp(T) =
# (dq/dt)/(m × (dT/dt))
# m is mass

def round_dec(x, base=1000000):
    '''Rounds to a specified base.
    base = 1000 --> thousandths place'''
    xr = round((x * base))/base
    return xr

tslope = 0.25/60 # Temperature ramp speed (slope) converted from K/min to K/s
# mW to W, then W/(g * K), 1 W = 1 J/s, ==> J/(s * g * K)
def cp_calc(mW, mass, dTdt):
    '''Calculates Cp from heatflow, mass, and Temp ramp'''
    Cp = np.abs((mW * 0.001)/(mass * dTdt))
    return Cp

def Cp_T(T, a, b, c, d, e):
    '''Shomate Formulation equation for fitting'''
    return a + (b * T) + (c * T**2) + (d * T**3) + (e /(T**2))

def plot(plots='all'):
    if plots=='all':
        samples = sample_data
    elif plots!='all':
        samples = {plots:sample_data[plots]}
    for sample in samples:
        df = sample_data[sample][0] # dataframe of sample
        m = sample_data[sample][1] # mass in grams
        
        fig, ax = plt.subplots(2,2,figsize=(15,10), dpi=500)#, tight_layout= True)
        fig.suptitle(f'{sample} wt.%', fontsize = 25)
        ax[0,0].set_title('Temperature Profile')
        ax[0,0].plot(df["Time(s)"], df["Sample Temperature(K)"], 'c-', label='Sample')
        ax[0,0].plot(df_blank["Time(s)"], df_blank["Sample Temperature(K)"], 'r', label='Blank')
        
        start_offset = 8500 # when ramp starts
        end_offset = 12100 # Ramp end (we set 3600 steps for the ramp)
        ax[0,0].plot(df_blank["Time(s)"][start_offset], df_blank["Sample Temperature(K)"][start_offset], 'yo')
        ax[0,0].plot(df_blank["Time(s)"][end_offset], df_blank["Sample Temperature(K)"][end_offset], 'go')
        ax[0,0].set_xlabel('Time (s)')
        ax[0,0].set_ylabel('Temperature (K)')
        ax[0,0].legend()
        
        #hf_corrected is the heatflow of the sample - heatflow of the blank
        hf_corrected = [df["HeatFlow(mW)"][i] - df_blank["HeatFlow(mW)"][i] for i in range(len(df["HeatFlow(mW)"]))]
        ax[0,1].plot(df["Sample Temperature(K)"], hf_corrected, 'm-', label='Corrected Sample')
        ax[0,1].axhline(0)
        
        # climb_temp is list of temps from the start to end of ramp
        climb_temp = list(df["Sample Temperature(K)"][start_offset:end_offset])
        climb_hf = hf_corrected[start_offset:end_offset]
        
        ax[0,1].plot(climb_temp, climb_hf, 'g--', label=r'Constant $\Delta$T ')
        ax[0,1].set_title(f"Heat Flow vs. {sample} wt.% Sample Temperature")
        ax[0,1].set_xlabel("Temperature (K)")
        ax[0,1].set_ylabel("Heat Flow (mW)")
        ax[0,1].legend()
        
        Cp = [cp_calc(i, m, tslope) for i in climb_hf] # Cp along constant temp ramp
        
        ax[1,0].plot(climb_temp, Cp, 'k')
        Cp_cuts = [100, 700, 1500, 2000, 3597, 3599] # Points along Cp curve that exclude peaks. Helps specify data points for later fit
        ax[1,0].plot([climb_temp[i] for i in Cp_cuts],[Cp[i] for i in Cp_cuts],'ro')
        ax[1,0].set_title("Heat Capacity Over Constant Temperature Ramp")
        ax[1,0].set_xlabel('Sample Temperature (K)')
        ax[1,0].set_ylabel(r'Heat Capacity $(\frac{J}{g * K})$')
        plt.title(f'Cp {sample} wt.%')
        
        # print(f'Steps of slope (# steps between yellow and cyan dots): {len(hf_corrected[start_offset:end_offset])}')
        # print(f'Time of slope (i.e., time in seconds from yellow to cyan dot): {len(hf_corrected[start_offset:end_offset]) * 11.6}')
        
        # Cp(T) = a + bT + cT2 + dT3 + e/T2
        
        climb_temp_np = climb_temp[Cp_cuts[0]:Cp_cuts[1]] + climb_temp[Cp_cuts[2]:Cp_cuts[3]] + climb_temp[Cp_cuts[4]:Cp_cuts[5]]
        Cp_np = Cp[Cp_cuts[0]:Cp_cuts[1]] + Cp[Cp_cuts[2]:Cp_cuts[3]] + Cp[Cp_cuts[4]:Cp_cuts[5]]
        
        # pCp, C_Cp = curve_fit(Cp_T, climb_temp_np, Cp_np, p0=[0.6, 0.0024, -0.000007, 0.00000008, -2000])
        pCp, C_Cp = curve_fit(Cp_T, climb_temp_np, Cp_np, p0=[0,0,0,0,0])
        Cp_pointfit= [Cp_T(i, *pCp) for i in climb_temp_np]
        
        temp_array = np.linspace(100, 300, 500)
        Cp_curve = Cp_T(temp_array, *pCp)
        
        ax[1,1].plot(climb_temp_np, Cp_np, 'k.', label = f"Cp Data ({sample} wt.%)")
        ax[1,1].plot(temp_array, Cp_curve, 'c', label = f"Fit {[round_dec(i,1000) for i in pCp]}")
        ax[1,1].set_title(r'${Cp\: Curve\: \: (Cp(T) = a + bT + cT^2 + dT^3 + \frac{e}{T^2})}$')
        ax[1,1].set_xlabel('Sample Temperature (K)')
        ax[1,1].set_ylabel(r'Heat Capacity $(\frac{J}{g * K})$')
        ax[1,1].legend()
        print(f"{sample} wt.% Fit Params [a, b, c, d, e]: {[round_dec(i) for i in pCp]}")

plot(10)