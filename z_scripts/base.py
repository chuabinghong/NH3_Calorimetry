import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

## EXPERIMENTAL CONSTANTS
import base

Q_err = 1e-7
m_err = 1e-4
T_err = 1e-4
## FUNCTIONS

def mol2mass(molfrac):
    H2O = 18.01528
    NH3 = 17.031
    mol_in_g = molfrac * NH3 + (1 - molfrac) * H2O
    wt = molfrac * NH3 / mol_in_g
    return(wt)

def splt_axis_label(fig,xlabel,ylabel):
    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def get_closest_df_value(df,column,value,length):
    idx = (df[column] - value).abs().argsort()[:length]
    out_df= df.iloc[idx]
    return(out_df)

def show_plot_max():
    plt.show(block=False)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

def Xl_calc(T):
    """Croft equation for liquidus of water-rich ammonia-water mixutre, up till peritectic (X<=0.329)
    input: current tempertaure
    output: ammonia weight%"""
    lqd_coeff = [56695, -46269, 11842, -1651.4, -53.07, 273.15 - T]  # CROFT 1987
    roots = np.roots(lqd_coeff) # get roots of the polynomial
    roots_real = np.real(roots[np.iscomplex(roots) == False]) #remove imaginary roots
    root_actual = roots_real[(0 < roots_real) & (roots_real <= 0.329)] #keep roots within peritectic
    if len(root_actual) == 0:
        root_actual = np.nan
    else:
        root_actual = root_actual[0]
    return root_actual

def T_calc(Xl):
    """Croft equation for liquidus of water-rich ammonia-water mixutre, up till peritectic (X<=0.329)
    input: ammonia weight%
    output: Temperature"""
    T = 273.15 - 53.07*Xl - 1651.4*Xl**2 + 11842*Xl**3 - 46269*Xl**4 + 56695*Xl**5
    return T

def Fc_calc(T,wt):
    """calculate crystal fraction using Croft 1988 equation"""
    Xl = pd.Series([Xl_calc(T) for z, T in T.iteritems()])

    Xsamp = wt / 100
    Xc = 0
    Fc = (Xl - Xsamp) / (Xl - Xc)
    Fc[Fc < 0] = np.nan
    return Fc

def shomate_eqn(T, a, b, c, d, e):
    """Shomate Formulation equation for fitting"""
    return a + (b * T) + (c * T**2) + (d * T**3) + (e /(T**2))

def quad_eqn(T, a, b, c):
    """Shomate Formulation equation for fitting"""
    return a + (b * T) + (c * T**2)


def plot_ci_manual(t, s_err, n, x, x2, y2, colour, ax=None):
    """Return an axes of confidence bands using a simple approach.
    """
    if ax is None:
        ax = plt.gca()

    ci = t * s_err * np.sqrt(1 / n + (x2 - np.mean(x)) ** 2 / np.sum((x - np.mean(x)) ** 2))
    ax.fill_between(x2, y2 + ci, y2 - ci, color=colour, linewidth=0,alpha=0.2)

    return ax

def plot_PI_melt(df,eqn,p,ax,colour):
    import scipy.stats as stats

    y_model = eqn(df['T(K)'], *p)
    x = df['T(K)']
    y = df['cp(J/gK)']

    n = df['T(K)'].size
    m = p.size
    dof = n - m
    t = stats.t.ppf(0.975, n - m)

    # Estimates of Error in Data/Model
    resid = y - y_model  # residuals; diff. actual data from predicted values
    chi2 = np.sum((resid / y_model) ** 2)  # chi-squared; estimates error in data
    chi2_red = chi2 / dof  # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid ** 2) / dof)  # standard deviation of the error

    x2 = np.linspace(184.5,233.5,500)
    y2 = eqn(x2,*p)

    pi = t * s_err * np.sqrt(1 + 1 / n + (x2 - np.mean(x)) ** 2 / np.sum((x - np.mean(x)) ** 2))
    ax.fill_between(x2, eqn(x2, *p) + pi, eqn(x2, *p) - pi, color=colour, alpha=.1, linewidth=0, zorder=-10)
    ax.plot(x2, eqn(x2, *p) - pi, "--", linewidth=.4, color=colour)
    ax.plot(x2, eqn(x2, *p) + pi, "--", linewidth=.4,color=colour)





## #TODO class
# noinspection PyUnresolvedReferences
class DataRaw:
    def __init__(self,mass,weight_pct,ramp_rate,t_end):
        self.m = mass
        self.wt = weight_pct
        self.dT = ramp_rate
        self.t_start = []
        self.t_end = t_end
        # elif self.dT == 0.25: #TODO add value for 0.25K/min
        self.df = []
        self.df2 = []
        self.liquidus = []
        self.df_mean = []

    def import_raw(self,dd):
        """import raw output of calorimeter, then crop data to within ramp sequence"""
        self.df =  pd.read_csv(dd + 'H2O-NH3_{}g.csv'.format(self.m), skiprows=8)
        self.t_start = self.df['Sample Temperature(K)'].idxmin()
        self.df2 = self.df[self.t_start:self.t_end].copy()

    def import_mean(self,dd,df_blank):
        """import output of calorimeter after mean averaging"""
        self.df2 =  pd.read_csv(dd + 'mean_H2O-NH3_{}g.csv'.format(self.m), header=0)
        if self.wt == 'blank':
            self.df['Q_Corrected'] = self.df["HeatFlow(mW)"]
        elif self.wt != 'blank':
            self.df['Q_Corrected'] = self.df["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"]
        self.df2 = self.df

    def correct_HF(self,df_blank):
        """subtract blank from sample data to calibrate."""
        if self.wt == 'blank':
            self.df['Q_Corrected'] = self.df["HeatFlow(mW)"]
        elif self.wt != 'blank':
            self.df2['Q_Corrected'] = self.df2["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"]

    def correct_HF_smooth(self, df_blank):
        """for smoothed data. Subtract blank from sample data to calibrate."""
        df = self.df_mean.set_index('T(K)',drop=False)
        self.df_mean = self.df_mean.set_index('T(K)',drop=False)
        df_blank = df_blank.set_index('T(K)',drop=False)
        df_data, df_blank = df.align(df_blank)
        self.df_mean['Q_Corrected'] = df_data["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"]
        # ADD 0.01 mW DUE TO -5 ML OF AIR DIFFERENCE
        # self.df_mean['Q_Corrected'] = df_data["HeatFlow(mW)"] - df_blank["HeatFlow(mW)"] + .01
        self.df_mean['err'] = self.df_mean['std'] + df_blank['std']
        self.df_mean['err%'] = np.abs(self.df_mean['err'] / self.df_mean['Q_Corrected'] *100)
        self.df2 = self.df_mean[self.df_mean['Q_Corrected'].notna()]
        self.df2 = self.df2.drop(['T(K)','T_Mean(K)','std','std_%'],axis=1)
        self.df2 = self.df2.reset_index(drop=True)

    def remove_kink(self,lb,ub):
        # l_idx = self.df2.iloc[(self.df2['Sample Temperature(K)'] - lb).abs().argsort()[:1]].index.tolist()[0]
        # u_idx = self.df2.iloc[(self.df2['Sample Temperature(K)'] - ub).abs().argsort()[:1]].index.tolist()[0]
        self.df2 = self.df2[(self.df2['Sample Temperature(K)']<lb) | (self.df2['Sample Temperature(K)']>ub)]

    def mean_smoothing(self,df,T_bin,sd,save=0):

        # calc average value within each degree
        df3 = df.copy()
        df3['T_bin'] = round(df3['Sample Temperature(K)'])
        df4 = df3['T_bin'].copy().drop_duplicates()

        df_mean = pd.DataFrame()
        for _,T_val in df4.iteritems():
            df5 = df3[(df3['T_bin'] < T_val + T_bin) & (df3['T_bin'] > T_val - T_bin)]
            sr = df5.mean()
            sr['T(K)'] = T_val
            sr['std'] = np.std(df5['HeatFlow(mW)'])
            sr['std_%'] = np.abs(sr['std']/sr['HeatFlow(mW)']*100)
            sr = sr.drop(['T_bin'])
            df_mean = pd.concat([df_mean,sr],axis=1)
        df_mean = df_mean.transpose()
        df_mean = df_mean.reset_index(drop=True)
        df_mean['T_Mean(K)'] = df_mean['Sample Temperature(K)'].copy()
        df_mean['Sample Temperature(K)'] = df_mean['T(K)'].copy()
        df_mean = df_mean[1:-1]

        self.df_mean = df_mean
        if save ==1:
            df_cp_mean.to_csv(sd + f'mean_{T_bin}K_H2O-NH3_{self.m}g.csv', index=False)

    def plot_Q_T(self,fd,plt_err=0,save=0):
        plt.figure()
        plt.scatter(self.df2['Sample Temperature(K)'],self.df2['Q_Corrected'],1)
        if plt_err==1:
            if self.m==0:
                err_name = 'std'
            else:
                err_name ='err'
            plt.errorbar(self.df2['Sample Temperature(K)'],self.df2['Q_Corrected'],
                         yerr=(self.df2[err_name]),
                     linestyle='', ecolor='gray', capsize=1, capthick=0.5, elinewidth=0.5,zorder=-1)
        plt.xlabel("Temperature (K)"); plt.ylabel('Heat Flow (mW)')
        plt.title('{} wt% Thermogram'.format(self.wt))
        if self.wt == 'blank':
            plt.ylim([-2,2])
        if save == 1:
            plt.savefig(fd+'{}wt%_thermogram_{}g.png'.format(self.wt,self.m))
            plt.close()
        elif save ==0:
            show_plot_max()

    def calc_Fc(self):
        """calculate crystal fraction using Croft 1988 equation"""
        Fc = Fc_calc(self.df2['Sample Temperature(K)'], self.wt)
        Fc_err = abs(Fc_calc(self.df2['Sample Temperature(K)'] + T_err, self.wt) -
                     Fc_calc(self.df2['Sample Temperature(K)'] - 1e-4, self.wt)) / 2
        # Fc.index = np.arange(self.t_start, self.t_end, 1)
        # Fc_err.index = np.arange(self.t_start, self.t_end, 1)
        self.df2['Fc'] = Fc
        self.df2['Fc_err'] = Fc_err

    def plot_Fc_T(self,fd,save=0):
         plt.figure()
         plt.plot(self.df2['Sample Temperature(K)'], self.df2['Fc'])
         plt.xlabel('Temperature (K)'); plt.ylabel('Crystal Fraction')
         plt.title('{} wt% Crystal Fraction'.format(self.wt))
         plt.ylim([0, 1])
         if save == 1:
             plt.savefig(fd + '{}wt%_crystal_fraction_{}g.png'.format(self.wt,self.m))
             plt.close()
         elif save == 0:
             plt.show()

    def calc_liquidus(self):
        self.liquidus = self.df2.loc[self.df2['Fc'].idxmin()]

    def save_data(self,od):
        saveFile = od + '{}wt%_hf_prepped_{}g.csv'.format(self.wt,self.m)
        col_names = ['index', 'time(s)', 'furnanceT(K)', 'sampleT(K)', 'Q(mW)', 'Q_corrected(mW)', 'Q_err', 'Q_err%','Fc', 'Fc_err']
        self.df2.to_csv(saveFile, index=False, header=col_names)

    def save_liquidus(self,od):
        saveFile = od + '{}wt%_liquidus_{}g.txt'.format(self.wt,self.m)
        with open(saveFile, 'w') as fileOut:
            lines = self.liquidus.to_string(index=True, header=True)
            fileOut.writelines(lines)

##  #TODO class 2

def cp_2_cpm(df):
    df['cpm(J/molK)'] = df['cp(J/gK)'] * (17.031 * df['X'] + 18.0528 * (1 - df['X']))

class DataPrep:
    def __init__(self,mass,weight_pct):
        self.m = mass
        self.wt = weight_pct
        self.df = []
        self.df_pure = []
        self.df_cp_p = []
        self.df_cp_m = []
        self.df_hb = []


    def import_data(self, dd):
        """import hf prepped files"""
        self.df = pd.read_csv(dd + '{}wt%_hf_prepped_{}g.csv'.format(self.wt,self.m), header=0)

    def calib_data(self, ld,cd):
        """apply calibration from Ice Cp offset from FS06. Based off HF"""

        # import FS06 to get corresponding HF with data
        df_FW_raw = pd.read_csv(ld + 'ice_Cp.csv', names=['T', 'Cp'])
        f_FW_Cp = interp1d(df_FW_raw['T'], df_FW_raw['Cp'] / 1000,fill_value="extrapolate")
        cp = f_FW_Cp(self.df['sampleT(K)'])
        self.df['FW_Q(mW)'] = -(self.m * 1 * cp / (10 * 60) * 1000)

        # import calibration function coefficient
        dev_coeff = (pd.read_csv(cd + 'calib_coeff.csv', header=None)).to_numpy()
        self.df['calib_factor'] = np.polyval(dev_coeff,self.df['FW_Q(mW)'])
        self.df['Q_corrected(mW)'] = self.df['Q_corrected(mW)'] - self.df['calib_factor']

    def calc_cp_melt_nomix(self, iceLh_f, iceCp_f, iceCp_err_interp):
        """calculate specific heat of liquid in melt phase"""

        # extract data from melting phase (based on Fc)
        df_melt = self.df[self.df["Fc"].first_valid_index():self.df["Fc"].last_valid_index()].copy()

        # arrays to plot heat budget
        c_melt = np.zeros([len(df_melt), 1])
        c_melt[:] = np.nan
        c_heat = np.zeros([len(df_melt), 1])
        c_heat[:] = np.nan
        l_heat = np.zeros([len(df_melt), 1])
        l_heat[:] = np.nan
        l_mix = np.zeros([len(df_melt), 1])
        l_mix[:] = np.nan
        j_total = np.zeros([len(df_melt), 1])
        j_total[:] = np.nan
        j_sum = np.zeros([len(df_melt), 1])
        j_sum[:] = np.nan

        # calculation
        Cp = np.zeros([len(df_melt), 1])
        Cp[:] = np.nan
        Cp_err = np.zeros([len(df_melt), 1])
        Cp_T = np.zeros([len(df_melt), 1])
        hf_err = np.zeros([len(df_melt), 1])
        for i in range(len(df_melt) - 1):
            Q_integral = -(
                        df_melt['Q_corrected(mW)'].iloc[i + 1] / 1000 + df_melt['Q_corrected(mW)'].iloc[i] / 1000) * (
                                 10 * 60) / 2
            dFc = df_melt['Fc'].iloc[i] - df_melt['Fc'].iloc[i + 1]
            dT = df_melt['sampleT(K)'].iloc[i + 1] - df_melt['sampleT(K)'].iloc[i]

            Cp_T[i] = (df_melt['sampleT(K)'].iloc[i + 1] + df_melt['sampleT(K)'].iloc[i]) / 2
            Fc_now = (df_melt['Fc'].iloc[i + 1] + df_melt['Fc'].iloc[i]) / 2
            # cp = (total heat - heat to melt ice - heat to warm ice) / (mass * dT)
            c_melt[i] = dFc * self.m * iceLh_f(Cp_T[i])
            c_heat[i] = Fc_now * self.m * iceCp_f(Cp_T[i]) * dT

            Cp[i] = (Q_integral - c_melt[i] - c_heat[i]) \
                    / (self.m * (1 - df_melt['Fc'].iloc[i + 1]) * dT)

            Qint_err = (df_melt['Q_err'].iloc[i + 1] / 1000 + df_melt['Q_err'].iloc[i] / 1000) * (10 * 60) / 2
            c_melt_err = np.sqrt(
                (df_melt['Fc_err'].iloc[i + 1] + df_melt['Fc_err'].iloc[i] / dFc) ** 2 + (m_err / self.m) ** 2) * \
                         c_melt[i]  # TODO missing Latent heat uncertainty
            c_heat_err = np.sqrt(
                (df_melt['Fc_err'].iloc[i + 1] / Fc_now) ** 2 + (m_err / self.m) ** 2 + iceCp_err_interp(
                    df_melt['sampleT(K)'].iloc[i]) ** 2) * c_heat[i]
            hf_err[i] = Qint_err + c_melt_err + c_heat_err

            # hf_err = Q_err + np.sqrt((m_err/self.m)**2 + (df_melt['Fc_err'].iloc[i]/dFc)**2) * c_melt[i] + \
            #          np.sqrt((df_melt['Fc_err'].iloc[i]/df_melt['Fc'].iloc[i])**2 + iceCp_err_interp(df_melt['sampleT(K)'].iloc[i])**2 +
            #                  m_err**2/self.m**2 + T_err**2/dT**2)*c_heat[i]
            # Cp_err[i] = np.sqrt((hf_err*17.7) ** 2 / Q_integral ** 2 + m_err ** 2 / self.m ** 2 + (2*T_err) ** 2 / dT ** 2)
            Cp_err[i] = np.sqrt((hf_err[i] / (Q_integral - c_melt[i] - c_heat[i])) ** 2 +
                                (m_err / self.m) ** 2 + (2 * T_err) ** 2 / dT ** 2) * 100

            # heat budget calculations
            j_total[i] = Q_integral
            l_heat[i] = self.m * (1 - df_melt['Fc'].iloc[i + 1]) * dT * Cp[i]
            j_sum[i] = c_melt[i] + c_heat[i] + l_heat[i]

        # calculate mass fraction and molality of data
        Xl_calc_v = np.vectorize(Xl_calc)
        Cp_X = Xl_calc_v(Cp_T)
        Cp_m = Cp_X / 17.031 / ((1 - Cp_X) / 1000)

        # output cp in a DF
        self.df_cp_m = pd.DataFrame(
            {'T(K)': np.ravel(Cp_T), 'Fc': np.ravel(df_melt['Fc']), 'X': np.ravel(Cp_X), 'm(mol/kg)': np.ravel(Cp_m),
             'cp(J/gK)': np.ravel(Cp), 'cp_err(%)': np.ravel(Cp_err)})
        self.df_cp_m = self.df_cp_m.dropna()
        cp_2_cpm(self.df_cp_m)

        # output heat budget in a DF
        self.df_hb = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'j_total': np.ravel(j_total), 'j_sum': np.ravel(j_sum),
                                   'c_melt': np.ravel(c_melt), 'c_heat': np.ravel(c_heat), 'l_heat': np.ravel(l_heat)})
        self.df_hb = self.df_hb.dropna()

    def calc_cp_melt(self,iceLh_f,iceCp_f,iceCp_err_interp,mix_coeff):
        """calculate specific heat of liquid in melt phase"""
        
        # extract data from melting phase (based on Fc)
        df_melt = self.df[self.df["Fc"].first_valid_index():self.df["Fc"].last_valid_index()].copy()

        # arrays to plot heat budget
        c_melt = np.zeros([len(df_melt), 1]); c_melt[:] = np.nan
        c_heat = np.zeros([len(df_melt), 1]); c_heat[:] = np.nan
        l_heat = np.zeros([len(df_melt), 1]); l_heat[:] = np.nan
        l_mix = np.zeros([len(df_melt), 1]); l_mix[:] = np.nan
        j_total = np.zeros([len(df_melt), 1]); j_total[:] = np.nan
        j_sum = np.zeros([len(df_melt), 1]); j_sum[:] = np.nan
        
        # calculation
        Cp = np.zeros([len(df_melt), 1]); Cp[:] = np.nan
        Cp_err = np.zeros([len(df_melt), 1])
        Cp_T = np.zeros([len(df_melt), 1])
        hf_err = np.zeros([len(df_melt), 1])
        for i in range(len(df_melt) - 1):
            Q_integral = -(df_melt['Q_corrected(mW)'].iloc[i + 1] / 1000 + df_melt['Q_corrected(mW)'].iloc[i] / 1000) * (
                    10*60) / 2
            dFc = df_melt['Fc'].iloc[i] - df_melt['Fc'].iloc[i + 1]
            dT = df_melt['sampleT(K)'].iloc[i + 1] - df_melt['sampleT(K)'].iloc[i]

            Cp_T[i] = (df_melt['sampleT(K)'].iloc[i + 1] + df_melt['sampleT(K)'].iloc[i]) / 2
            Fc_now = (df_melt['Fc'].iloc[i+1] + df_melt['Fc'].iloc[i])/2
            # cp = (total heat - heat to melt ice - heat to warm ice - heat to mix) / (mass * dT)
            c_melt[i] = dFc * self.m * iceLh_f(Cp_T[i])
            c_heat[i] = Fc_now  * self.m * iceCp_f(Cp_T[i]) * dT

            Hmix = np.polyval(mix_coeff, Xl_calc(df_melt['sampleT(K)'].iloc[i + 1])) - np.polyval(mix_coeff, Xl_calc(df_melt['sampleT(K)'].iloc[i]))

            l_mix[i] = dFc * self.m * Hmix

            Cp[i] = (Q_integral - c_melt[i] - c_heat[i] - l_mix[i]) \
                    / (self.m * (1 - df_melt['Fc'].iloc[i + 1]) * dT)

            Qint_err = (df_melt['Q_err'].iloc[i+1] /1000 + df_melt['Q_err'].iloc[i]/1000)*(10*60)/2
            c_melt_err = np.sqrt((df_melt['Fc_err'].iloc[i+1]+df_melt['Fc_err'].iloc[i]/dFc)**2 + (m_err/self.m)**2) * c_melt[i] #TODO missing Latent heat uncertainty
            c_heat_err = np.sqrt((df_melt['Fc_err'].iloc[i+1]/Fc_now)**2 + (m_err/self.m)**2 + iceCp_err_interp(df_melt['sampleT(K)'].iloc[i])**2) * c_heat[i]
            hf_err[i] = Qint_err + c_melt_err + c_heat_err

            # hf_err = Q_err + np.sqrt((m_err/self.m)**2 + (df_melt['Fc_err'].iloc[i]/dFc)**2) * c_melt[i] + \
            #          np.sqrt((df_melt['Fc_err'].iloc[i]/df_melt['Fc'].iloc[i])**2 + iceCp_err_interp(df_melt['sampleT(K)'].iloc[i])**2 +
            #                  m_err**2/self.m**2 + T_err**2/dT**2)*c_heat[i]
            # Cp_err[i] = np.sqrt((hf_err*17.7) ** 2 / Q_integral ** 2 + m_err ** 2 / self.m ** 2 + (2*T_err) ** 2 / dT ** 2)
            Cp_err[i] = np.sqrt((hf_err[i]/(Q_integral - c_melt[i] - c_heat[i]))**2 +
                                (m_err/self.m) ** 2 + (2*T_err) ** 2 / dT ** 2) * 100

            # heat budget calculations
            j_total[i] = Q_integral
            l_heat[i] = self.m * (1 - df_melt['Fc'].iloc[i + 1]) * dT * Cp[i]
            j_sum[i] = c_melt[i] + c_heat[i] + l_heat[i] + l_mix[i]

        # calculate mass fraction and molality of data
        Xl_calc_v = np.vectorize(Xl_calc)
        Cp_X = Xl_calc_v(Cp_T)
        Cp_m = Cp_X / 17.031 / ((1 - Cp_X) / 1000)

        # output cp in a DF
        self.df_cp_m = pd.DataFrame(
            {'T(K)': np.ravel(Cp_T), 'Fc': np.ravel(df_melt['Fc']), 'X': np.ravel(Cp_X), 'm(mol/kg)': np.ravel(Cp_m),
             'cp(J/gK)': np.ravel(Cp), 'cp_err(%)': np.ravel(Cp_err)})
        self.df_cp_m = self.df_cp_m.dropna()

        # apply calibration correction #TODO
        # self.df_cp_m['cp(J/gK)'] = self.df_cp_m['cp(J/gK)'] * (1-calib)
        # self.df_cp_m['cp_err(%)'] = self.df_cp_m['cp_err(%)'] * (1-calib)
        # self.df_cp_m = self.df_cp_m.merge(df_calib,on='T(K)')
        # self.df_cp_m['cp(J/gK)'] = self.df_cp_m['cp(J/gK)'] * (1 - self.df_cp_m['dev%']/100)
        # self.df_cp_m['cp_err(%)'] = self.df_cp_m['cp_err(%)'] * (1 - self.df_cp_m['dev%']/100)



        cp_2_cpm(self.df_cp_m)

        # output heat budget in a DF
        self.df_hb = pd.DataFrame({'T(K)': np.ravel(Cp_T),'j_total': np.ravel(j_total), 'j_sum': np.ravel(j_sum),
                                   'c_melt':np.ravel(c_melt),'c_heat':np.ravel(c_heat),'l_heat':np.ravel(l_heat),'l_mix':np.ravel(l_mix)})
        self.df_hb = self.df_hb.dropna()

    def calc_cp_pure(self,range_arg):
        """calculate specific heat of liquid in pure liquid phase"""

        # for wt% > 0, only calculate in region after melt.
        # for wt% = 0, calculate whole segment (pure ice)
        if self.wt != 0:
            df_pure = self.df[self.df["Fc"].last_valid_index()+1:self.df.index[-1]].copy()
        elif self.wt == 0:
            df_pure = self.df.copy()
        
        Cp = np.zeros([len(df_pure), 1]); Cp[:] = np.nan
        Cp_err = np.zeros([len(df_pure), 1])
        Cp_T = np.zeros([len(df_pure), 1])
        for i in range(len(df_pure) - 1):
            # heatflow per temperature step (HARDCODED TO 0.1K/MIN SETTING BECAUSE INPUT IS 1 K STEPS)
            Q_integral = -(df_pure['Q_corrected(mW)'].iloc[i + 1] / 1000 + df_pure['Q_corrected(mW)'].iloc[i] / 1000) * (
                    10*60) / 2
            # temperature step
            dT = df_pure['sampleT(K)'].iloc[i + 1] - df_pure['sampleT(K)'].iloc[i]

            Cp[i] = Q_integral / (self.m * dT) # calculate cp
            Cp_err[i] = np.sqrt(((df_pure['Q_err'].iloc[i+1] /1000 + df_pure['Q_err'].iloc[i]/1000)*(10*60)/2)**2/Q_integral**2 +
                                m_err**2/self.m**2 + (2*T_err)**2/dT**2) *100

            # Cp_err[i] = np.sqrt(
            #     ((1/self.m*dT)**2)*(Q_err*(df_pure['time(s)'].iloc[i + 1] - df_pure['time(s)'].iloc[i]))**2 +
            # ((Q_integral/(self.m**2)*dT)**2*m_err**2) +
            # (Q_integral/self.m*dT**2)**2*(2*T_err)**2) / Cp[i] *100

            Cp_T[i] = (df_pure['sampleT(K)'].iloc[i + 1] + df_pure['sampleT(K)'].iloc[i]) / 2 #set T as in between step

        Cp_m = self.wt / 100 / 17.031 / ((1 - self.wt / 100) / 1000)
        self.df_cp_p = pd.DataFrame({'T(K)': np.ravel(Cp_T), 'X': self.wt / 100, 'm': Cp_m, 'cp(J/gK)': np.ravel(Cp),
                                     'cp_err(%)': np.ravel(Cp_err)})
        self.df_cp_p = self.df_cp_p.dropna()

        # apply calibration correction #TODO
        # self.df_cp_p['cp(J/gK)'] = self.df_cp_p['cp(J/gK)'] * (1 - calib)
        # self.df_cp_p['cp_err(%)'] = self.df_cp_p['cp_err(%)'] * (1 - calib)
        # self.df_cp_p = self.df_cp_p.merge(df_calib,on='T(K)')
        # self.df_cp_p['cp(J/gK)'] = self.df_cp_p['cp(J/gK)'] * (1 - self.df_cp_p['dev%']/100)
        # self.df_cp_p['cp_err(%)'] = self.df_cp_p['cp_err(%)'] * (1 - self.df_cp_p['dev%']/100)

        # remove in between calculations due to kink removed
        if range_arg==1:
            self.df_cp_p = self.df_cp_p[(self.df_cp_p['T(K)'] < 294) | (self.df_cp_p['T(K)'] > 297)]
            if self.m == 3.7107:
                self.df_cp_p = self.df_cp_p[(self.df_cp_p['T(K)'] < 290) | (self.df_cp_p['T(K)'] > 297)]
            if self.m == 3.8153:
                self.df_cp_p = self.df_cp_p[(self.df_cp_p['T(K)'] < 294) | (self.df_cp_p['T(K)'] > 304)]
        cp_2_cpm(self.df_cp_p)

        # manual cutting for ice cp
        if self.wt == 0:
            # self.df_cp_p = self.df_cp_p[0:-750].copy()
            self.df_cp_p = self.df_cp_p[9:159].copy()


    def save_cp_melt(self,od):
        saveFile = od + '{}wt%_cp_melt_{}g.csv'.format(self.wt,self.m)
        self.df_cp_m.to_csv(saveFile, index=False)

    def save_cp_pure(self, od):
        saveFile = od + '{}wt%_cp_pure_{}g.csv'.format(self.wt,self.m)
        self.df_cp_p.to_csv(saveFile, index=False)

    def plot_cp(self, cp_type, df_melt, df_pure, fd, save=0, ylim=None):
        if cp_type == 'cp':
            cp_col = 'cp(J/gK)'
            y_label = 'Specific Heat (J $g^{-1}$ $K^{-1}$)'
        elif cp_type == 'cpm':
            cp_col = 'cpm(J/molK)'
            y_label = 'Molar Specific Heat (J $mol^{-1}$ $K^{-1}$)'
        plt.figure()
        plt.scatter(df_melt['T(K)'], df_melt[cp_col], 0.1, label='Melt Phase')
        plt.scatter(df_pure['T(K)'], df_pure[cp_col], 0.1, label='Liquid Phase')
        plt.xlabel('T (K)'); plt.ylabel(y_label)
        plt.title('{}wt% {}'.format(self.wt,cp_type))
        plt.ylim(ylim)
        plt.legend(markerscale=5)
        if save == 1:
            if ylim:
                plt.savefig(fd + '{}wt%_{}-T_zoom_{}g.png'.format(self.wt,cp_type,self.m))
            else:
                plt.savefig(fd + '{}wt%_{}-T_{}g.png'.format(self.wt,cp_type,self.m))
            plt.close()
        elif save == 0:
            show_plot_max()

    def plot_hb(self,fd,od,l_bound,u_bound,save=0):

        l_idx = self.df_hb.iloc[(self.df_hb['T(K)']-l_bound).abs().argsort()[:1]].index.tolist()[0]
        u_idx = self.df_hb.iloc[(self.df_hb['T(K)']-u_bound).abs().argsort()[:1]].index.tolist()[0]
        df_hb_c = self.df_hb[l_idx:u_idx+1].copy()

        if self.m == 4.5858 and len(self.df_hb)>1800:
            df_hb_c = df_hb_c.drop(range(1600,1800))


        #normal heat budget
        plt.figure()
        plt.axhline(y=0, color='k',linestyle='--')
        plt.plot(df_hb_c['T(K)'], df_hb_c['j_total'], 'k', label='total recorded heat')
        plt.plot(df_hb_c['T(K)'], df_hb_c['c_melt'], 'r', label='melt ice')
        plt.plot(df_hb_c['T(K)'], df_hb_c['c_heat'], 'm', label='heat ice fraction')
        plt.plot(df_hb_c['T(K)'], df_hb_c['l_mix'], 'g', label='mix molten ice')
        plt.plot(df_hb_c['T(K)'], df_hb_c['l_heat'], 'b', label='heat liquid fraction')
        plt.plot(df_hb_c['T(K)'], df_hb_c['j_sum'], 'y--', label='sum of heat processes')
        # plt.ylim([-1, 2])
        plt.xlabel('Temperature (K)')
        plt.ylabel(r'Energy (J)')
        plt.title('{} wt% Energy Budget'.format(self.wt))
        plt.legend(fontsize=5)
        if save == 1:
            saveFig = fd + '{}wt%_energy_budget_{}g.png'.format(self.wt,self.m)
            plt.savefig(saveFig)
            plt.close()
        elif save == 0:
            plt.show()

        # Stacked heat budget
        plt.figure()
        plt.stackplot(df_hb_c['T(K)'], df_hb_c['c_melt'] / df_hb_c['j_total'],
                      df_hb_c['c_heat'] / df_hb_c['j_total'],
                      df_hb_c['l_mix'] / df_hb_c['j_total'],
                      df_hb_c['l_heat'] / df_hb_c['j_total'],
                      labels=['melt ice', 'heat ice fraction', 'mix molten ice', 'heat liquid fraction'])
        plt.xlabel('T (K)')
        # plt.ylim([0, 1])
        plt.title('{} wt% Energy Budget (Normalized)'.format(self.wt))
        plt.legend(fontsize=5, frameon=1, loc='lower right')
        if save == 1:
            saveFig = fd + '{}wt%_cp_energy_budget_normalized_{}g.png'.format(self.wt,self.m)
            df_hb_c.to_csv(od +'{}wt%_hb_{}g.csv'.format(self.wt, self.m),index=False)
            plt.savefig(saveFig)
            plt.close()
        elif save == 0:
            plt.show()

    def cut_cp(self,df,l_bound,u_bound,od,fd,plt_err=1,test=1):

        def plot_fx(plt_err):
            plt.scatter(df['T(K)'], df['cp(J/gK)'], 1, label = 'Cut')
            plt.scatter(df_c['T(K)'], df_c['cp(J/gK)'], 1, label = 'Kept')
            if plt_err == 1:
                plt.errorbar(df_c['T(K)'], df_c['cp(J/gK)'], yerr=(df_c['cp_err(%)']/100*df_c['cp(J/gK)']),
                             linestyle='', ecolor='gray', capsize=1, capthick=0.5, elinewidth=0.5,zorder=-1)

            plt.xlabel('Temperature (K)')
            plt.ylabel(r'Specific Heat (J $g^{-1}$ $K^{-1}$)')
            if df.equals(self.df_cp_p):
                plt.ylim([3.4, 4.5])
            elif df.equals(self.df_cp_m):
                plt.ylim([2.4, 4.3])
            # plt.ylim([])

        l_idx = df.iloc[(df['T(K)']-l_bound).abs().argsort()[:1]].index.tolist()[0]
        u_idx = df.iloc[(df['T(K)']-u_bound).abs().argsort()[:1]].index.tolist()[0]
        df_c = df[l_idx:u_idx+1].copy()

        if self.m == 4.5858 and len(df)>1800:
            df_c = df_c.drop(range(1600,1800))

        plt.figure()
        plt.subplot(1,2,1)
        plot_fx(plt_err)
        plt.title('{} wt% lower'.format(self.wt))
        # plt.xlim([176, 195])
        plt.xlim([270, 284])

        plt.subplot(1,2,2)
        plot_fx(plt_err)
        plt.title('{} wt% upper'.format(self.wt))
        # plt.xlim([df['T(K)'].iloc[-500],df['T(K)'].iloc[-1]])

        if test ==1:
            show_plot_max()
            print('ye')
        else:
            if df['X'].iloc[0]*100 < self.wt + 0.01:
                cp_type = 'liquid_phase'
            else:
                cp_type = 'along_liquidus'
            plt.figure()
            plot_fx(plt_err)
            plt.title('{} wt% Specific Heat'.format(self.wt))
            plt.legend()

            plt.savefig(fd +'{}wt%_cp_{}_{}g.png'.format(self.wt, cp_type,self.m))
            plt.close()

            df_c.to_csv(od +'{}wt%_cp_{}_{}g.csv'.format(self.wt, cp_type,self.m),index=False)

    def cut_savgol(self,df,l_bound,od,fd,plot=0):
        """Use savgol filter to smoothen thermogram. Then, calculate derivative of filtered data, and apply another
         round of filtering. then, using the last index closest to zero to determine turning point of thermogram, as
         upper bound"""

        def plot_fx():
            plt.scatter(df['T(K)'], df['cp(J/gK)'], 0.1, label = 'Cut')
            plt.scatter(df_cut['T(K)'], df_cut['cp(J/gK)'], 0.1, label = 'Kept')
            plt.xlabel('Temperature (K)')
            plt.ylabel(r'Specific Heat (J $g^{-1}$ $K^{-1}$)')
            plt.ylim([2, 5])

        l_idx = df.iloc[(df['T(K)']-l_bound).abs().argsort()[:1]].index.tolist()[0]
        # u_idx = df.iloc[(df['T(K)']-u_bound).abs().argsort()[:1]].index.tolist()[0]
        df_c = df[l_idx:].copy()

        if self.m == 4.5858 and len(df)>1800:
            df_c = df_c.drop(range(1600,1800))

        # spline and sma testing stuff
        # spl = UnivariateSpline(df_c['T(K)'], df_c['cp(J/gK)'], k=2)
        # sma_idx_step = round(sma_step / (df['T(K)'][1] - df['T(K)'][0]))
        # sma = df_c['cp(J/gK)'].rolling(sma_idx_step, center=False).mean()
        # plt.plot(df_c['T(K)'], spl(df_c['T(K)']), c='tab:green', label='Spline')
        # plt.plot(df_c['T(K)'], sma, c='tab:red', label='SMA')
        # plt.plot(df_c['T(K)'].iloc[0:-1]+.5,np.diff(sma,n=1))

        savgol = savgol_filter(df_c['cp(J/gK)'],11,1)
        savgol_diff = savgol_filter(np.diff(savgol, n=1), 11, 2)
        sg_diff = pd.Series(savgol_diff)

        # find last index of value closest to zero but not negative
        sg_a = sg_diff - 0
        sg_b = sg_a[sg_a>0].index[-1]

        df_cut = df[l_idx:sg_b+l_idx].copy()

        if self.m == 4.5858:
            df_cut = df_cut.reset_index(drop=True)
            df_cut = df_cut.drop(df_cut.index[36:])

        if df['X'].iloc[0] * 100 < self.wt + 0.01:
            cp_type = 'pure'
        else:
            cp_type = 'melt'

        plt.figure()
        plot_fx()
        plt.title('{} wt% data cut savgol'.format(self.wt))
        plt.legend(markerscale=5)

        plt.savefig(fd + '{}wt%_cp_cut_{}_savgol_{}g.png'.format(self.wt, cp_type,self.m))
        plt.close()

        df_cut.to_csv(od + '{}wt%_cp_cut_{}_{}g.csv'.format(self.wt, cp_type,self.m), index=False)

        if plot ==1:
            plt.subplot(2, 1, 1)
            plt.scatter(df_c['T(K)'], df_c['cp(J/gK)'], 0.1, label='Data')
            plt.plot(df_c['T(K)'], savgol, c='tab:purple', label='Savgol')
            plt.axvline(x=df_c['T(K)'].iloc[sg_b], c='k')
            plt.ylabel(r'Cp (J $g^{-1}$ $K^{-1}$)')
            plt.ylim([3, 5])
            plt.xlim([176, 270])
            plt.title('{} wt% savgol fit and 1st order derivative'.format(self.wt))

            plt.subplot(2, 1, 2)
            plt.plot(df_c['T(K)'].iloc[0:-1] + .5, savgol_diff)
            plt.axhline(y=0, c='k')
            plt.axvline(x=df_c['T(K)'].iloc[sg_b], c='k')
            plt.xlabel('Temperature (K)')
            plt.ylabel(r'dSMA')
            plt.xlim([176, 270])
            plt.ylim([-.1, .1])

            # show_plot_max()
            plt.savefig('t_spline_cut/' + '{}wt%_cp_melt_savgol_derivative_{}g.png'.format(self.wt,self.m))
            plt.close()

##  #TODO class 3

class DataCP:
    def __init__(self,m,weight_pct):
        self.m = m
        self.wt = weight_pct
        self.df_m = []
        self.df_p = []
        self.df_m_mean = []
        self.df_p_mean = []
        self.spl_m = []
        self.spl_p = []
        self.savgol_m = []
        self.savgol_p = []
        self.linear_m = []
        self.coef_m = []
        self.linear_p = []
        self.shomate_m = []
        self.shomate_p = []
        self.sd_m = []
        self.sd_p = []

    def import_data(self, dd, mean=0):
        """import cp cut files"""
        self.df_m = pd.read_csv(dd + '{}wt%_cp_cut_melt_{}g.csv'.format(self.wt,self.m), header=0)
        self.df_p = pd.read_csv(dd + '{}wt%_cp_cut_pure_{}g.csv'.format(self.wt,self.m), header=0)
        if mean == 1:
            self.df_m_mean = pd.read_csv(dd + '{}wt%_cp_cut_melt_mean_{}g.csv'.format(self.wt,self.m), header=0)
            self.df_p_mean = pd.read_csv(dd + '{}wt%_cp_cut_pure_mean_{}g.csv'.format(self.wt,self.m), header=0)

    def mean_cp(self,regime,T_bin,od):

        if regime == 'melt':
            df = self.df_m
        elif regime == 'pure':
            df = self.df_p

        # calc average value within each degree
        df2 = df.copy()
        df2['T_bin'] = round(df2['T(K)'])
        df3 = df2['T_bin'].copy().drop_duplicates()

        df_cp_mean = pd.DataFrame()
        for _,T_val in df3.iteritems():
            df4 = df2[(df2['T_bin'] < T_val+T_bin) & (df2['T_bin'] > T_val-T_bin)]
            sr = df4.mean()
            sr['T(K)'] = T_val
            sr['std'] = np.std(df4['cp(J/gK)'])
            sr['std_%'] = sr['std']/sr['cp(J/gK)']*100
            if regime == 'melt':
                sr['X'] = Xl_calc(T_val)
            sr = sr.drop(['T_bin'])
            df_cp_mean = pd.concat([df_cp_mean,sr],axis=1)
        df_cp_mean = df_cp_mean.transpose()
        df_cp_mean = df_cp_mean.reset_index()
        df_cp_mean = df_cp_mean.drop(columns=['index'])
        df_cp_mean = df_cp_mean[1:-1]
        df_cp_mean = df_cp_mean[(df_cp_mean['T(K)'] <294) | (df_cp_mean['T(K)'] >297)]

        if regime == 'melt':
            self.df_m_mean = df_cp_mean
        elif regime == 'pure':
            self.df_p_mean = df_cp_mean

        df_cp_mean.to_csv(od + '{}wt%_cp_cut_{}_mean_{}g.csv'.format(self.wt, regime,self.m), index=False)

    def calc_SMA(self,regime,T_bin,od):
        if regime == 'melt':
            df = self.df_m
        elif regime == 'pure':
            df = self.df_p

        # calc SMA based on T_bin
        df2 = df.iloc[1:100].copy()
        # determines number of indexes to form one Tbin (eg. 1C, 5C), to be used as rolling average bin width
        T_bin_idx = abs(abs(df2['T(K)'].iloc[0] - df2['T(K)'])-1).argsort()[1]
        df[f'cp_SMA_{T_bin}K'] = df['cp(J/gK)'].rolling(T_bin_idx, center=True).mean()
        df.to_csv(od + '{}wt%_cp_cut_{}_{}g.csv'.format(self.wt, regime,self.m), index=False)

    def calc_spline(self):
        self.spl_m = UnivariateSpline(self.df_m['T(K)'], self.df_m['cp(J/gK)'], k=3)
        self.spl_p = UnivariateSpline(self.df_p['T(K)'], self.df_p['cp(J/gK)'], k=3)

    def calc_savgol(self):
        self.savgol_m = savgol_filter(self.df_m['cp(J/gK)'], 301, 1)
        self.savgol_p = savgol_filter(self.df_p['cp(J/gK)'], 301, 1)

    def calc_linear(self):
        self.coef_m,cov = np.polyfit(self.df_m['T(K)'],self.df_m['cp(J/gK)'],1,cov=True)
        self.linear_m = np.poly1d(self.coef_m)
        self.sd_m = np.sqrt(np.diag(cov))

        coef = np.polyfit(self.df_p['T(K)'], self.df_p['cp(J/gK)'], 1)
        self.linear_p = np.poly1d(coef)

    def calc_shomate(self, data):
        if data == 'mean':
            T_m = self.df_m_mean['T(K)']
            df_m = self.df_m_mean['cp(J/gK)']
            T_p = self.df_p_mean['T(K)']
            df_p = self.df_p_mean['cp(J/gK)']
        elif data == 'SMA':
            T_m = self.df_m['T(K)'].loc[np.isnan(self.df_m['cp_SMA_1K']) == False]
            df_m = self.df_m['cp_SMA_1K'].loc[np.isnan(self.df_m['cp_SMA_1K']) == False]
            T_p = self.df_p['T(K)'].loc[np.isnan(self.df_p['cp_SMA_1K']) == False]
            df_p = self.df_p['cp_SMA_1K'].loc[np.isnan(self.df_p['cp_SMA_1K']) == False]
        elif data == 'all':
            T_m = self.df_m['T(K)']
            df_m = self.df_m['cp(J/gK)']
            sigma_m = self.df_m['cp_err(%)']/100 * self.df_m['cp(J/gK)']
            T_p = self.df_p['T(K)']
            df_p = self.df_p['cp(J/gK)']
            sigma_p = self.df_p['cp_err(%)']/100 * self.df_p['cp(J/gK)']

            self.shomate_m, pcov = curve_fit(base.quad_eqn,T_m,df_m,sigma=sigma_m,p0=[0,0,0])
            self.sd_m = np.sqrt(np.diag(pcov))
            self.shomate_p, pcov = curve_fit(base.shomate_eqn, T_p, df_p, sigma=sigma_p,p0=[0, 0, 0, 0, 0])
            self.sd_p = np.sqrt(np.diag(pcov))


    def plot_data(self,ax, colour,marker,idx,raw_data=1, shomate=1,errors=1,ci=1,melt=1,pure=1):
        # manually control where labels dependent on melt/pure arguments
        # fig = plt.figure()
        # ax = plt.gca()
        def plot_them(shom_label=0,ls='-'):
            if raw_data == 1:
                ax.scatter(df['T(K)'], df['cp(J/gK)'], 5, marker=marker, c=colour, linewidths=0, zorder=15)
            if shomate == 1:
                Trange = np.linspace(df['T(K)'][0], df['T(K)'].iloc[-1], 500)
                linestuff = ax.plot(Trange, eqn(Trange, *coeff), ls, linewidth=1, c=colour, zorder=10,alpha=1)
                # linestuff = ax.plot(Trange, shomalite_eqn(Trange, *coeff), ls, linewidth=0.75, c=colour, zorder=10)
                if raw_data == 0 and shom_label ==1:
                    linestuff[0].set_label(f'{self.wt} wt%')
            if errors == 1 and raw_data==1:
                ax.errorbar(df['T(K)'], df['cp(J/gK)'], yerr=(df['cp_err(%)'] / 100 * df['cp(J/gK)']),
                             linestyle='', ecolor=colour, capsize=1, capthick=0.5, elinewidth=0.5,
                             errorevery=(idx * 2, 16), zorder=5)

        if melt==1 and pure ==1:
            df = self.df_p
            eqn = shomate_eqn
            coeff = self.shomate_p
            plot_them()

        else:
            if melt == 1:
                df = self.df_m
                eqn = quad_eqn
                coeff = self.shomate_m
                plot_them()
            if pure ==1:
                df = self.df_p
                eqn = shomate_eqn
                coeff = self.shomate_p
                plot_them()

    def plot_PI(self,ax,colour,melt=1,pure=1,for_bf=0):
        import scipy.stats as stats

        def plot_it():
            y_model = eqn(df['T(K)'], *p)
            x = df['T(K)']
            y = df['cp(J/gK)']

            n = df['T(K)'].size
            m = p.size
            dof = n - m
            t = stats.t.ppf(0.975, n - m)

            # Estimates of Error in Data/Model
            resid = y - y_model  # residuals; diff. actual data from predicted values
            chi2 = np.sum((resid / y_model) ** 2)  # chi-squared; estimates error in data
            chi2_red = chi2 / dof  # reduced chi-squared; measures goodness of fit
            s_err = np.sqrt(np.sum(resid ** 2) / dof)  # standard deviation of the error

            x2 = np.linspace(df['T(K)'][0], df['T(K)'].iloc[-1], 500)
            y2 = eqn(x2, *p)
            # plot_ci_manual(t, s_err, n, x, x2, y2, colour,ax=ax)


            pi = t * s_err * np.sqrt(1 + 1 / n + (x2 - np.mean(x)) ** 2 / np.sum((x - np.mean(x)) ** 2))
            ax.fill_between(x2, eqn(x2, *p) + pi, eqn(x2, *p) - pi, color=colour, alpha=.2, linewidth=0,zorder=-10)
            # ax.plot(x2, eqn(x2, *p) - pi, "--", color=colour)
            # ax.plot(x2, eqn(x2, *p) + pi, "--", color=colour)

        if pure ==1:
            df = self.df_p
            p = self.shomate_p
            eqn = shomate_eqn
            plot_it()
        if melt ==1:
            df = self.df_m
            p = self.shomate_m
            eqn = quad_eqn
            plot_it()

    def plot_fit(self,colour,melt=1,pure=1):
        # manually control where labels are
        if melt ==1:
            plt.scatter(self.df_m['T(K)'], self.df_m['cp(J/gK)'], 0.05, c=colour, marker='x', alpha=0.05)
            # plt.plot(self.df_m['T(K)'], self.linear_m(self.df_m['T(K)']), c=colour, label='{} wt%'.format(self.wt))
            plt.plot(self.df_m['T(K)'], shomate_eqn(self.df_m['T(K)'],*self.shomate_m), c=colour, label='{} wt%'.format(self.wt))
            # plt.plot(self.df_m['T(K)'], self.spl_m(self.df_m['T(K)']), c=colour, label='{} wt%'.format(self.wt))
            # plt.plot(self.df_m['T(K)'], self.savgol_m, '--',  c=colour)
        if pure ==1:
            plt.scatter(self.df_p['T(K)'], self.df_p['cp(J/gK)'], 0.05, c=colour, marker='x', alpha=0.05)
            # plt.scatter(self.df_p['T(K)'], self.df_p['cp(J/gK)'], 0.05, c=colour, marker='x', alpha=0.5,label='{} wt%'.format(self.wt))
            # plt.plot(self.df_p['T(K)'], self.linear_p(self.df_p['T(K)']), c=colour)#,label='{} wt%'.format(self.wt))
            plt.plot(self.df_p['T(K)'], shomate_eqn(self.df_p['T(K)'], *self.shomate_p), c=colour, label='{} wt%'.format(self.wt))
            # plt.plot(self.df_p['T(K)'], self.spl_p(self.df_p['T(K)']), c=colour)#,label='{} wt%'.format(self.wt))
            # plt.plot(self.df_p['T(K)'], self.savgol_p, '--', c=colour)

    def plot_spl(self,colour,melt=1,pure=1):
        # manually control where labels are
        if melt ==1:
            plt.scatter(self.df_m['T(K)'], self.df_m['cp(J/gK)'], 0.05, c=colour, marker='x', alpha=0.05)
            plt.plot(self.df_m['T(K)'], self.spl_m(self.df_m['T(K)']), c=colour, label='{} wt%'.format(self.wt))
        if pure ==1:
            plt.scatter(self.df_p['T(K)'], self.df_p['cp(J/gK)'], 0.05, c=colour, marker='x', alpha=0.05)
            plt.plot(self.df_p['T(K)'], self.spl_p(self.df_p['T(K)']), c=colour)#,label='{} wt%'.format(self.wt))

    def plot_spl2(self,melt=1,pure=1):
        # manually control where labels are
        if melt ==1:
            plt.scatter(self.df_m['T(K)'], self.df_m['cp(J/gK)'], 0.05, marker='x', alpha=0.05)
            plt.plot(self.df_m['T(K)'], self.spl_m(self.df_m['T(K)']), label='{:.1f} wt%'.format(self.wt))
        if pure ==1:
            plt.scatter(self.df_p['T(K)'], self.df_p['cp(J/gK)'], 0.05, marker='x', alpha=0.05)
            plt.plot(self.df_p['T(K)'], self.spl_p(self.df_p['T(K)']))#,lab