"""
created by Bing Hong CHUA 29Sep22

script objective:
from cut specific heat data, apply various spline and moving average fits for further analysis
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import base

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
fd = r"t_all_cp/"

## SCRIPT ARGUMENTS
ramp_rate = 0.1
wt_ls = [5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912] # based off liquidus alignment - 0.5K

mass_ls = [4.5386,4.1943,4.5858,4.5202,3.8153,3.7107,3.7778]

# colour_ls = ['#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84','#081d58']
# colour_ls = ['#d73027','#fc8d59','#fee090','#f2f26d','#e0f3f8','#91bfdb','#4575b4']
colour_ls = ['#d73027','#f46d43','#fdae61','#e8e884','#abd9e9','#74add1','#4575b4']
marker_ls = ['o','s','d','P','X','^','v']
## PLOT ALL CP ON SAME PLOT. OPTIONS FOR SF AND TRF EOS


def overall_fun(axis,pure,lit,water,inset,melt,best_fit,comp,data_arg,conc):
    pure_arg = pure
    lit_arg = lit
    water_arg = water
    inset_arg = inset
    melt_arg = melt
    best_fit_arg = best_fit
    comp_arg = comp
    conc_arg = conc
    axis = axis

    # plt.close()
    # fig = plt.figure('combined')
    # axis = fig.add_subplot(111)

    def fun(ax):
        coef_num = 3
        coef_arr = np.empty((len(wt_ls),coef_num))
        sd_arr = np.empty((len(wt_ls),coef_num))
        dfbf = pd.DataFrame()

        for idx, wt in enumerate(wt_ls):
            # idx = 1 #TODO for single data debugging
            # wt = wt_ls[idx] #TODO for single data debugging
            colour = colour_ls[idx]
            marker = marker_ls[idx]
            m = mass_ls[idx]

            data = base.DataCP(m,wt)
            data.import_data(dd,mean=0)

            if conc_arg == 1:
                data.df_m['T(K)'] = data.df_m['X']*100

            data.calc_shomate('all')
            # data.calc_shomate('all',eqn=base.shomalite_eqn)

            # save coef and sd of polynomial for best fit calculation of all fits
            if best_fit_arg==1:
                coef_arr[idx] = data.shomate_m
                sd_arr[idx] = data.sd_m
            elif best_fit_arg==2:
                dfbf = pd.concat ([dfbf,data.df_m])

            # plt.figure('combined')
            data.plot_data(ax, colour, marker, idx, raw_data=data_arg, shomate=1,errors=1,melt=melt_arg, pure=pure_arg)

            if not (melt == 1 and pure == 1) and data_arg==0:
                data.plot_PI(ax,colour,melt_arg,pure_arg)

            if lit_arg == 1:
                Trange=np.arange(210,320,.1)
                from scipy.interpolate import interp1d
                df_St = pd.read_csv(ld + 'Baptiste_cp.csv', header=0)
                TRF_interp = interp1d(df_St['T(K)'], df_St[str(wt)])
                ax.plot(Trange,TRF_interp(Trange),linestyle=(0, (3, 1, 1, 1)),c=colour,zorder=-5)

        if lit_arg ==1:
            # IMPORT AND PLOT WREWSKY AND CHERNENKAYA
            df_wrewsky = pd.read_csv(ld + 'Baptiste/' + 'wrewsky_1924.csv', header=0)
            df_wrewsky = df_wrewsky.copy().iloc[:-3]
            df_chernenkaya = pd.read_csv(ld + 'Baptiste/' + 'chernenkaya_1971.csv', header=0)

            c_id = [1.4,5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912,40.1]
            c_val = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#e8e884', '#abd9e9', '#74add1', '#4575b4','#313695']

            norm = plt.Normalize(min(c_id), max(c_id))
            tuples = list(zip(map(norm, c_id), c_val))
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)


            plt.scatter(df_wrewsky['T(K)'], df_wrewsky['cp(J/kgK)'] / 1000, 8, marker='P', edgecolors='k',
                        linewidths=.2, c=df_wrewsky['massfrac']*100, cmap=cmap,
                        label='Wrewsky & Kaigorodoff 1924')
            plt.scatter(df_chernenkaya['T(K)'], df_chernenkaya['cp(J/kgK)'] / 1000, 8, marker='X', edgecolors='k',
                        linewidths=.2,c=df_chernenkaya['massfrac']*100, cmap=cmap,
                        label='Chernenkaya 1971')

        if best_fit_arg==2:
            from scipy.optimize import curve_fit
            bf_coeff, pcov = curve_fit(base.quad_eqn, dfbf['T(K)'], dfbf['cp(J/gK)'], sigma=dfbf['cp_err(%)']/100 * dfbf['cp(J/gK)'], p0=[0, 0, 0])
            T_range = np.arange(184.5, 233.5, 0.1)

            fit, = ax.plot(T_range, base.quad_eqn(T_range, *bf_coeff), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=1,
                           zorder=5, label='overall fit')

            base.plot_PI_melt(dfbf,base.quad_eqn,bf_coeff,ax,'k')
            coef_mean = bf_coeff

        elif best_fit_arg==1:
            # calculate overall best fit
            weights_arr = (1 / sd_arr ** 2)
            coef_mean = np.sum(coef_arr * weights_arr, 0) / np.sum(weights_arr, 0)
            sd_mean = np.sqrt(1 / np.sum(weights_arr, 0))
            # plot the overall best fit
            if conc_arg ==1:
                T_range = np.linspace(0.31231*100, 0.20561*100, 500)
            else:
                T_range = np.arange(184.5, 233.5,0.1)

            # df_bfcp = pd.DataFrame()
            # df_bfcp['T(K)'] = T_range
            # func = np.vectorize(base.Xl_calc)
            # df_bfcp['massfrac'] = func(T_range)
            # df_bfcp['cp(J/gK)'] = base.quad_eqn(T_range,*coef_mean)
            # df_bfcp.to_csv(od + 'z_liquidus_cp.csv', index=False)

            fit, = ax.plot(T_range, base.quad_eqn(T_range,*coef_mean),'k',linestyle=(0, (5, 1)),linewidth=1, zorder=5,label = 'overall fit')

        if melt_arg ==1 and pure_arg ==1:
            fit.set_label('$C_p$ along liquidus')

            # ADD IN CHAN !!!!
            df_chan = pd.read_csv(ld + 'Baptiste/' + 'chan_1964.csv', header=0)
            this = df_chan.copy().iloc[0]

            c_id = [1.4, 5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912, 40.1]
            c_val = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#e8e884', '#abd9e9', '#74add1', '#4575b4',
                     '#313695']

            norm = plt.Normalize(min(c_id), max(c_id))
            tuples = list(zip(map(norm, c_id), c_val))
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

            plt.scatter(this['T(K)'], this['cp(J/kgK)'] / 1000, 8, c=matplotlib.colors.rgb2hex(cmap(this['massfrac']*100)), marker='s', edgecolors='k',
                        linewidths=.2,label='33 wt% Chan & Giauque 1964')

            # HARDCODE TICK MARKS AT 30, 27, 24, 21 WT%
            val_ls = [30, 27, 24, 21]
            for val in val_ls:
                func = np.vectorize(base.Xl_calc)
                df = pd.DataFrame(func(T_range))
                out_df = base.get_closest_df_value(df, 0, val / 100, 1)
                T_set = T_range[out_df.index[:]]
                length = .03
                plt.vlines(T_set, base.quad_eqn(T_set, *coef_mean) - length, base.quad_eqn(T_set, *coef_mean) + length,
                           linewidth=.5, color='k')
                plt.text(T_set-3, base.quad_eqn(T_set, *coef_mean) + length * 1.1, f'{val} wt%', fontsize=3.2,
                         horizontalalignment = 'center', verticalalignment='bottom')
        if water_arg ==1:
            df_water = pd.read_csv(ld + 'water_Cp_IAPWS.csv', names=['T(K)', 'Cp'])
            ax.plot(df_water['T(K)'], df_water['Cp'], linewidth=0.75, c='k', linestyle=(0, (5, 1)), zorder=-20,label='$H_2O$ (0 wt%) IAPWS-95')
            df_water = pd.read_csv(ld + 'water_Cp_Bol.csv',names=['T(K)','Cp'])
            ax.plot(df_water['T(K)'],df_water['Cp'],linewidth=0.75,c='gray',linestyle=(0, (5, 1)),zorder=-20,label='$H_2O$ (0 wt%) B19')

        if lit_arg==1:
            ax.plot(0, 0, 'k', linestyle=(0, (3, 1, 1, 1)), label='Tillner-Roth & Friend 1998')


        if comp_arg == 1 and pure_arg==0:
            # plot secondary axis of composition

            def Xl_calc100(T):
                """Croft equation for liquidus of water-rich ammonia-water mixutre, up till peritectic (X<=0.329)
                input: current tempertaure
                output: ammonia weight%"""
                lqd_coeff = [56695, -46269, 11842, -1651.4, -53.07, 273.15 - T]  # CROFT 1987
                roots = np.roots(lqd_coeff)  # get roots of the polynomial
                roots_real = np.real(roots[np.iscomplex(roots) == False])  # remove imaginary roots
                root_actual = roots_real[(0 < roots_real) & (roots_real <= 0.329)]  # keep roots within peritectic
                if len(root_actual) == 0:
                    root_actual = np.nan
                else:
                    root_actual = root_actual[0]
                return root_actual*100

            def T_calc100(Xl):
                """Croft equation for liquidus of water-rich ammonia-water mixutre, up till peritectic (X<=0.329)
                input: ammonia weight%
                output: Temperature"""
                T = 273.15 - 53.07 * (Xl/100) - 1651.4 * (Xl/100) ** 2 + 11842 * (Xl/100) ** 3 - 46269 * (Xl/100) ** 4 + 56695 * (Xl/100) ** 5
                return T


            Xl_calc_v = np.vectorize(Xl_calc100)
            plt.tick_params(axis='x', which='both', top=False)
            secax = ax.secondary_xaxis('top', functions=(Xl_calc_v, T_calc100))
            secax.set_xlabel('$NH_3$ Mass Fraction along liquidus (wt%)')
            plt.xlim([180,240])

    fun(axis)

    # plot the rest
    if melt_arg == 1 and pure_arg ==0:
        plt.xlabel('Temperature along liquidus (K)');
    else:
        a=1
        # plt.xlabel('Temperature (K)');
    # plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')

    # legend = plt.legend(prop={'size': 4},title_fontsize =4,loc='lower right')

    # if data_arg ==1:
    #     # legend.set_title("starting mass fraction")

    # if  lit_arg ==1 or (melt_arg ==1 and pure_arg ==1):
    #     legend = plt.legend(prop={'size': 4},title="starting mass fraction",title_fontsize =4,loc='lower right')

    if melt_arg == 1 and pure_arg == 1:
        figname = 'all'
        # plt.ylim([3.12, 4.65])
        plt.ylim([2.8, 4.4])
        plt.xlim([180, 315])
    elif melt_arg == 1 and pure_arg != 1:
        figname = 'melt'
        plt.ylim([2.58, 3.95])
        plt.xlim([182,236])
    elif melt_arg != 1 and pure_arg == 1:

        figname = 'pure'
        if lit_arg ==1:
            figname = figname + '_lit'
        plt.ylim([3.4, 4.53])
        plt.xlim([210, 315])

    if inset_arg == 1 and melt_arg ==0:
        # inset axes....
        # axins = axis.inset_axes([0.57, 0.1, 0.4, 0.4]) #bottom right
        axins = axis.inset_axes([0.09, 0.56, 0.315, 0.315]) #top left
        fun(axins)
        # subregion of the original image
        x1, x2, y1, y2 = 266, 293, 4.16, 4.245
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.tick_params(axis='y', labelsize=4)
        axins.tick_params(axis='x', labelsize=4)
        axins.xaxis.set_major_locator(plt.MaxNLocator(3))
        axins.yaxis.set_major_locator(plt.MaxNLocator(2))

        in_w = 0.3

        indi = axis.indicate_inset_zoom(axins, edgecolor="black",linewidth=in_w)
        for a,b in enumerate(indi[1]):
            b.set(linewidth=in_w)


    # base.show_plot_max()

    print('Finished running script')

## ABOVE LIQUIDUS

fig, ax = plt.subplots(2,1,sharey=True,sharex=True)

plt.tight_layout()

overall_fun(ax[0],1,0,1,1,0,0,0,1,0)
overall_fun(ax[1],1,1,1,0,0,0,0,0,0)

# label and legend
base.splt_axis_label(fig,"Temperature (K)",'Specific Heat (J $g^{-1}$ $K^{-1}$)')
ax[0].text(0.025, 0.95, 'a)', horizontalalignment='left',verticalalignment='top', transform=ax[0].transAxes)
ax[1].text(0.025, 0.95, 'b)', horizontalalignment='left',verticalalignment='top', transform=ax[1].transAxes)

ax[0].plot(250, 4, linewidth=0.75, c=colour_ls[0], label='Shomate Fit')
# ax[0].scatter(250,4, 8, marker='P', edgecolors='k',linewidths=.2, c='#a50026',label='Wrewsky & Kaigorodoff 1924',alpha=0)
# ax[0].scatter(250,4, 8, marker='X', edgecolors='k',linewidths=.2, c='#a50026',label='Chernenkaya 1971',alpha=0)
# ax[0].plot(250, 4, linestyle=(0, (3, 1, 1, 1)),c=colour_ls[0], label='Tillner-Roth & Friend 1998')
ax[0].scatter(250,4, 8, marker='P', edgecolors='k',linewidths=.2, c='#a50026',label='WK24',alpha=0)
ax[0].scatter(250,4, 8, marker='X', edgecolors='k',linewidths=.2, c='#a50026',label='C71',alpha=0)
ax[0].plot(250, 4, linestyle=(0, (3, 1, 1, 1)),c=colour_ls[0], label='TF98')

handles, labels = ax[0].get_legend_handles_labels()
order = [0,1,2,3,4,5,6,9,7,8,12,10,11]

leg = ax[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],
bbox_to_anchor=(1.01, 1),loc='upper left',prop={'size': 4})
for lh in leg.legendHandles:
    lh.set_alpha(1)

base.show_plot_max()
plt.savefig(fd + 'above_liquidus.png')

##
# overall_fun(0,0,0,0,1,2,1,1,0)
# overall_fun(1,0,1,0,1,2,0,0,0)

# overall_fun(0,0,0,0,1,1,0,1,1)
## debug end line
print('Finished running script')
