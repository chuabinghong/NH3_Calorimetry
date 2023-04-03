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
fd = r"o_finalPlots/"

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

            # if not (melt == 1 and pure == 1):
            #     data.plot_PI(ax,colour,melt_arg,pure_arg)

            if lit_arg == 1:
                Trange=np.arange(210,320,.1)
                from scipy.interpolate import interp1d
                df_St = pd.read_csv(ld + 'TF98_cp_selected.csv', header=0)
                TRF_interp = interp1d(df_St['T(K)'], df_St[str(wt)])
                ax.plot(Trange,TRF_interp(Trange),linestyle=(0, (3, 1, 1, 1)),c=colour,zorder=-5)

        if lit_arg ==1:
            # IMPORT AND PLOT WREWSKY AND CHERNENKAYA
            df_wrewsky = pd.read_csv(ld + 'lit_cp/' + 'wrewsky_1924.csv', header=0)
            df_wrewsky = df_wrewsky.copy().iloc[:-3]
            df_chernenkaya = pd.read_csv(ld + 'lit_cp/' + 'chernenkaya_1971.csv', header=0)
            df_prieto = pd.read_csv(ld + 'lit_cp/' + 'prieto_2022_noLiq.csv', header=0)

            import matplotlib as mpl
            c_id = [1.4,5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912,40.1]
            c_val = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#e8e884', '#abd9e9', '#74add1', '#4575b4','#313695']

            norm = plt.Normalize(min(c_id), max(c_id))
            tuples = list(zip(map(norm, c_id), c_val))
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

            plt.scatter(df_wrewsky['T(K)'], df_wrewsky['cp(J/kgK)'] / 1000, 15, marker='P', edgecolors='k',
                        linewidths=.5, c=df_wrewsky['massfrac']*100, cmap=cmap,vmin=c_id[0], vmax=c_id[-1])
            plt.scatter(df_chernenkaya['T(K)'], df_chernenkaya['cp(J/kgK)'] / 1000, 15, marker='X', edgecolors='k',
                        linewidths=.5,c=df_chernenkaya['massfrac']*100, cmap=cmap,vmin=c_id[0], vmax=c_id[-1])
            plt.scatter(df_prieto['T(K)'], df_prieto['cp(J/kgK)'] / 1000, 15, marker='p', edgecolors='k',
                        linewidths=.5, c=df_prieto['massfrac'] * 100, cmap=cmap,vmin=c_id[0], vmax=c_id[-1])


        if best_fit_arg==2:
            from scipy.optimize import curve_fit
            bf_coeff, pcov = curve_fit(base.quad_eqn, dfbf['T(K)'], dfbf['cp(J/gK)'], sigma=dfbf['cp_err(%)']/100 * dfbf['cp(J/gK)'], p0=[0, 0, 0])
            T_range = np.arange(184.5, 233.5, 0.1)

            fit, = ax.plot(T_range, base.quad_eqn(T_range, *bf_coeff), 'k', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=1,
                           zorder=5, label='Overall Fit')

            base.plot_PI_melt(dfbf,base.quad_eqn,bf_coeff,ax,'gray')


            # df_bfcp = pd.DataFrame()
            # df_bfcp['T(K)'] = T_range
            # func = np.vectorize(base.Xl_calc)
            # df_bfcp['massfrac'] = func(T_range)
            # df_bfcp['cp(J/gK)'] = base.quad_eqn(T_range,*bf_coeff)
            # df_bfcp.to_csv(od + 'z_liquidus_cp_2.csv', index=False)

            if pure_arg==0:
                # ADD IN CHAN !!!!
                df_chan = pd.read_csv(ld + 'lit_cp/' + 'chan_1964.csv', header=0)
                this = df_chan.copy().iloc[0]

                c_id = [1.4, 5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912, 40.1]
                c_val = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#e8e884', '#abd9e9', '#74add1', '#4575b4',
                         '#313695']

                norm = plt.Normalize(min(c_id), max(c_id))
                tuples = list(zip(map(norm, c_id), c_val))
                cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

                plt.scatter(this['T(K)'], this['cp(J/kgK)'] / 1000, 15, c=matplotlib.colors.rgb2hex(cmap(this['massfrac']*100)), marker='s', edgecolors='k',
                            linewidths=.5,label='33 wt% CG64')

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



            fit, = ax.plot(T_range, base.quad_eqn(T_range, *coef_mean), 'k', linestyle=(0, (5, 1)), linewidth=1,
                           zorder=5, label='overall fit')

        if melt_arg ==1 and pure_arg ==1:
            fit.set_label('$C_p$ along liquidus')

            # ADD IN CHAN !!!!
            df_chan = pd.read_csv(ld + 'lit_cp/' + 'chan_1964.csv', header=0)
            this = df_chan.copy().iloc[0]

            c_id = [1.4, 5.2, 8.2, 8.4, 10.0, 14.3, 20.07, 26.912, 40.1]
            c_val = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#e8e884', '#abd9e9', '#74add1', '#4575b4',
                     '#313695']

            norm = plt.Normalize(min(c_id), max(c_id))
            tuples = list(zip(map(norm, c_id), c_val))
            cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

            # plt.scatter(this['T(K)'], this['cp(J/kgK)'] / 1000, 8, c=matplotlib.colors.rgb2hex(cmap(this['massfrac']*100)), marker='s', edgecolors='k',
            #             linewidths=.2)

            # HARDCODE TICK MARKS AT 30, 27, 24, 21 WT%
            val_ls = [30, 27, 24, 21]
            coef_mean = bf_coeff
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
            ax.plot(df_water['T(K)'], df_water['Cp'], linewidth=1, c='k', linestyle=(0, (5, 1)), zorder=-20,label='0 wt% IAPWS-95')
            df_water = pd.read_csv(ld + 'water_Cp_Bol.csv',names=['T(K)','Cp'])
            ax.plot(df_water['T(K)'],df_water['Cp'],linewidth=1,c='gray',linestyle=(0, (5, 1)),zorder=-20,label='0 wt% B19')


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
        plt.ylim([2.55, 4.53])
        plt.xlim([180, 316])
    elif melt_arg == 1 and pure_arg != 1:
        figname = 'melt'
        plt.ylim([2.49, 3.87])
        plt.xlim([182,236])
    elif melt_arg != 1 and pure_arg == 1:

        figname = 'pure'
        if lit_arg ==1:
            figname = figname + '_lit'
        plt.ylim([3.4, 4.53])
        plt.xlim([210, 316])

    if inset_arg == 1 and melt_arg ==0:
        # inset axes....
        axins = axis.inset_axes([0.57, 0.13, 0.38, 0.35]) #bottom right
        # axins = axis.inset_axes([0.09, 0.56, 0.315, 0.315]) #top left
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

    print('Finished one plot')


def label_one(layout):
    # labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Temperature (K)")
    plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
    ax[0].text(0.025, 0.98, 'a)', horizontalalignment='left', verticalalignment='top', transform=ax[0].transAxes)
    ax[1].text(0.025, 0.98, 'b)', horizontalalignment='left', verticalalignment='top', transform=ax[1].transAxes)

    # legend
    for i, wt in enumerate(wt_ls):
        wt_str = wt
        if wt == 20.07:
            wt_str = 20.1
        elif wt == 26.912:
            wt_str = 26.9
        plt.scatter(240, 4, 5,marker=marker_ls[i], c=colour_ls[i], linewidths=0, alpha=0, label=f'{wt_str} wt%')
    plt.plot(250, 4, linewidth=1, c=colour_ls[0], label='Shomate Fit')
    plt.plot(250, 4, linewidth=1, c='k', linestyle=(0, (5, 1)), label='0 wt% IAPWS-95')
    plt.plot(250, 4, linewidth=1, c='gray', linestyle=(0, (5, 1)), label='0 wt% B19')
    # plt.scatter(250,4, 15, marker='P', edgecolors='k',linewidths=.2, c='#a50026',label='Wrewsky & Kaigorodoff 1924',alpha=0)
    # plt.scatter(250,4, 8, marker='X', edgecolors='k',linewidths=.2, c='#a50026',label='Chernenkaya 1971',alpha=0)
    # plt.plot(250, 4, linestyle=(0, (3, 1, 1, 1)),c=colour_ls[0], l qabel='Tillner-Roth & Friend 1998')
    plt.scatter(250, 4, 15, marker='P', edgecolors='k', linewidths=.5, c='#a50026', label='WK24', alpha=0)
    plt.scatter(250, 4, 15, marker='X', edgecolors='k', linewidths=.5, c='#a50026', label='C71', alpha=0)
    plt.scatter(250, 4, 15, marker='p', edgecolors='k', linewidths=.5, c='#a50026', label='P22', alpha=0)
    plt.plot(250, 4, linestyle=(0, (3, 1, 1, 1)), c=colour_ls[0], label='TF98')

    handles, labels = plt.gca().get_legend_handles_labels()
    # order = [0,1,2,3,4,5,6,9,7,8,12,10,11]
    # order = [2,3,4,5,6,7,8,9,0,1,12,10,11]

    if layout == 'vert':
        bbox = (.5, -.16)
        ncol = 5
        loc = 'upper center'

    if layout == 'hori':
        bbox = (1.02,.5)
        ncol = 1
        loc = 'center left'
    leg = plt.legend(bbox_to_anchor=bbox, loc=loc, prop={'size': 4}, ncol=ncol,
                     frameon=True, edgecolor='k', framealpha=1)
    leg.get_frame().set_linewidth(0.5)

    for lh in leg.legendHandles:
        lh.set_alpha(1)


def label_two():
    plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
    for i, wt in enumerate(wt_ls):
        wt_str = wt
        if wt == 20.07:
            wt_str = 20.1
        elif wt == 26.912:
            wt_str = 26.9
        plt.scatter(240, 4, marker=marker_ls[i], c=colour_ls[i], linewidths=0, alpha=0, label=f'{wt_str} wt%')
    plt.plot(250, 4, linewidth=1, c=colour_ls[0], label='Quadratic Fit')
    plt.plot(250, 4, "--", linewidth=.5, color='gray', label='Prediction Limit')

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 10]
    # order = [2,3,4,5,6,7,8,9,0,1,12,10,11]

    leg = plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                     prop={'size': 4}, title='Starting Mass Fraction', title_fontsize=5, ncols=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
## ABOVE LIQUIDUS

layout = 'hori'

if layout == 'vert':
    w = 4
    fig, ax = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(w, w * 2 / 3))
elif layout == 'hori':
    w = 3.3
    fig, ax = plt.subplots(2, 1, sharey=True, sharex=True, figsize=(w, w * 3 / 3))
# fig, ax = plt.subplots(1,2,sharey=True,sharex=True,figsize=(3.20666667, 1.91166667))
plt.tight_layout()

overall_fun(ax[0],1,0,1,1,0,0,0,1,0)
overall_fun(ax[1],1,1,1,0,0,0,0,0,0)
label_one(layout)
plt.tight_layout()

# base.show_plot_max()
plt.savefig(fd + f'cp_liquidPhase.png')

##
w = 3.3
fig, ax = plt.subplots(1,1,sharey=True,sharex=True,figsize=(w,w*2/3))
overall_fun(ax,0,0,0,0,1,2,1,1,0,)
label_two()
# base.show_plot_max()
plt.savefig(fd + 'cp_alongLiquidus.png')

##

w = 3.3
fig, ax = plt.subplots(1,1,sharey=True,sharex=True,figsize=(w,w*2/3))
overall_fun(ax,1,0,1,0,1,2,0,0,0,)

plt.xlabel("Temperature (K)")
plt.ylabel('Specific Heat (J $g^{-1}$ $K^{-1}$)')
for i, wt in enumerate(wt_ls):
    wt_str = wt
    if wt == 20.07:
        wt_str = 20.1
    elif wt == 26.912:
        wt_str = 26.9
    plt.plot(240, 4, linewidth=1, c=colour_ls[i],alpha=0, label=f'{wt_str} wt%')

plt.plot(250, 4, "--", linewidth=.5, color='gray', label='Prediction Limit')

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,10,3,4,5,6,7,8,9,1,2]

leg = plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                     prop={'size': 4}, ncols=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)
plt.savefig(fd + 'cp_all.png')





## debug end line
print('Finished running script')
