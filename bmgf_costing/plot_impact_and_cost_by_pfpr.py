"""
Jaline Gerardin
Nov 2018

Plot impact of ATSB, IRS, and next-gen ITN interventions in 8-site setup.
Expects output generated by atsb_llin_impact_analyzer.

Plotters:
- Cases averted relative to 1) no interventions and 2) diminished efficacy pyrthroid nets, as a function of initial
PfPR2-10.
- Cost per case averted relative to diminished efficacy pyrethroid nets. User can supply multiple per-device costs.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
mpl.rcParams['pdf.fonttype'] = 42

from plotting.colors import load_color_palette

projectdir = os.path.join('E:/', 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'atsb')
datadir = os.path.join(projectdir, 'sim_data')
plotdir = os.path.join(projectdir, 'sim_plots')


def load_sim_df(data_fname):

    df = pd.read_csv(data_fname)
    df = df.groupby(['Site_Name', 'x_Temporary_Larval_Habitat', 'intervention']).agg(np.mean).reset_index()
    df['cases per 1000'] = df['New_Clinical_Cases']/df['Population']*1000
    df['infections per 1000'] = df['New_Clinical_Cases']/df['Population']*1000

    df = df.sort_values(by=['Site_Name', 'x_Temporary_Larval_Habitat'])

    return df


def plot_cases_averted(df, name, baseline_channel, ignore_list, savename) :

    num_interventions = len(df['intervention'].unique()) - len(ignore_list)

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    fig = plt.figure(name, figsize=(15,5*(num_interventions/3)))
    axes = [fig.add_subplot(num_interventions, 3, d + 1) for d in range(3*(num_interventions))]
    fig.subplots_adjust(left=0.05, right=0.98)
    palette = load_color_palette()

    for s, (site, sdf) in enumerate(df.groupby('Site_Name')):
        sdf = sdf.sort_values(by='PfPR2to10')
        for d, datachannel in enumerate(['cases per 1000', 'infections per 1000', 'Blood_Smear_Parasite_Prevalence']):
            baseline = sdf[sdf['intervention'] == baseline_channel]
            for i, (intervention, idf) in enumerate(sdf[~(sdf['intervention'].isin(ignore_list))].groupby('intervention')) :
                diff = [x - y for x, y in zip(baseline[datachannel], idf[datachannel])]
                ax = axes[i*3+d]
                xvar = [y for x, y in zip(diff, baseline['PfPR2to10'].values) if (y > 0 and x >= 0)]
                yvar = [x for x, y in zip(diff, baseline['PfPR2to10'].values) if (y > 0 and x >= 0)]
                ys = lowess(yvar, xvar, frac=0.2)[:,1]
                ax.plot(xvar, ys, '-', color=palette[s], label=site)
                # ax.scatter(xvar, yvar, 5, color=palette[s])
                ax.set_title('%s %s' % (intervention, name))
                ax.set_ylabel('%s averted' % datachannel)
                if i == num_interventions-1 :
                    ax.set_xlabel('initial PfPR2to10')
                if 'per' in datachannel :
                    ax.set_ylim(-50, 4000)
                else :
                    ax.set_ylim(-0.002, 0.4)

    axes[-1].legend()
    fig.savefig(os.path.join(plotdir, '%s.png' % savename))
    fig.savefig(os.path.join(plotdir, '%s.pdf' % savename), format='PDF')
    plt.close(fig)


def plot_cost_per_case_averted(df, datachannel, savename, basechannel):

    costs = { 'itn' : [1.85, 2.13, 2.3], # per person
              'llin' : [2.25, 2.5, 3], # per person
              'irs' : [3, 4.8, 10], # per person
              'atsb' : [1, 3, 5, 7],
              } # per device
    palettes = ['Blues', 'Reds', 'Greens', 'Purples']

    # accounts for population, coverage, and distribution frequency. assume household size=5
    cost_scale_factors = {
        'itn' : 2000*0.6/2,
        'llin' : 2000*0.6/2,
        'irs' : 2000*0.6*3,
        'atsb' : 2000/5*12*0.6,
    }

    fig = plt.figure(figsize=(10, 8))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.07, top=0.95)

    for s, (site, sdf) in enumerate(df.groupby('Site_Name')):

        sdf = sdf.sort_values(by='PfPR2to10')
        ax = fig.add_subplot(3,3,s+1)
        baseline = sdf[sdf['intervention'] == basechannel]

        sdf = sdf[~(sdf['intervention'].isin(['none', 'itn']))]
        for i, (intervention, idf) in enumerate(sdf.groupby('intervention')):
            intervention_type = intervention.split('_')[0]

            minidf = pd.DataFrame( { 'baseline PfPR2to10' : baseline['PfPR2to10'].values,
                                     'baseline_%s' % datachannel : baseline[datachannel].values,
                                     '%s_%s' % (intervention, datachannel) : idf[datachannel].values})
            minidf['diff'] = minidf['baseline_%s' % datachannel] - minidf['%s_%s' % (intervention, datachannel)]
            mdf = minidf[(minidf['diff'] >= 0) & (minidf['baseline PfPR2to10'] > 0)]
            palette = sns.color_palette(palettes[i], len(costs[intervention_type]))

            for c, single_cost in enumerate(costs[intervention_type]):
                cost = single_cost*cost_scale_factors[intervention_type]
                xvar = mdf['baseline PfPR2to10'].values
                yvar = [cost/x for x in mdf['diff'].values]
                ys = lowess(yvar, xvar, frac=0.2)[:,1]
                ax.plot(xvar, ys,
                        '-', color=palette[c], label='%s %.2f' % (intervention, single_cost))
        ax.set_yscale('log')
        ax.set_ylim(1e-1, 1e3)
        if s > 4:
            ax.set_xlabel('initial PfPR2to10')
        if s == 3:
            ax.set_ylabel('dollars per %s averted' % datachannel)
        ax.set_title(site)
        if s == 7:
            ax.legend()
    fig.savefig(os.path.join(plotdir, '%s_cost_per_%s_averted_v_%s.png' % (savename, datachannel, basechannel)))
    fig.savefig(os.path.join(plotdir, '%s_cost_per_%s_averted_v_%s.pdf' % (savename, datachannel, basechannel)),
                format='PDF')
    plt.close(fig)


if __name__ == '__main__' :

    expt_name = "atsb_llin_HS_v1v2"
    data_fname = os.path.join(datadir, "%s.csv" % expt_name)

    df = load_sim_df(data_fname)

    savename = '%s_cases_averted_by_site' % expt_name
    plot_cases_averted(df, 'baseline', 'none', ['none'], '%s_v_baseline' % savename)
    plot_cases_averted(df, 'itn', 'itn', ['none', 'itn'], '%s_v_itn' % savename)

    basechannel = 'itn'
    interventions = ['none', 'itn', 'atsb_cdc', 'atsb_hlc', 'irs_180']
    sdf = df[df['intervention'].isin(interventions)]
    plot_cost_per_case_averted(sdf, 'New_Clinical_Cases', '%s_part1' % expt_name, basechannel)
    plot_cost_per_case_averted(sdf, 'New_Infections', '%s_part1' % expt_name, basechannel)

    interventions = ['none', 'itn', 'llin_no_disc', 'llin']
    sdf = df[df['intervention'].isin(interventions)]
    plot_cost_per_case_averted(sdf, 'New_Clinical_Cases', '%s_part2' % expt_name, basechannel)
    plot_cost_per_case_averted(sdf, 'New_Infections', '%s_part2' % expt_name, basechannel)

    basechannel = 'none'
    interventions = ['none', 'itn', 'llin_no_disc', 'llin']
    sdf = df[df['intervention'].isin(interventions)]
    plot_cost_per_case_averted(sdf, 'New_Clinical_Cases', '%s_part2' % expt_name, basechannel)
    plot_cost_per_case_averted(sdf, 'New_Infections', '%s_part2' % expt_name, basechannel)

    interventions = ['none', 'itn', 'atsb_alone_hlc', 'atsb_alone_cdc']
    sdf = df[df['intervention'].isin(interventions)]
    plot_cost_per_case_averted(sdf, 'New_Clinical_Cases', '%s_part3' % expt_name, basechannel)
    plot_cost_per_case_averted(sdf, 'New_Infections', '%s_part3' % expt_name, basechannel)

    plt.show()
