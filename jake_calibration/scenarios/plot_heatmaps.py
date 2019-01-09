import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from scipy import interpolate
import os
import time
import copy
mpl.rcParams['pdf.fonttype'] = 42

userpath = 'E:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
plotdir = os.path.join(wdir, 'sim_plots')
datadir = os.path.join(wdir, 'sim_data')


def old_heatmap_plotter(df, channels, xvar, yvar, expname) :

    for channel in channels :

        df = df.groupby([xvar, yvar, 'Run_Number'])[channel].agg(np.mean).reset_index()
        df = df.groupby([xvar, yvar])[channel].agg(np.mean).reset_index()
        rdf = df[df[xvar] == 0]
        rdf = rdf.rename(columns={channel: '%s_ref' % channel})

        df = pd.merge(left=df, right=rdf[[yvar, '%s_ref' % channel]], on=yvar)
        if channel == 'Blood Smear Parasite Prevalence' :
            df = df[(df['%s_ref' % channel] >= 0.1) & (df['%s_ref' % channel] <= 0.47)]
        df['%s_reduction' % channel] = (df['%s_ref' % channel] - df[channel]) / df['%s_ref' % channel]
        df.loc[df['%s_reduction' % channel] < 0, '%s_reduction' % channel] = 0

        plot_heatmaps(df, channel, xvar, '%s_ref' % channel, '%s_reduction' % channel, expname)


def plot_heatmaps(df, channel, xvar, yvar, zvar, plotname, ymin=0.1, ymax=0.47) :

    df = df[(df[yvar] >= ymin) & (df[yvar] <= ymax)]

    x = np.linspace(np.min(df[xvar]), np.max(df[xvar]), 1000)
    y = np.linspace(np.min(df[yvar]), np.max(df[yvar]), 1000)
    xx, yy = np.meshgrid(x, y)

    sns.set_style('white', {'axes.linewidth': 0.5})
    fig = plt.figure('ATSB Infection Reduction Heatmap', figsize=(10,5))
    ax = fig.add_subplot(1,2,1)

    ax.set_xlabel(xvar)
    ax.set_ylabel(yvar)

    sdf = df.sample(n=1000)

    inter = interpolate.Rbf(sdf[xvar],
                            sdf[yvar],
                            sdf[zvar],
                            function='linear', smooth=0.05)
    zz = inter(xx, yy)
    cmap = 'RdYlBu'
    plt.imshow(zz, cmap = plt.get_cmap(cmap, 10), vmin=0, vmax=1)
    ax.invert_yaxis()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.colorbar()
    ax.set_title(channel)

    ax = fig.add_subplot(1,2,2)
    ax.set_xlabel(xvar)
    ax.set_ylabel(yvar)
    palette = sns.color_palette(cmap, 101)
    c = plt.contour(xx, yy, zz, cmap=plt.get_cmap((cmap)), vmin=0, vmax=1,
                    levels=np.arange(0,1.05,0.1))
    # plt.clabel(c, inline=1, fontsize=8, fmt='%.2f')
    ax.scatter(sdf[xvar], sdf[yvar], 20, color=[palette[int(x * 100)] for x in sdf[zvar].values], alpha=0.3)
    ax.set_xlim(0,0.25)
    ax.set_ylim(ymin, ymax)

    plt.savefig(os.path.join(plotdir, '%s %s.png' % (plotname, channel)))
    plt.savefig(os.path.join(plotdir, '%s %s.pdf' % (plotname, channel)), format='PDF')
    # plt.show()
    plt.close('all')


def plot_summary(summary_fname):

    xvar = 'killing'
    yvar = 'x_Temporary_Larval_Habitat'
    binchannel = "PfPR by Parasitemia and Age Bin"
    binnames = ['true', 'PCR', 'RDT']
    ybounds = { 'RDT PfPR by Parasitemia and Age Bin' : [0.05, 0.4],
                'true PfPR by Parasitemia and Age Bin' : [0.05, 1],
                'PCR PfPR by Parasitemia and Age Bin' : [0.05, 0.85],
                "PfPR by Age Bin" : [0.1, 0.47],
                "Annual Clinical Incidence by Age Bin" : [0.5, 7]}
    # plotchannels =['%s %s' % (x, binchannel) for x in binnames] + ["PfPR by Age Bin",
    #                                                                 "Annual Clinical Incidence by Age Bin"]
    plotchannels = ['RDT PfPR by Parasitemia and Age Bin']

    df = pd.read_csv(summary_fname)
    df.sort_values(by=['pop', xvar, yvar, 'month', 'Run_Number'])
    df['year'] = df['month'].apply(lambda x: int(x / 12) + 1)
    df = df[df['year'] == 2]

    for channel in plotchannels:
        for pop in ['Cohort', 'All']:

            sdf = df[df['pop'] == pop]
            if 'Clinical' in channel:
                sdf = sdf.groupby([xvar, yvar, 'year', 'Run_Number'])[channel].agg(np.mean).reset_index()
                sdf = sdf.groupby([xvar, yvar, 'year'])[channel].agg(np.sum).reset_index()
            else:
                sdf = sdf.groupby([xvar, yvar, 'year'])[channel].agg(np.mean).reset_index()

            rdf = sdf[sdf[xvar] == 0]
            rdf = rdf.rename(columns={channel: '%s_ref' % channel})

            sdf = pd.merge(left=sdf, right=rdf[[yvar, '%s_ref' % channel]],
                           on=[yvar], how='left')
            sdf['%s_reduction' % channel] = (sdf['%s_ref' % channel] - sdf[channel]) / sdf['%s_ref' % channel]
            sdf.loc[sdf['%s_reduction' % channel] < 0, '%s_reduction' % channel] = 0

            plot_heatmaps(sdf, channel, xvar, '%s_ref' % channel, '%s_reduction' % channel, 'summary %s' % pop,
                          ymin=ybounds[channel][0], ymax=ybounds[channel][1])


def plot_inset(inset_fname):

    xvar = 'killing'
    yvar = 'x_Temporary_Larval_Habitat'

    adf = pd.read_csv(inset_fname)
    inset_channels = [#'Blood Smear Parasite Prevalence', 'Infected', 'New Clinical Cases',
                      'PCR Parasite Prevalence', 'True Prevalence']

    ybounds = { 'Blood Smear Parasite Prevalence' : [0.1, 0.45],
                'Infected' : [0.1, 1],
                'PCR Parasite Prevalence': [0.1, 1],
                'True Prevalence': [0.1, 1],
                'New Clinical Cases' : [0, 20]}

    for channel in inset_channels:

        if 'Clinical' in channel:
            df = adf.groupby([xvar, yvar, 'Run_Number'])[channel].agg(np.mean).reset_index()
            df = df.groupby([xvar, yvar])[channel].agg(np.sum).reset_index()
        else:
            df = adf.groupby([xvar, yvar])[channel].agg(np.mean).reset_index()
        rdf = df[df[xvar] == 0]
        rdf = rdf.rename(columns={channel: '%s_ref' % channel})

        idf = pd.merge(left=df, right=rdf[[yvar, '%s_ref' % channel]], on=yvar)
        if channel == 'Blood Smear Parasite Prevalence':
            idf = idf[(idf['%s_ref' % channel] >= 0.1) & (idf['%s_ref' % channel] <= 0.47)]
        idf['%s_reduction' % channel] = (idf['%s_ref' % channel] - idf[channel]) / idf['%s_ref' % channel]
        idf.loc[idf['%s_reduction' % channel] < 0, '%s_reduction' % channel] = 0

        plot_heatmaps(idf, channel, xvar, '%s_ref' % channel, '%s_reduction' % channel, 'inset',
                      ymin=ybounds[channel][0], ymax=ybounds[channel][1])

if __name__ == '__main__' :

    # plot_summary(os.path.join(datadir, 'atsb_killing_hab_summary_50k.csv'))
    plot_inset(os.path.join(datadir, 'atsb_killing_hab_inset_50k.csv'))
