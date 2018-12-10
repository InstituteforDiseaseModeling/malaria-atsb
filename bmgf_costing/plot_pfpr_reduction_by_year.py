"""
Jaline Gerardin
Dec 2018

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
    df = df.groupby(['Site_Name', 'x_Temporary_Larval_Habitat', 'intervention', 'year']).agg(np.mean).reset_index()

    df = df.sort_values(by=['Site_Name', 'x_Temporary_Larval_Habitat', 'year'])

    return df


def plot_pfpr(df) :

    baseline_channel = 'itn'
    ignore_list = ['none', 'itn']
    datachannel = 'PfPR2to10'

    num_interventions = len(df['intervention'].unique()) - len(ignore_list)
    num_years = len(df['year'].unique())

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    fig = plt.figure('PfPR', figsize=(15,10))
    axes = [fig.add_subplot(num_interventions, num_years, d + 1) for d in range(num_years*num_interventions)]
    palette = load_color_palette()

    for s, (site, sdf) in enumerate(df.groupby('Site_Name')):
        sdf = sdf.sort_values(by='PfPR2to10')
        for d, (year, ydf) in enumerate(sdf.groupby('year')):
            baseline = sdf[(sdf['intervention'] == baseline_channel) & (sdf['year'] == year)]
            for i, (intervention, idf) in enumerate(ydf[~(ydf['intervention'].isin(ignore_list))].groupby('intervention')) :
                ax = axes[i*num_years+d]
                xvar = [y for x, y in zip(idf[datachannel].values, baseline[datachannel].values) if (y > 0)]
                yvar = [x for x, y in zip(idf[datachannel].values, baseline[datachannel].values) if (y > 0)]
                ys = lowess(yvar, xvar, frac=0.2)[:,1]
                ax.plot(xvar, ys, '-', color=palette[s], label=site)
                ax.set_title('%s year %d' % (intervention, year))
                if i == num_interventions-1 :
                    ax.set_xlabel('%s no intervention' % datachannel)
                if d%num_years == 0 :
                    ax.set_ylabel('%s with intervention' % datachannel)

    axes[-1].legend()
    fig.savefig(os.path.join(plotdir, '%s.png' % savename))
    fig.savefig(os.path.join(plotdir, '%s.pdf' % savename), format='PDF')
    plt.close(fig)


if __name__ == '__main__' :

    expt_name = "atsb_llin_v3"
    data_fname = os.path.join(datadir, "%s_pfpr.csv" % expt_name)

    df = load_sim_df(data_fname)

    savename = '%s_pfpr_by_site' % expt_name
    plot_pfpr(df)

    plt.show()
