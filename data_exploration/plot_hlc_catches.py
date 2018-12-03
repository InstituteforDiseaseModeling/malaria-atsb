import os
import json
import pandas as pd
import numpy as np
import datetime

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from plotting.colors import load_color_palette
mpl.rcParams['pdf.fonttype'] = 42


rootdir = 'C:/Users/jkurlander/Dropbox (IDM)'
datadir = os.path.join(rootdir, 'Malaria Team Folder/data/Novel VC/ATSB/phase 2 data/merged_2018_07')
plotdir = os.path.join(rootdir, 'Malaria Team Folder/projects/atsb/data_exploration_plots')

hlc_fname = os.path.join(datadir, 'HLC_raw_all_dates.csv')

atsb_start = datetime.date(2017, 6, 1)
intervention_villages = ['Krekrelo', 'Sirakele', 'Trekourou', 'Farabale', 'Kignele', 'Tiko', 'Sambadani']


def add_date_column() :

    df = pd.read_csv(hlc_fname)

    import calendar

    cal = {v: k for k,v in enumerate(calendar.month_abbr)}
    df['date'] = df.apply(lambda x : datetime.date(x['Year'], cal[x['Month'][:3]], 1), axis=1)
    df.to_csv(hlc_fname, index=False)

def load_df() :

    df = pd.read_csv(hlc_fname)
    df['date'] = pd.to_datetime(df['date'])
    df['arm'] = 'control'
    df.loc[df['Study Sites'].isin(intervention_villages), 'arm'] = 'intervention'
    return df

def aggregate_counts(df, sep_by_position=True) :

    grouplist = ['arm', 'Study Sites', 'date', 'Position']
    if not sep_by_position :
        grouplist.pop()

    grouped = df.groupby(grouplist)
    sdf = grouped['mosquito ID'].agg(len).reset_index()
    sdf.rename(columns={'mosquito ID': 'vector count'}, inplace=True)
    df = sdf.sort_values(by=['Study Sites', 'date'])

    return df

def plot_hlc_by_site() :

    df = load_df()

    palette = load_color_palette()

    fig = plt.figure('HLC vector counts', figsize=(16,10))
    axes = [fig.add_subplot(4, 4, x + 1) for x in range(14)]
    fig.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.05)
    linestyles = ['-o', '--o']
    for locindex, (location, ldf) in enumerate(df.groupby('Position')) :
        grouped = ldf.groupby(['Study Sites', 'date'])
        sdf = grouped['mosquito ID'].agg(len).reset_index()
        sdf.rename(columns={'mosquito ID' : 'vector count'}, inplace=True)
        sdf.sort_values(by=['Study Sites', 'date'], inplace=True)

        for s, (site, ssdf) in enumerate(sdf.groupby('Study Sites')) :
            ax = axes[s]
            ax.plot(ssdf['date'], ssdf['vector count'], linestyles[locindex], color=palette[s],
                    label=location if s == 0 else '')
            ax.set_title(site)
            ax.set_xlim(np.min(df['date']), np.max(df['date']))
            if s < 10 :
                ax.set_xticklabels([])
            else :
                ax.set_xticklabels(['Jun \'16', '', '', 'Dec \'16', '', '', 'Jun \'17', '', '', 'Dec \'17', ''])
            if s == 0 and locindex > 0 :
                ax.legend()
            if s%4 == 0 :
                ax.set_ylabel('mosquito count')
            if locindex == 1 and site in intervention_villages :
                ymin, ymax = ax.get_ylim()
                ax.fill_between([atsb_start, np.max(df['date'])], [ymin,ymin], [ymax, ymax],
                                color='k', alpha=0.25, linewidth=0)
                ax.set_ylim(ymin, ymax)
    plt.savefig(os.path.join(plotdir, 'HLC_by_site.png'))
    plt.savefig(os.path.join(plotdir, 'HLC_by_site.pdf'), format='PDF')


def calculate_change_by_month(df) :

    change_by_sites = {'month' : list(range(6,13))}
    for site, sdf in df.groupby('Study Sites'):
        changes = []
        for month in range(6,13) :
            old = sdf[sdf['date'] == datetime.date(2016, month, 1)]
            new = sdf[sdf['date'] == datetime.date(2017, month, 1)]
            if len(old) == 0 or len(new) == 0 :
                changes.append(-1)
                continue
            fold_change = (np.sum(new['vector count']) - np.sum(old['vector count']))/np.sum(old['vector count'])
            changes.append(fold_change)
        change_by_sites[site] = changes
    return pd.DataFrame(change_by_sites)

def plot_year_on_year_fold_change_by_site(df) :

    df = aggregate_counts(df)
    fold_changes = calculate_change_by_month(df)

    palette = load_color_palette()
    fig = plt.figure()
    ax = fig.gca()
    for i, village in enumerate(df['Study Sites'].unique()) :
        if village in intervention_villages :
            linestyle = '-o'
            arm = 'int'
        else :
            linestyle = '--o'
            arm = 'con'
        sdf = fold_changes[fold_changes[village] != -1]
        ax.plot(sdf['month'], sdf[village]*100, linestyle, color=palette[i], label='%s %s' % (village, arm))

    ax.set_ylabel('year on year pct change in vector count')
    ax.set_xticklabels(['', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    ax.legend()

    plt.savefig(os.path.join(plotdir, 'HLC_year_on_year_pct_change_by_site.png'))
    plt.show()

def plot_int_vs_control() :

    df = aggregate_counts(load_df())
    meandf = df.groupby(['arm', 'date', 'Position'])['vector count'].agg(np.mean).reset_index()

    sns.set_style('whitegrid')
    fig = plt.figure()
    ax = plt.gca()
    for pos, pdf in meandf.groupby('Position') :
        data = []
        for month in range(6, 13):
            date = datetime.date(2017, month, 1)
            control = pdf[(pdf['arm'] == 'control') & (pdf['date'] == date)]
            intervention = pdf[(pdf['arm'] == 'intervention') & (pdf['date'] == date)]
            c = control['vector count'].values[0]
            i = intervention['vector count'].values[0]
            data.append((i - c)/c*100)
        ax.plot(list(range(6,13)), data, '-o', label=pos)

    ax.set_xticklabels(['', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    ax.set_ylabel('2017 percent diff in vector count in int vs control')
    ax.set_ylim(0,-100)
    ax.legend()
    plt.savefig(os.path.join(plotdir, 'HLC_pct_change_int_vs_control.png'))
    plt.show()

if __name__ == '__main__' :

    # plot_hlc_by_site()
    plot_int_vs_control()