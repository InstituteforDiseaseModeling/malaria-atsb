import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from datetime import date

dropbox_path = 'E:/Dropbox (IDM)/Malaria Team Folder'
datadir = os.path.join(dropbox_path, 'data/Novel VC/ATSB/phase 2 data')
plotdir = os.path.join(dropbox_path, 'projects/atsb/data_exploration_plots')


def plot_hlc(ax, color) :

    hlc_data_fname = os.path.join(datadir, 'merged_2018_07/HLC_raw_all_dates.csv')
    col = 'greater than 3 GC'

    adf = pd.read_csv(hlc_data_fname)

    adf['date'] = pd.to_datetime(adf['date'])
    adf.sort_values(by='date')

    df = adf[adf[col].isin(['0', '1'])]
    df[col] = df[col].astype(int)

    def frac_one(x) :
        return np.sum(x)/len(x)

    gdf = df.groupby('date')[col].agg(frac_one).reset_index()
    ax.plot(gdf['date'], gdf[col], '-o', color=color, label='HLC')


def plot_traps(ax, sheet_name, color) :

    data_fname = os.path.join(datadir, 'data/Data year 2/ATSB trial 2017/DB CDC Malaise PSC catches 2017.xls')
    year = 2017
    col_gr, col_ls = '>3 gon. Cyc.', '< 3 gon. Cyc.'
    col = 'fraction >3 GC'

    adf = pd.read_excel(data_fname, sheet_name=sheet_name)
    df = adf.groupby('Month')[[col_gr, col_ls]].agg(np.sum).reset_index()

    df[col] = df[col_gr]/(df[col_gr] + df[col_ls])
    df['date'] = df['Month'].apply(lambda x : date(year, int(x), 1))

    ax.plot(df['date'], df[col], '-o', color=color, label=sheet_name)


if __name__ == '__main__' :

    fig = plt.figure(figsize=(9,4))
    ax = fig.gca()
    palette = sns.color_palette('Set2')
    plot_hlc(ax, palette[0])

    for i, trap in enumerate(['CDC', 'Malaise', 'PSC']) :
        plot_traps(ax, trap, palette[i+1])

    ax.set_ylabel('fraction >3 GC')
    ax.set_ylim(0,)
    ax.legend()

    fig.savefig(os.path.join(plotdir, 'fraction_greater_than_3_GC_2017.png'))

    plt.show()
