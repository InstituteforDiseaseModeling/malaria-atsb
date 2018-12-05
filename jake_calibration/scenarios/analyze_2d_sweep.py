import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

import os
from plotting.colors import load_color_palette

mpl.rcParams['pdf.fonttype'] = 42
if not SetupParser.initialized:
    SetupParser.init('HPC')

userpath = 'D:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
plotdir = os.path.join(wdir, 'sim_plots')
datadir = os.path.join(wdir, 'sim_data')


class Sweep2DAnalyzer(BaseAnalyzer) :

    def __init__(self, sweep_vars):
        super(Sweep2DAnalyzer, self).__init__()
        self.sweep_variables = sweep_vars
        self.channel = 'Adult Vectors'
        self.filenames = ['output/InsetChart.json']

    def select_simulation_data(self, data, simulation):
        channeldata = data[self.filenames[0]]['Channels'][self.channel]['Data']
        simdata = pd.DataFrame( { self.channel : channeldata,
                                  'time' : list(range(len(channeldata)))})
        for tag in self.sweep_variables + ['Run_Number']:
            simdata[tag] = simulation.tags[tag]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim,data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        d = pd.concat(selected).reset_index(drop=True)
        d['year'] = d['time'].apply(lambda x : int(x/365))

        d = d[d['time'] >= 100]
        summed = d.groupby(self.sweep_variables + ['Run_Number'])[self.channel].agg(np.sum).reset_index()
        df = summed.groupby(self.sweep_variables)[self.channel].agg(np.mean).reset_index()

        print(df)

        sns.set_style('white', {'axes.linewidth' : 0.5})
        fig = plt.figure()
        ax = fig.gca()
        sns.heatmap(df.pivot(index=self.sweep_variables[0], columns=self.sweep_variables[1], values=self.channel),
                    ax=ax)
        ax.set_ylabel(self.sweep_variables[0])
        ax.set_xlabel(self.sweep_variables[1])
        plt.show()
        plt.close('all')


if __name__ == '__main__' :

    expid = '63903283-1f7b-e811-a2c0-c4346bcb7275'
    am = AnalyzeManager(expid,
                        analyzers=Sweep2DAnalyzer(['coverage', 'initial_killing']))
    am.analyze()
