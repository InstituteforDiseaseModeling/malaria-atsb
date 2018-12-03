import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from VectorCountAnalyzer import VectorCountAnalyzer
import os

mpl.rcParams['pdf.fonttype'] = 42
if not SetupParser.initialized:
    SetupParser.init('HPC')

userpath = 'C:\\Users\\jkurlander'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
plotdir = os.path.join(wdir, 'sim_plots')
datadir = os.path.join(wdir, 'sim_data')


class EliminationAnalyzer(VectorCountAnalyzer):
    def __init__(self):
        super(EliminationAnalyzer, self).__init__()
        self.channel = 'Cumulative Infections'
        self.filenames = ['output/InsetChart.json']

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)

        sns.set_style('white', {'axes.linewidth': 0.5})
        fig = plt.figure()
        ax = fig.gca()
        eliminated = {}
        infections = {}
        totalPlots = 0
        coverages = df['coverage'].unique()
        repetitions = df['repetitions'].unique()
        for c in coverages:
            for r in repetitions:
                j = 1000
                sdf = df.loc[(df['coverage'] == c) & (df['repetitions'] == r)].reset_index()
                totalPlots += 1
                eliminated[str(c) + ' / ' + str(r)] = 0
                infect = [0.0] * 1460
                infections[str(c) + ' / ' + str(r)] = 0
                while j + 460 < len(sdf[self.channel]):
                    if sdf[self.channel][j] == sdf[self.channel][j + 365]:
                        eliminated[str(c) + ' / ' + str(r)] += 1
                    #    infections[str(c) + ' / ' + str(r)] = sdf[self.channel]
                    for num in range(1460):
                        infect[num] += sdf[self.channel][j - 1000 + num]
                    j += 4 * 365
                for num in range(1460):
                    infect[num] /= 100.0
                # ax.plot(range(1460), infect, linewidth = 1, label='cov %.2f' % (c))
        plotNames = [''] * totalPlots
        j = 0
        for (c, r), sdf in df.groupby(['coverage', 'repetitions']):
            plotNames[j] = str(c) + ' / ' + str(r)
            j += 1
        ax.legend()
        width = 0.6
        plotValues = [0] * totalPlots
        for p in range(totalPlots):
            plotValues[p] = eliminated[plotNames[p]] / 10.0
            # plotValues[p] = infections[plotNames[p]] / (100-eliminated[plotNames[p]])
        ax.bar(range(totalPlots), plotValues, width)
        plt.ylim((0, 100))

        # ax.set_ylabel('Infected ')
        ax.set_ylabel('% Chance of elimination')
        ax.set_xlabel('Coverage / Repetitions')
        ax.set_xticks(range(totalPlots))
        ax.set_xticklabels(plotNames)
        # plt.savefig(os.path.join(plotdir, '123myplot123.pdf'), format='PDF')
        plt.title('ATSB Vector Control Scenarios')
        # plt.title('ATSB Infection Reduction (CDC)')
        plt.show()
        plt.close('all')




if __name__ == '__main__':
    expid = '9c1e70d3-699a-e811-a2c0-c4346bcb7275'
    am = AnalyzeManager(expid,
                        analyzers=EliminationAnalyzer())
    am.analyze()
