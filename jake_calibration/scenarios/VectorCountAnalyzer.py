import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

import os

mpl.rcParams['pdf.fonttype'] = 42
if not SetupParser.initialized:
    SetupParser.init('HPC')

userpath = 'C:\\Users\\jkurlander'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
plotdir = os.path.join(wdir, 'sim_plots')
datadir = os.path.join(wdir, 'sim_data')

class VectorCountAnalyzer(BaseAnalyzer):
    def __init__(self):
        super(VectorCountAnalyzer, self).__init__()
        self.channel = 'Infected'
        self.filenames = ['output/InsetChart.json']

    def select_simulation_data(self, data, simulation):
        channeldata = data[self.filenames[0]]['Channels'][self.channel]['Data']
        simdata = pd.DataFrame({self.channel: channeldata,
                                'time': list(range(len(channeldata))),
                                'simLength': [[data[self.filenames[0]]['Header']['Timesteps']]]*len(channeldata)})
        for tag in ['coverage', 'duration', 'repetitions']:
           simdata[tag] = simulation.tags[tag]
        return simdata

    def finalize(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)

        simLength = df['simLength'][1][0]

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
                totalPlots +=1
                eliminated[str(c) + ' / ' + str(r)] = 0
                infect = [0.0]*1460
                infections[str(c) + ' / ' + str(r)] = 0
                while j+460 < len(sdf[self.channel]):
                    if sdf[self.channel][j] == sdf[self.channel][j+365]:
                        eliminated[str(c) + ' / ' + str(r)] += 1
                    #    infections[str(c) + ' / ' + str(r)] = sdf[self.channel]
                    for num in range (1460):
                        infect[num] += sdf[self.channel][j-1000+num]
                    j += 4*365
                for num in range(1460):
                    infect[num]
                ax.plot(range(1460), infect, linewidth = 1, label='cov %.2f' % (c))
        ax.legend()
        ax.set_xlabel('Time')
        ax.set_ylabel('Infected')
        #plt.savefig(os.path.join(plotdir, '123myplot123.pdf'), format='PDF')
        plt.title('ATSB Infection Reduction (CDC)')
        plt.show()
        plt.close('all')








if __name__ == '__main__' :

    expid = '905a91f4-2e9b-e811-a2c0-c4346bcb7275'
    am = AnalyzeManager(expid,
                        analyzers=VectorCountAnalyzer())
    am.analyze()
