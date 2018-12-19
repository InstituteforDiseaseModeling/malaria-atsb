import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
import os

if not SetupParser.initialized:
    SetupParser.init('HPC')

userpath = 'E:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
plotdir = os.path.join(wdir, 'sim_plots')
datadir = os.path.join(wdir, 'sim_data')


class InsetAnalyzer(BaseAnalyzer):
    def __init__(self, expname, channels=None):
        super(InsetAnalyzer, self).__init__()
        self.filenames = ['output/InsetChart.json']
        self.channels = channels if channels else ['PCR Parasite Prevalence', 'True Prevalence',
                                                   'Blood Smear Parasite Prevalence', 'Infected',
                                                   'New Clinical Cases']
        self.xvar = 'killing'
        self.yvar = 'x_Temporary_Larval_Habitat'
        self.expname = expname

    def select_simulation_data(self, data, simulation):
        simdata = pd.DataFrame({ x : data[self.filenames[0]]['Channels'][x]['Data'][-578:-213] for x in self.channels })
        simdata['time'] = simdata.index

        for tag in simulation.tags:
            if tag in ['killing', 'x_Temporary_Larval_Habitat', 'Run_Number'] :
               simdata[tag] = simulation.tags[tag]
        return simdata

    def finalize(self, all_data):
        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)
        # df = df.groupby(['killing', 'x_Temporary_Larval_Habitat', 'Run_Number']).agg(np.mean).reset_index(drop=True)
        df.to_csv(os.path.join(self.working_dir, '%s.csv' % self.expname), index=False)


class SummaryReportAnalyzer(BaseAnalyzer):

    def __init__(self, expname,
                 channels=["PfPR by Age Bin",
                           "Annual Clinical Incidence by Age Bin",
                           "PfPR by Parasitemia and Age Bin"],
                 force=False, inset_csv_fname='atsb_killing_hab_inset'):
        super(SummaryReportAnalyzer, self).__init__()
        self.channels = channels
        self.descriptions = ['All', 'Cohort']
        self.filenames = ['output/MalariaSummaryReport_%sMonthly.json' % x for x in self.descriptions]
        self.xvar = 'killing'
        self.yvar = 'x_Temporary_Larval_Habitat'
        self.detection_channels = ['true', 'PCR', 'RDT']
        self.expname = expname
        self.force = force
        self.inset_csv_fname = os.path.join(datadir, '%s.csv' % inset_csv_fname)
        if not self.force :
            try :
                self.df = pd.read_csv(os.path.join(datadir, '%s.csv' % self.expname))
            except FileNotFoundError :
                self.force = True

    def select_simulation_data(self, data, simulation):

        if self.force :
            df = pd.DataFrame()
            for pop, fname in zip(self.descriptions, self.filenames) :

                d = data[fname]['DataByTimeAndAgeBins']
                time_age_data = { key : val for key, val in d.items() if key in self.channels }
                d = data[fname]['DataByTimeAndPfPRBinsAndAgeBins']
                time_age_par_data = { key : val for key, val in d.items() if key in self.channels}

                for key, val in time_age_data.items() :
                    time_age_data[key] = [x[0] for x in val]

                simdata = pd.DataFrame(time_age_data)

                simdata['month'] = simdata.index - 1
                simdata['pop'] = pop

                t = {}

                for key, val in time_age_par_data.items() :
                    d1 = [[x[i] for x in val] for i in range(3)]
                    d2 = [[x[0] for x in d] for d in d1]
                    for i, name in enumerate(self.detection_channels) :
                        t['%s %s' % (name, key)] = d2[i]

                simdata1 = pd.DataFrame(t)
                simdata1['PCR %s' % key] = simdata1['PCR %s' % key] + simdata1['RDT %s' % key]
                simdata1['true %s' % key] = simdata1['true %s' % key] + simdata1['PCR %s' % key]
                simdata1['month'] = simdata.index - 1

                simdata = pd.merge(left=simdata, right=simdata1, on='month')
                df = pd.concat([df, simdata])

            for tag in simulation.tags:
                if tag in ['killing', 'x_Temporary_Larval_Habitat', 'Run_Number'] :
                   df[tag] = simulation.tags[tag]
            return df
        else :
            return pd.DataFrame()

    def finalize(self, all_data):
        if self.force :

            selected = [data for sim, data in all_data.items()]
            if len(selected) == 0:
                print("No data have been returned... Exiting...")
                return

            df = pd.concat(selected).reset_index(drop=True)
            df.to_csv(os.path.join(datadir, '%s.csv' % self.expname), index=False)
        else :
            df = self.df

        refchannel = 'Blood Smear Parasite Prevalence'

        idf = pd.read_csv(self.inset_csv_fname)
        idf = idf.groupby([self.xvar, self.yvar, 'Run_Number'])[refchannel].agg(np.mean).reset_index()
        idf = idf.groupby([self.xvar, self.yvar])[refchannel].agg(np.mean).reset_index()
        idf = idf[idf[self.xvar] == 0]
        idf = idf.rename(columns={refchannel: '%s_ref' % refchannel})

        pop = self.descriptions[0]
        channel = self.channels[0]

        sdf = df[df['pop'] == pop]
        sdf = sdf.groupby([self.xvar, self.yvar, 'month'])[channel].agg(np.mean).reset_index()

        sdf = pd.merge(left=sdf, right=idf, on=[self.xvar, self.yvar], how='left')
        pass


if __name__ == '__main__' :

    expid = '8da601be-cda4-e811-a2c0-c4346bcb7275'
    am = AnalyzeManager(expid,
                        analyzers=[InsetAnalyzer(expname='atsb_killing_hab_inset_50k',
                                                 channels=['PCR Parasite Prevalence', 'True Prevalence',
                                                           'Blood Smear Parasite Prevalence', 'Infected',
                                                           'New Clinical Cases'], force=False),
                                   # SummaryReportAnalyzer(expname='atsb_killing_hab_summary_50k', force=False)
                                   ]
                        )
    am.analyze()
