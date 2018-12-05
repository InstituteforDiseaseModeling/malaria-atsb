"""
Jaline Gerardin
Nov 2018


"""

import os
import pandas as pd
import numpy as np
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


projectdir = os.path.join('E:/', 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'atsb')


class PrevalenceAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, report_names=["AnnualAverage"], sweep_variables=None, working_dir="."):
        super(PrevalenceAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.sitenames=report_names
        self.expt_name = expt_name
        self.reference = pd.read_csv('site_details.csv')

    def select_simulation_data(self, data, simulation):
        simdata = []
        for site_name in self.sitenames:

            channeldata = data["output/MalariaSummaryReport_{name}.json".format(name=site_name)]["DataByTime"]["PfPR_2to10"]

            sitedata = pd.DataFrame({'PfPR2to10': channeldata,
                                    "Site_Name": site_name})
            sitedata = sitedata[:-1]
            simdata.append(sitedata)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected).reset_index(drop=True)
        df.to_csv(os.path.join(self.working_dir, '%s_PfPR.csv' % self.expt_name), index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    out_dir = os.path.join(projectdir, 'sim_data')

    sites = pd.read_csv("site_details.csv")

    experiments = {
                   "atsb_llin_v2" :"31c65386-86e7-e811-a2bd-c4346bcb1555"
                   }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id, analyzers=[PrevalenceAnalyzer(working_dir=out_dir,
                                                                     expt_name=expt_name,
                                                                     report_names = sites["name"].tolist(),
                                                                      sweep_variables=["Run_Number",
                                                                                       "x_Temporary_Larval_Habitat",
                                                                                       "intervention"
                                                                                       ])],
                            force_analyze=True)

        print(am.experiments)
        am.analyze()

