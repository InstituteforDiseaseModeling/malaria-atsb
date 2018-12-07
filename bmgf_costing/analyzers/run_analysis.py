import os
import sys
import pandas as pd
sys.path.append(os.path.dirname(__file__))

from simtools.Analysis.AnalyzeManager import AnalyzeManager
from prevalence_reduction_analyzer import PrevalenceAnalyzer
from atsb_llin_impact_analyzer import ATSBAnalyzer

if __name__ == "__main__":

    sites = pd.read_csv("site_details.csv")

    experiments = {
                   "atsb_llin_v3" :"8ae7cd0e-eaf8-e811-a2bd-c4346bcb1555"
                   }

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id, analyzers=[PrevalenceAnalyzer(expt_name=expt_name,
                                                                     report_names = sites["name"].tolist(),
                                                                      sweep_variables=["Run_Number",
                                                                                       "x_Temporary_Larval_Habitat",
                                                                                       "intervention"
                                                                                       ]),
                                                        ATSBAnalyzer(expt_name=expt_name,
                                                                     report_names=sites["name"].tolist(),
                                                                     sweep_variables=["Run_Number",
                                                                                      "x_Temporary_Larval_Habitat",
                                                                                      "intervention"
                                                                                      ])
                                                        ],
                            force_analyze=True)

    am.analyze()