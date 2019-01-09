import os
import sys
sys.path.append(os.path.dirname(__file__))

from simtools.Analysis.AnalyzeManager import AnalyzeManager
from HeatmapAnalyzer import InsetAnalyzer, SummaryReportAnalyzer

if __name__ == '__main__' :

    expid = '8da601be-cda4-e811-a2c0-c4346bcb7275'
    am = AnalyzeManager(expid,
                        analyzers=[InsetAnalyzer(expname='atsb_killing_hab_inset_50k',
                                                 channels=['PCR Parasite Prevalence', 'True Prevalence',
                                                           'Blood Smear Parasite Prevalence', 'Infected',
                                                           'New Clinical Cases']),
                                   # SummaryReportAnalyzer(expname='atsb_killing_hab_summary_50k', force=False)
                                   ]
                        )
    am.analyze()
