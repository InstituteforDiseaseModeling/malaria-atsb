import logging
import os
import numpy as np
import pandas as pd
import calendar
from Helpers import multi_year_ento_data, multi_year_ento_data_clustered

from EntomologyCalibSite import EntomologyCalibSite
from ChannelByMultiYearSeasonCohortInsetAnalyzer import ChannelByMultiYearSeasonCohortInsetAnalyzer
from calibtool.study_sites.site_setup_functions import *

logger = logging.getLogger(__name__)


class MultiYearEntoCalibSite(EntomologyCalibSite):

    homedir = 'C:\\Users\\jkurlander\\Dropbox (IDM)'

    def __init__(self, spec, reference_fname,
                 itn_fname='', irs_fname='', hfca=None, **kwargs):
        self.metadata = {
        'village': 'Kangaba',
        'months': [calendar.month_abbr[i] for i in range(1, 13)],
        'species': [spec],
        'HFCA': hfca,
        'csvfilename' : reference_fname
         }
        if 'throwaway' in kwargs:
            self.throwaway = kwargs['throwaway']
        self.duration = int(max(self.get_reference_data('entomology_by_season').reset_index()['Month']) / 12) + 1
        self.irs_fname = irs_fname
        self.itn_fname = itn_fname
        super(MultiYearEntoCalibSite, self).__init__(hfca)

    def get_reference_data(self, reference_type):
        super(MultiYearEntoCalibSite, self).get_reference_data(reference_type)

        # Load the Parasitology CSV
        reference_data = multi_year_ento_data_clustered(self.metadata)

        return reference_data


    def get_setup_functions(self):

        def format_vc_df(fn, hfca):
            irs_events = pd.read_csv(fn)
            irs_events = irs_events[irs_events['cluster_name'] == hfca]
            irs_events['fulldate'] = pd.to_datetime(irs_events['fulldate'])
            irs_events['Year'] = irs_events['fulldate'].apply(lambda x: int(x.strftime('%y')))
            irs_events['Year'] = irs_events['Year'].apply(lambda x: x % min(irs_events['Year']))
            irs_events['simday'] = irs_events['fulldate'].apply(lambda x: int(x.strftime('%j')))
            irs_events['simday'] += irs_events['Year'] * 365

            return irs_events

        irs, itn = [], []
        if self.irs_fname :
            irs_df = format_vc_df(self.irs_fname, self.metadata['HFCA'])
            irs = [ add_irs_fn(start=int(row['simday']),
                        coverage=float(row['cov_all']),
                        waning={"Killing_Config": {
                "Initial_Effect": float(row['killing']),
                "Decay_Time_Constant": float(row['exp_duration']),
                "Box_Duration" :float(row['box_duration']),
                "class": "WaningEffectBoxExponential"
            }
            }) for r, row in irs_df.iterrows() ]
        if self.itn_fname :
            itn_df = format_vc_df(self.itn_fname, self.metadata['HFCA'])
            itn = [ add_itn_age_season_fn(start=int(row['simday']),
                                   age_dep={'youth_cov': float(row['age_cov']), 'youth_min_age': 5,
                                            'youth_max_age': 20},
                                   coverage=float(row['cov_all']),
                                   # as_birth=False, # no infrustructure for birth-triggered yet
                                   seasonal_dep={'min_cov': float(row['min_season_cov']), 'max_day': 1},
                                   discard={'halflife1': 260, 'halflife2': 2106,
                                            'fraction1': float(row['fast_fraction'])}) for r, row in itn_df.iterrows()]
        return irs + itn

    def get_analyzers(self):

        return [
            ChannelByMultiYearSeasonCohortInsetAnalyzer(site=self, duration=self.duration)]

