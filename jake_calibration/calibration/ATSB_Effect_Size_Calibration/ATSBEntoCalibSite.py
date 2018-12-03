import logging
import os
import numpy as np
import pandas as pd
import calendar
from ATSB_Ento_Calibration.Helpers import multi_year_ento_data, multi_year_ento_data_clustered
from dtk.vector.species import update_species_param
from dtk.interventions.novel_vector_control import add_ATSB
from ATSB_Ento_Calibration.EntomologyCalibSite import EntomologyCalibSite
from ATSB_Ento_Calibration.ChannelByMultiYearSeasonCohortInsetAnalyzer import ChannelByMultiYearSeasonCohortInsetAnalyzer
from malaria.study_sites.site_setup_functions import *

logger = logging.getLogger(__name__)


class ATSBEntoCalibSite(EntomologyCalibSite):

    homedir = 'C:\\Users\\jkurlander\\Dropbox (IDM)'

    def __init__(self, spec, reference_fname,
                 itn_fname='', irs_fname='', atsb_fname='', hfca=None, **kwargs):
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
        self.atsb_fname = atsb_fname

        super(ATSBEntoCalibSite, self).__init__(hfca)

    def get_reference_data(self, reference_type):
        super(ATSBEntoCalibSite, self).get_reference_data(reference_type)

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

    # To change the population numbers, leave the first and last values in "values".
        # Replace the middle 12 with the 12-month habitat numbers.

    # Jake's 2018 Mali CDC Habitats 0.001022239,   0.001720589, 0.001294523, 0.032008412, 0.044822194, 0.014206219,
        # 0.617550763, 3.537562285, 3.907661784, 1.487261934, 0.579478307, 0.068954299
        # Funestus, anthropophily .8, life expectancy 20, CDC Control larval habitat 7.08  max

    # Jake's 2018 Mali HLC Habitats 0.001, 0.001236706, 0.001933152, 0.056693638, 0.057953358, 0.015, 0.95, 2.159928736,
        # 3.205076212, 0.43290933, 0.391090655, 0.138816133
        # Funestus, anthropophily .8, life expectancy 20, HLC Control larval habitat, 7.44 max

        hab = {'Capacity_Distribution_Number_Of_Years': 2, 'Capacity_Distribution_Over_Time': {
            'Times': [0, 365.0, 395.4166666666667, 425.8333333333333, 456.25, 486.6666666666667, 517.0833333333334,
                      547.5, 577.9166666666666, 608.3333333333334, 638.75, 669.1666666666667, 699.5833333333333, 729],
            'Values': [0.01, 0.001, 0.001236706, 0.001933152, 0.056693638, 0.057953358, 0.015, 0.95, 2.159928736,
                       3.205076212, 0.43290933, 0.391090655, 0.138816133, 0.01]}, 'Max_Larval_Capacity': 2.75e7}

        setupFunctions = [species_param_fn('funestus', 'Indoor_Feeding_Fraction', 0.9),
                        larval_habitat_fn('gambiae', {'LINEAR_SPLINE' : hab})]

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
        # if self.atsb_fname :
        #     atsb_df = format_vc_df(self.atsb_fname, self.metadata['HFCA'])
        #     atsb = [add_ATSB(
        #         i dont think this is right so commented out
        #         ) for x in atsb_df.iterrows()]
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
        return setupFunctions + irs + itn

    def get_analyzers(self):

        return [
            ChannelByMultiYearSeasonCohortInsetAnalyzer(site=self, duration=self.duration)]

