# Execute directly: 'python example_optimization.py'
# or via the calibtool.py script: 'calibtool run example_optimization.py'
import copy
import math
import random
import pandas as pd
import numpy as np
import os
from dtk.interventions.irs import add_IRS
from dtk.vector.species import update_species_param
from dtk.interventions.mosquito_release import add_mosquito_release
from dtk.interventions.itn_age_season import add_ITN_age_season

from malaria.reports.MalariaReport import add_habitat_report
from scipy.special import gammaln

from calibtool.CalibManager import CalibManager
from calibtool.algorithms.OptimTool import OptimTool
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from calibtool.plotters.OptimToolPlotter import OptimToolPlotter
from calibtool.plotters.SiteDataPlotter import SiteDataPlotter
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports.CustomReport import BaseVectorStatsReport
from malaria.reports.MalariaReport import add_filtered_report
from dtk.vector.species import set_larval_habitat
from simtools.SetupParser import SetupParser
from ATSB_Effect_Size_Calibration.spline_functions import get_spline_values
from dtk.interventions.novel_vector_control import add_ATSB
from ATSB_Effect_Size_Calibration.ATSBEntoCalibSite import ATSBEntoCalibSite


# Which simtools.ini block to use for this calibration
SetupParser.default_block = 'HPC'

# Start from a base MALARIA_SIM config builder
# This config builder will be modify by the different sites defined below
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.5)
update_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 0.9)
update_species_param(cb, 'gambiae', 'Vector_Sugar_Feeding_Frequency', 'VECTOR_SUGAR_FEEDING_EVERY_DAY', overwrite=False)
update_species_param(cb, 'gambiae', 'Anthropophily', 0.8)
update_species_param(cb, 'gambiae', 'Adult_Life_Expectancy', 20)

datadir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)')
# List of sites we want to calibrate on
throwaway = 1
species  = 'gambiae'
reference_fname = 'cluster_mosquito_counts_per_house_by_month.csv'
reference_spline = 'Multi_year_calibration_by_HFCA_180404/best_180409/Panjane_funestus.csv'
irs_fn = os.path.join(datadir,
                      'Malaria Team Folder/projects/Mozambique/entomology_calibration/cluster_all_irs_events_pyrethroid_2014.csv')
itn_fn = os.path.join(datadir,
                      'Malaria Team Folder/projects/Mozambique/entomology_calibration/cluster_all_itn_events.csv')

referenceData = os.path.join(datadir, 'Malaria Team Folder', 'projects', 'atsb', 'Reference_data', 'HLCIntervention.csv')
hfca = 'Kangaba'
sites = [ATSBEntoCalibSite(species, hfca=hfca, throwaway=throwaway,
                          reference_fname = referenceData)]

# WaningEffectMapPiecewise (Seasonal), WaningEffectConstant, WaningEffectExponential, WaningEffectBoxExponential
Class = "WaningEffectConstant"
expname = 'Mali ATSB Calibration ' + str(Class) + ' HLC'
max_iterations = 5
max_diff = 0.5
max_frac_diff = 0.05

# The default plotters used in an Optimization with OptimTool
# plotters = [LikelihoodPlotter(combine_sites=True),
plotters = [SiteDataPlotter(num_to_plot=20, combine_sites=True),
            OptimToolPlotter()  # OTP must be last because it calls gc.collect()
]


#############################################################################################################
ref = sites[0].analyzers[0].reference.copy()
duration = 1#int(max(ref.reset_index()['Month'])/12) + 1

params = [

      {
            #% of mosquitoes killed per day
           'Name': 'KillingProfile',
           'Dynamic': True,
           'Guess': 0.0337,
           'Min': 0.0,
           'Max': 0.5,
      },
      {
            # Only for exponential / box exponential
            'Name': 'HalfLife',
            'Dynamic': False,
            'Guess': 213,
            'Min': 1,
            'Max': 36500
      },
      {
            # Only for box / box exponential
           'Name': 'BoxDuration',
           'Dynamic': False,
           'Guess': 900,
           'Min': 1,
           'Max': 365
      },
      {
          # Map Piecewise data point (Currently set up with 2, but could be made with more
          'Name': 'WetSeason',
          'Dynamic': False,
          'Guess': 0.05,
          'Min': 0.0,
          'Max': 0.5
      },
      {
          # Map piecewise data point 2
          'Name': 'DrySeason',
          'Dynamic': False,
          'Guess': 0.11,
          'Min': 0.0,
          'Max': 1
      },
        # Map piecewise time variable
      {
          'Name': 'WetSeasonDuration',
          'Dynamic': False,
          'Guess': 100,
          'Min': 60,
          'Max': 120
      }
]

spline = get_spline_values(ref, species)

def map_sample_to_model_input(cb, sample):

    tags = {}
    # Updating values
    killRate = sample.pop('KillingProfile')
    tags.update({'KillingProfile': killRate})
    dry = sample.pop('DrySeason')
    tags.update({'DrySeason': dry})
    wet = sample.pop('WetSeason')
    tags.update({'WetSeason': wet})
    wetDur = sample.pop('WetSeasonDuration')
    tags.update({'WetSeasonDuration': wetDur})
    time = sample.pop('HalfLife')
    tags.update({'HalfLife': time})
    box = sample.pop('BoxDuration')
    tags.update({'BoxDuration': box})
    add_ATSB(cb, start = 517, coverage = 1.0, kill_cfg = {'Species': species,
                                'Killing_Config': {"class": Class,
                                                   # This intial effect is the killing rate, but if you're running a
                                                   # seasonal/piecewise model, the durability map values are multiplied
                                                   # by this, so set it to 1.0 instead of the killRate value.
                                                   "Initial_Effect": killRate,
                                                   #"Initial_Effect": 1.0
                                                    "Box_Duration": box,
                                                   "Decay_Time_Constant": time,
                                                   "Durability_Map": {
                                                        "Times": [0, wetDur],
                                                        "Values": [wet, dry]
                                                    }
                                }},
                  duration = 365)

    # For testing only, the duration should be handled by the site !! Please remove before running in prod!
    tags.update(cb.set_param("Simulation_Duration", (throwaway+duration)*365))
    tags.update(cb.set_param('Run_Number', 0))
    return tags

cb.update_params({'Demographics_Filenames': ['single_node_no_malaria_demographics.json'],

                      'Base_Population_Scale_Factor': 1,

                      "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
                      "Climate_Model": 'CLIMATE_BY_DATA',

                      "Air_Temperature_Filename": 'Burkina Faso_30arcsec_air_temperature_daily.bin',
                      "Land_Temperature_Filename": 'Burkina Faso_30arcsec_air_temperature_daily.bin',
                      "Rainfall_Filename": 'Burkina Faso_30arcsec_rainfall_daily.bin',
                      "Relative_Humidity_Filename": 'Burkina Faso_30arcsec_relative_humidity_daily.bin',

                      "Death_Rate_Dependence": "NONDISEASE_MORTALITY_BY_AGE_AND_GENDER",
                      "Disable_IP_Whitelist": 1,

                      "Enable_Vital_Dynamics": 1,
                      "Enable_Birth": 1,
                      'Enable_Default_Reporting': 1,
                      'Enable_Demographics_Other': 1,
                      # 'Enable_Property_Output': 1,

                      'x_Temporary_Larval_Habitat': 1,
                      'logLevel_default': 'ERROR',

                      "Vector_Species_Names": [species]
                      })

# Just for fun, let the numerical derivative baseline scale with the number of dimensions
volume_fraction = 0.01   # desired fraction of N-sphere area to unit cube area for numerical derivative (automatic radius scaling with N)
num_params = len([p for p in params if p['Dynamic']])

#Change radius
r = math.exp(1/float(num_params)*(math.log(volume_fraction) + gammaln(num_params/2.+1) - num_params/2.*math.log(math.pi)))
#r = 0.5
#r *= 0.5


optimtool = OptimTool(params,
    mu_r = r,           # <-- radius for numerical derivatve.  CAREFUL not to go too small with integer parameters
    sigma_r = r/10.,    # <-- stdev of radius
    center_repeats=1, # <-- Number of times to replicate the center (current guess).  Nice to compare intrinsic to extrinsic noise
    samples_per_iteration=2 # was 32  # <-- Samples per iteration, includes center repeats.  Actual number of sims run is this number times number of sites.
)


# cb.add_reports(BaseVectorStatsReport(type='ReportVectorStats', stratify_by_species=1))
add_filtered_report(cb, start=365*throwaway)

calib_manager = CalibManager(name=expname,    # <-- Please customize this name
                             config_builder=cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=1,  # <-- Replicates
                             max_iterations=max_iterations,   # <-- Iterations
                             plotters=plotters)


run_calib_args = {
    "calib_manager":calib_manager
}

if __name__ == "__main__":
    SetupParser.init()
    cm = run_calib_args["calib_manager"]
    cm.run_calibration()
