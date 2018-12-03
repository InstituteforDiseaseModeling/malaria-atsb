import os
import json
import pandas as pd
import numpy as np

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser

from simtools.ModBuilder import ModBuilder, ModFn

from dtk.vector.study_sites import configure_site
from dtk.vector.species import update_species_param, set_larval_habitat
from dtk.interventions.novel_vector_control import add_ATSB
from dtk.interventions.outbreakindividual import recurring_outbreak
from dtk.interventions.property_change import change_individual_property
from dtk.interventions.biting_risk import change_biting_risk

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign

exp_name = 'ATSB_heatmap'
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

userpath = 'D:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
sim_setup_dir = os.path.join(wdir, 'sim_data')

species = 'gambiae'
years = 3
# Start of June in Mali
start_day = (years - 2) * 365 + 152
numseeds = 5

cb.update_params( {
    'Config_Name' : 'ATSB_Kangaba',
    'Simulation_Duration' : 365*years,
    'Vector_Species_Names' : [species],
    'Demographics_Filenames': ['demographics.json', 'demographics_cohort_IP_overlay.json'],

    'Base_Population_Scale_Factor': 1,

    "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
    "Climate_Model": 'CLIMATE_BY_DATA',

    "Enable_Default_Reporting" : 1,
    "Enable_Property_Output" : 0,
    "Disable_IP_Whitelist" : 1,
    "Disable_NP_Whitelist": 1,

    "Air_Temperature_Filename": 'Burkina Faso_30arcsec_air_temperature_daily.bin',
    "Land_Temperature_Filename": 'Burkina Faso_30arcsec_air_temperature_daily.bin',
    "Rainfall_Filename": 'Burkina Faso_30arcsec_rainfall_daily.bin',
    "Relative_Humidity_Filename": 'Burkina Faso_30arcsec_relative_humidity_daily.bin',
    # "Serialization_Time_Steps": [365*years],

    # 50-year burn-in with a 12-year climate cycle based on Burkina Faso
    #"Serialized_Population_Path": '\\\\internal.idm.ctr\\IDM\\Home\\jkurlander\\output\\50 year burn-in_20180801_171145\\b43\\765\\05a\\b4376505-ae95-e811-a2c0-c4346bcb7275\\output',
    "Serialized_Population_Filenames": ['state-18250.dtk'],
    "logLevel_VectorHabitat": 'ERROR',
})

def update_vector_params(cb) :

    update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.5)
    update_species_param(cb, 'gambiae', 'Vector_Sugar_Feeding_Frequency', 'VECTOR_SUGAR_FEEDING_EVERY_DAY', overwrite=False)
    update_species_param(cb, 'gambiae', 'Anthropophily', 0.8)
    update_species_param(cb, 'gambiae', 'Adult_Life_Expectancy', 20)

    # Habitat is based on a 1-year calibration of the control site of the Mali study.
    #CDC Habitats 0.001022239,  0.001720589,     0.001294523,   0.032008412,    0.044822194,    0.014206219,        0.617550763,    3.537562285,    3.907661784,    1.487261934,    0.579478307,    0.068954299     Gambiae, anthropophily .8, life expectancy 20, CDC Control larval habitat  7.08 Gambiae max
    #HLC Habitats 0.001,        0.001236706,    0.001933152,    0.056693638,    0.057953358,    0.015,              0.95,           2.159928736,    3.205076212,    0.43290933,     0.391090655,    0.138816133     Gambiae, anthropophily .8, life expectancy 20, HLC Control larval habitat, 7.44 Gambiae_max
    hab = {'Capacity_Distribution_Number_Of_Years': 1, 'Capacity_Distribution_Over_Time': {
        'Times': [0, 30.4166666666667, 60.8333333333333, 91.25, 121.6666666666667, 152.0833333333334, 182.5,
                  213.9166666666666, 243.3333333333334, 273.75, 304.5833333333333, 335],
        'Values': [0.001, 0.001236706, 0.001933152, 0.056693638, 0.057953358, 0.015, 0.95, 2.159928736, 3.205076212,
                   0.43290933, 0.391090655, 0.138816133]}, 'Max_Larval_Capacity': 2.75e7}
    set_larval_habitat(cb, {species: {'LINEAR_SPLINE': hab,
                                      'CONSTANT': 1e6 }})

def atsb_fn(cb, killing):

    add_ATSB(cb, start = start_day,
             coverage = 1.0, kill_cfg = {'Species': species,
                                    'Killing_Config': {"class": 'WaningEffectConstant',
                                                       "Initial_Effect": killing,
                                    }},
             duration = 10000)
    # add_ATSB(cb, coverage=coverage, start=100, duration=365, kill_cfg=killing_cfg[1])
    return {'killing': killing}


df = pd.read_csv(os.path.join(sim_setup_dir, 'LH_burnin.csv'))
# habs_to_sample = [float(x) for x in df['larval hab'].values]

update_vector_params(cb)
# change_biting_risk(cb, start_day=0,
#                    risk_config={'Risk_Distribution_Type': 'EXPONENTIAL_DURATION', 'Exponential_Mean': 1})
change_biting_risk(cb, start_day=0, trigger='Birth',
                   risk_config={'Risk_Distribution_Type': 'EXPONENTIAL_DURATION', 'Exponential_Mean': 1})
recurring_outbreak(cb, outbreak_fraction=0.01, start_day=152, repetitions=years, tsteps_btwn= 365)

change_individual_property(cb, 'CohortStatus', 'InCohort', start_day=start_day-2,
                           target={'agemin' : 0.5, 'agemax' : 17}, coverage=0.1, revert=365)
add_drug_campaign(cb, 'MDA', 'AL', start_days=[start_day+1],
                  ind_property_restrictions=[{"CohortStatus": "InCohort"}], repetitions=1)

add_summary_report(cb, start=start_day-365, interval=30, age_bins=[125],
                   parasitemia_bins=[0.1, 40, 1000000],
                   ipfilter='CohortStatus:InCohort', description='CohortMonthly')
add_summary_report(cb, start=start_day-365, interval=30, age_bins=[125],
                   parasitemia_bins=[0.1, 40, 1000000],
                   description='AllMonthly')
# add_event_counter_report(cb, ['Received_Campaign_Drugs'])



if __name__ == "__main__":

    SetupParser.init('HPC')

    builder = ModBuilder.from_list([[
        ModFn(atsb_fn, killing),
        ModFn(DTKConfigBuilder.update_params, {
            'x_Temporary_Larval_Habitat': float(row['x_Temporary_Larval_Habitat']),
            'Run_Number': s,
            'Serialized_Population_Path': os.path.join(row['outpath'], 'output')
        }),
    ] for r, row in df.iterrows() for s in range(numseeds) for killing in np.linspace(0, 0.25, 100)
    ])

    run_sim_args = {'config_builder': cb,
                    'exp_name': exp_name,
                    'exp_builder': builder
                    }

    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
