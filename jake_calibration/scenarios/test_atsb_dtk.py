import os
import json
import pandas as pd
import numpy as np

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser

from simtools.ModBuilder import ModBuilder, ModFn

from dtk.vector.study_sites import configure_site
from dtk.interventions.novel_vector_control import add_ATSB
from dtk.vector.species import update_species_param


exp_name = 'atsb_campaign_class_test'
cb = DTKConfigBuilder.from_defaults('VECTOR_SIM')
# cb.campaign.Use_Defaults = False
configure_site(cb, 'Namawala')

species = 'funestus'

cb.update_params( {
    'Config_Name' : 'ATSB_Test_Namawala',
    'Simulation_Duration' : 365,
    'Vector_Species_Names' : [species, 'gambiae']
})

update_species_param(cb, 'funestus', 'Vector_Sugar_Feeding_Frequency', 'VECTOR_SUGAR_FEEDING_EVERY_DAY', overwrite=False)
update_species_param(cb, 'gambiae', 'Vector_Sugar_Feeding_Frequency', 'VECTOR_SUGAR_FEEDING_EVERY_DAY', overwrite=False)

def atsb_fn(cb, coverage, killing) :

    killing_cfg = [{
            "Species" : 'funestus',
            "Killing_Config" : {
                "class": "WaningEffectMapLinearSeasonal",
                "Initial_Effect" : 1.0,
                "Durability_Map" :
                {
                    "Times"  : [ 0.0, 90.0, 91.0, 300.0, 301.0, 365.0 ],
                    "Values" : [ 1.0, 1.0,   0.01,  0.01,  1.0,   1.0 ]
                }
            }
    },
        {
        "Species" : 'gambiae',
            "Killing_Config" : {
                "class": "WaningEffectBoxExponential",
                "Initial_Effect": killing,
                "Box_Duration": 180,
                "Decay_Time_Constant": 30
            }
        }
        ]

    add_ATSB(cb, coverage=coverage, start=100, duration=365,
             kill_cfg=killing_cfg
             )
    return { 'coverage' : coverage,
             'initial_killing' : killing}


# add_ATSB(cb, start=100, coverage=1, duration=365)

builder = ModBuilder.from_list( [ [
    ModFn(atsb_fn, cov, kill),
    ModFn(DTKConfigBuilder.set_param, 'Run_Number', x)
    ] for cov in [0.1, 1] for kill in [0.15] for x in range(1)
])


run_sim_args = {'config_builder': cb,
                'exp_name': exp_name,
                'exp_builder': builder
                }


if __name__ == "__main__":
    SetupParser.init('HPC')
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    # Wait for the simulations to be done
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
