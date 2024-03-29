import pandas as pd
import json
import math

import os
import pdb

from dtk.vector.species import set_params_by_species, set_species_param
from dtk.interventions.habitat_scale import scale_larval_habitats
from malaria.interventions.health_seeking import add_health_seeking
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property
from dtk.interventions.novel_vector_control import add_ATSB

from malaria.interventions.malaria_drug_campaigns import add_drug_campaign

def assign_net_ip(cb, hates_net_prop):
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop),
    change_individual_property(cb, "NetUsage", "HatesNets", coverage=hates_net_prop,
          trigger_condition_list=["Births"])
    return {"Hates_Nets": hates_net_prop}

def simulation_setup(cb, species_details, site_vector_props, max_larval_capacity=4e8):

    site_dir = os.path.join("sites", "all")

    # directories
    cb.update_params({
                    "Demographics_Filenames": [os.path.join(site_dir, "demographics.json"),
                                               os.path.join(site_dir, "demographics_net_overlay.json")],
                    "Air_Temperature_Filename": os.path.join(site_dir,
                                                             "air_temperature_daily.bin"),
                    "Land_Temperature_Filename": os.path.join(site_dir,
                                                              "air_temperature_daily.bin"),
                    "Rainfall_Filename": os.path.join(site_dir,
                                                      "rainfall_daily.bin"),
                    "Relative_Humidity_Filename": os.path.join(site_dir,
                                                               "relative_humidity_daily.bin")
                    }
    )

    # Find vector proportions for each vector
    set_params_by_species(cb.params, [name for name in species_details.keys()])

    larval_habs_per_site = {"NodeID": site_vector_props["node_id"]}

    for species_name, species_modifications in species_details.items():
        set_species_param(cb, species_name, "Adult_Life_Expectancy", 20)
        set_species_param(cb, species_name, 'Vector_Sugar_Feeding_Frequency', 'VECTOR_SUGAR_FEEDING_EVERY_DAY')

        for param, val in species_modifications.items():
            if param == "habitat_split":
                new_vals = {hab: hab_prop * max_larval_capacity for hab, hab_prop in val.items()}
                set_species_param(cb, species_name, "Larval_Habitat_Types", new_vals)
                larval_habs_per_site.update({".".join([species_name, hab]): site_vector_props[species_name]
                                             for hab in val.keys()})
            else:
                set_species_param(cb, species_name, param, val)

    scale_larval_habitats(cb, pd.DataFrame(larval_habs_per_site))


# itns
def add_annual_itns(cb, year_count=1, coverage=0.8, start_day=0, initial_killing=0.3, discard_time=270, IP=[]):

    for year in range(year_count):
        add_ITN_age_season(cb,
                           coverage_all=coverage,
                           discard={"halflife": discard_time},
                           waning={'kill_initial': initial_killing},
                           start=(365 * year) + start_day,
                           ind_property_restrictions=IP)

    return {"ITN_Coverage": coverage,
            "ITN_Start": start_day,
            "ITN_killing" : initial_killing,
            'ITN_discard' : discard_time}


# irs
def add_irs_group(cb, coverage=1.0, start_days=[0], decay=270):

    waning = {
        "Killing_Config": {
            "Initial_Effect": 0.6,
            "Decay_Time_Constant": decay,
            "class": "WaningEffectExponential"
        },
        "Blocking_Config": {
            "Initial_Effect": 0.0,
            "Decay_Time_Constant": 730,
            "class": "WaningEffectExponential"
        }}

    for start in start_days:
        add_IRS(cb, start, [{"min": 0, "max": 200, "coverage": coverage}],
                waning=waning)

    return {"IRS_Start": start_days[0], "IRS_Coverage": coverage}


# atsb
def add_atsb_by_coverage(cb, coverage=1, killing = 0.0337, species_list=[]):

    add_ATSB(cb, start = 5,
             coverage = coverage,
             kill_cfg = [{ 'Species' : sp,
                           'Killing_Config' : {"class": 'WaningEffectConstant',
                                               "Initial_Effect": killing }
                           } for sp in species_list],
             duration=3*365)
    return {'atsb_coverage': coverage}


# act
def add_healthseeking_by_coverage(cb, coverage=1.0, rate=0.15, drugname="AL"):
    drugs = {"AL": ["Artemether", "Lumefantrine"],
             "DP": ["DHA", "Piperaquine"]}

    add_health_seeking(cb,
                       targets=[{"trigger": "NewClinicalCase",
                                 "coverage": coverage,
                                 "agemin": 0,
                                 "agemax": 100,
                                 "seek": 1.0,
                                 "rate": rate}],
                       drug=drugs[drugname],
                       dosing="FullTreatmentNewDetectionTech",
                       nodes={"class": "NodeSetAll"},
                       repetitions=1,
                       tsteps_btwn_repetitions=365,
                       broadcast_event_name="Received_Treatment")

    return {"CM_Coverage": coverage, "CM_Daily_Prob": rate, "CM_Drug": drugname}

