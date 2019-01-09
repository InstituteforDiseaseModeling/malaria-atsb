import os
import pandas as pd
import numpy as np
import json

from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import get_asset_collection

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.Utilities.COMPSUtilities import COMPS_login

from malaria.reports.MalariaReport import add_summary_report, add_event_counter_report
from simtools.Utilities.Experiments import retrieve_experiment
from malaria.interventions.health_seeking import add_health_seeking

from sweep_functions import *

# variables
run_type = "intervention"  # set to "burnin" or "intervention"
burnin_id = "96e9c858-a8ce-e811-a2bd-c4346bcb1555"
asset_exp_id = "66d8416c-9fce-e811-a2bd-c4346bcb1555"

intervention_coverages = [100]
interventions = ["llin_no_disc"]
num_runs = 40
# hs_daily_probs = [0.15, 0.3, 0.7]

hates_net_prop = 0.1 # based on expert opinion from Caitlin
new_inputs = False

# Serialization
print("setting up")
if run_type == "burnin":
    years = 15
    sweep_name = "MAP_II_New_Sites_Burnin"
    serialize = True
    pull_from_serialization = False
elif run_type == "intervention":
    years = 3
    sweep_name = "ATSB_LLIN_multisite_withHS_v2"
    serialize = False
    pull_from_serialization = True
else:
    raise ValueError("Unknown run type " + run_type)

# setup
location = "HPC"
SetupParser.default_block = location


cb = DTKConfigBuilder.from_defaults("MALARIA_SIM",
                                    Simulation_Duration=int(365*years),
                                    Config_Name=sweep_name,
                                    Birth_Rate_Dependence="FIXED_BIRTH_RATE",
                                    Age_Initialization_Distribution_Type= "DISTRIBUTION_COMPLEX",
                                    Num_Cores=1,

                                    # interventions
                                    Enable_Default_Reporting=0,
                                    Enable_Demographics_Risk=1,
                                    Enable_Vector_Species_Report=0,
                                    )

cb.update_params({"Disable_IP_Whitelist": 1,
                  "Enable_Property_Output": 0,
                  "Enable_Spatial_Output": 1,
                  "Spatial_Output_Channels": ["Population", "Blood_Smear_Parasite_Prevalence", 'New_Infections',
                                              'New_Clinical_Cases']
                  })

if serialize:
    cb.update_params({"Serialization_Time_Steps": [365*years]})

assign_net_ip(cb, hates_net_prop)
add_health_seeking(cb, start_day=0,
                   drug=['Artemether', 'Lumefantrine'],
                   targets=[
                       {'trigger': 'NewClinicalCase', 'coverage': 0.5, 'agemin': 0, 'agemax': 5, 'seek': 1,
                        'rate': 0.3},
                       {'trigger': 'NewClinicalCase', 'coverage': 0.3, 'agemin': 5, 'agemax': 100, 'seek': 1,
                        'rate': 0.3},
                       {'trigger': 'NewSevereCase', 'coverage': 0.7, 'agemin': 0, 'agemax': 100, 'seek': 1,
                        'rate': 0.5}]
                   )


def add_intervention(cb, intervention, species_details) :

    if intervention == 'itn' :
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.3,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
          )
    elif intervention == 'llin' :
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.8,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
          )
    elif intervention == 'llin_no_disc':
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.8,
                        discard_time=365*50,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
                        )
    elif intervention == 'irs_180' :
        add_irs_group(cb, coverage=60 / 100, decay=180,
                      start_days=[365 * start for start in range(years)])
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.3,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
          )
    elif intervention == 'atsb_cdc' :
        add_atsb_by_coverage(cb, 60 / 100.,
                             killing=0.115,
                             species_list=list(species_details.keys()))
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.3,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
          )
    elif intervention == 'atsb_hlc':
        add_atsb_by_coverage(cb, 60 / 100.,
                             killing=0.0337,
                             species_list=list(species_details.keys()))
        add_annual_itns(cb, year_count=1,
                        coverage=60. / 100,
                        initial_killing=0.3,
                        start_day=5,
                        IP=[{"NetUsage": "LovesNets"}]
                        )
    elif intervention == 'atsb_alone_cdc' :
        add_atsb_by_coverage(cb, 60 / 100.,
                             killing=0.115,
                             species_list=list(species_details.keys()))
    elif intervention == 'atsb_alone_hlc':
        add_atsb_by_coverage(cb, 60 / 100.,
                             killing = 0.0337,
                             species_list = list(species_details.keys()))

    return {'intervention' : intervention}


if __name__=="__main__":

    SetupParser.init()

    # collect site-specific data to pass to builder functions
    COMPS_login("https://comps.idmod.org")
    sites = pd.read_csv("site_details.csv")

    print("finding collection ids and vector details")
    site_input_dir = os.path.join("sites", "all")

    with open("species_details.json") as f:
        species_details = json.loads(f.read())

    if asset_exp_id:
        print("retrieving asset experiment")
        asset_expt = retrieve_experiment(asset_exp_id)
        template_asset = asset_expt.simulations[0].tags
        cb.set_exe_collection(template_asset["exe_collection_id"])
        cb.set_dll_collection(template_asset["dll_collection_id"])
        cb.set_input_collection(template_asset["input_collection_id"])

    if new_inputs:
        print("generating input files")
        generate_input_files(site_input_dir, pop=2000, overwrite=True)

    # Find vector proportions for each vector in our site
    site_vectors = pd.read_csv(os.path.join(site_input_dir, "vector_proportions.csv"))
    simulation_setup(cb, species_details, site_vectors)

    # reporting
    for idx, row in site_vectors.iterrows():
        add_summary_report(cb,
                           age_bins = list(range(10, 130, 10)),
                           nodes={
                               "class": "NodeSetNodeList",
                               "Node_List": [int(row["node_id"])]
                           },
                           description = row["name"])
    # add_event_counter_report(cb, ["Bednet_Using", "Received_Vaccine"])
    # add_vector_stats_report(cb)

    if pull_from_serialization:
        print("building from pickup")

        # serialization
        print("retrieving burnin")
        expt = retrieve_experiment(burnin_id)

        df = pd.DataFrame([x.tags for x in expt.simulations])
        df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])

        df = df[df['Run_Number'] == 0]

        from_burnin_list = [
            ModFn(DTKConfigBuilder.update_params, {
                "Serialized_Population_Path": os.path.join(df["outpath"][x], "output"),
                "Serialized_Population_Filenames": [name for name in os.listdir(os.path.join(df["outpath"][x], "output")) if "state" in name],
                "Run_Number": y,
                "x_Temporary_Larval_Habitat": df["x_Temporary_Larval_Habitat"][x]})

            for x in df.index for y in range(num_runs)]

        builder = ModBuilder.from_list([
            [burnin_fn,
             ModFn(add_intervention, intervention, species_details)]
            for burnin_fn in from_burnin_list
            for intervention in interventions
        ])

    else:
        print("building burnin")
        builder = ModBuilder.from_list([[
            ModFn(DTKConfigBuilder.update_params, {
                "Run_Number": run_num,
                "x_Temporary_Larval_Habitat":10 ** hab_exp}),
        ]
            for run_num in range(10)
            for hab_exp in np.concatenate((np.arange(-3.75, -2, 0.25), np.arange(-2, 2.25, 0.1)))
            # for hab_exp in [0, 1, 2]
        ])

    run_sim_args = {"config_builder": cb,
                    "exp_name": sweep_name,
                    "exp_builder": builder}

    em = ExperimentManagerFactory.from_cb(cb)
    em.run_simulations(**run_sim_args)

