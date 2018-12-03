import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from input_file_generation.add_properties_to_demographics import generate_demographics_properties
from sim_output_processing.createSimDirectoryMap import createSimDirectoryMap

userpath = 'D:/'

wdir = os.path.join(userpath, 'Dropbox (IDM)', 'Malaria Team Folder/projects/atsb')
datadir = os.path.join(wdir, 'sim_data')

def sample_LHC(paramdict, numsamples, output_fname='', force=False) :

    if output_fname and not force :
        try :
            df = pd.read_csv(output_fname)
        except FileNotFoundError :
            force = True

    if force :
        df = pd.DataFrame()

        for key, val in paramdict.items() :
            if 'log' in val and val['log'] :
                samples = np.logspace(np.log10(val['min']), np.log10(val['max']), numsamples)
            else :
                samples = np.linspace(val['min'], val['max'], numsamples)
            np.random.shuffle(samples)
            df[key] = samples

        if output_fname :
            df.to_csv(output_fname, index=False)
    return df

def generate_samples(output_fname) :

    paramdict = {
        'feeding rate' : { 'min' : 0,
                           'max' : 0.5,
                           'log' : False},
        'larval hab' : { 'min' : 0.05,
                         'max' : 100,
                         'log' : True}
    }

    df = sample_LHC(paramdict, 1000, output_fname=output_fname)

def set_up_cohort_IPs() :

    inputdir = os.path.join(wdir, 'Input')
    demo_fname = os.path.join(inputdir, 'demographics.json')
    overlay_fname = os.path.join(inputdir, 'demographics_cohort_IP_overlay.json')

    generate_demographics_properties(demo_fname, overlay_fname, as_overlay=True,
                                     df=pd.read_csv(os.path.join(datadir, 'cohort_IP_setup.csv')))


def add_serialization_paths_to_hab_csv(fname, expid) :

    df = pd.read_csv(fname)
    df = df.round({'larval hab' : 6})
    sdf = createSimDirectoryMap(expid, is_name=False)
    sdf = sdf.drop(columns=[x for x in sdf.columns.values if 'id' in x])
    sdf = sdf.round({'x_Temporary_Larval_Habitat' : 6})
    df = pd.merge(left=df, right=sdf, left_on='larval hab', right_on='x_Temporary_Larval_Habitat')
    df = df.rename(columns={'outpath' : 'serialized file dir'})
    df.to_csv(fname, index=False)

if __name__ == '__main__' :

    expid = '58f94e3f-9aa0-e811-a2c0-c4346bcb7275'

    output_fname = os.path.join(datadir, 'LH_burnin.csv')
    # add_serialization_paths_to_hab_csv(output_fname, expid)

    df = createSimDirectoryMap(expid, is_name=False)
    df.to_csv(output_fname, index=False)