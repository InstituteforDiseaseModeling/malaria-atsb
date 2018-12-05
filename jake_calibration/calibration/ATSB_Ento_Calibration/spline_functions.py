import pandas as pd
import os
import json
import numpy as np

def get_spline_values(ref, spec, year=-1):
    ref = ref.reset_index()
    ref['Counts_normalized'] = ref.groupby('Channel')['Counts'].transform(lambda x: x/ x.sum())
    ref = ref[ref['Channel']==spec]
    temp_list =list(ref['Counts_normalized'])
    ref['Counts_normalized'] = temp_list
    ref['Min'] = [max([0.001, 0.01*i]) for i in temp_list]
    ref['Max'] = [max([0.01, 50*i]) for i in temp_list]
    ref['Guess'] = [0.5*i for i in temp_list]
    ref['Dynamic'] = True

    ref.reset_index(drop=True, inplace=True)
    return ref


def check_if_run(param, expname, max_diff=0.1, max_frac_diff=0.01) :

    if param['Dynamic'] == False or 'max' in param['Name']:
        return param

    try :
        LLdf = pd.read_csv(os.path.join(expname, '_plots/LL_all.csv'))
        sample = LLdf.loc[0, 'sample']
        iteration = LLdf.loc[0, 'iteration']
        with open(os.path.join(expname, 'iter%i' % iteration, 'IterationState.json')) as fin:
            iteration_state = json.loads(fin.read())['analyzers']
        if iteration_state :
            analyzer_data = next(iter(iteration_state.values()))
            df = pd.DataFrame(analyzer_data['samples'][sample])
            data = np.mean(df[df['Month'] == param['Month']]['Counts'])
            df = pd.DataFrame(analyzer_data['ref'])
            ref = np.mean(df[df['Month'] == param['Month']]['Counts'])

            if abs(ref - data) <= max_diff or (ref > 0 and abs(ref - data)/ref <= max_frac_diff) :
                print('Setting %s to Static' % param['Name'])
                param['Dynamic'] = False

    except IOError :
        pass
    return param


def get_spline_values_from_previous_run(fname, species, oldspline) :

    df = pd.read_csv(fname)
    df = df.iloc[1:]
    df = df.iloc[:-1]
    df.reset_index(drop=True, inplace=True)
    df['Month'] = oldspline['Month']
    df = df.rename(columns={'Values' : 'Guess'})
    df['Min'] = df['Guess']*0.1
    df['Max'] = df['Guess']*10
    df['Dynamic'] = True
    return df

def count_static_params(params) :

    num_params = len(params)
    num_static = 0
    for p in params :
        if not p['Dynamic'] :
            num_static += 1    
    print('%i of %i params now static' % (num_static, num_params))

def combine_counts() :
    wdir = 'C:/Users/jgerardin/Dropbox/Malaria Team Folder/projects/Mozambique/entomology_calibration/'
    fname = os.path.join(wdir, 'cluster_mosquito_counts_per_house_by_month.csv')
    adf = pd.read_csv(fname)
    catchments = adf['cluster_name'].unique()
    t = pd.DataFrame()
    for species in ['gambiae', 'funestus']:
        gdf = adf.groupby('month')[species].agg(np.sum).reset_index()
        if t.empty :
            t = gdf
        else :
            t = pd.merge(left=t, right=gdf, on='month')
    gdf = t
    gdf.sort_values(by='month', inplace=True)
    df = pd.DataFrame()
    for c in catchments :
        gdf['cluster_name'] = c
        df = pd.concat([df, gdf])
    df.to_csv(os.path.join(wdir, 'combined_%s_count.csv' % species), index=False)

# combine_counts()
