import itertools
import random
from datetime import date, datetime
import calendar
import logging
from collections import OrderedDict
import numpy.ma as ma
import json

import pandas as pd
import numpy as np


logger = logging.getLogger(__name__)



def convert_annualized(s, reporting_interval=None, start_day=None):
    """
    Use Time index level to revert annualized rate channels,
    e.g. Annual Incidence --> Incidence per Reporting Interval
         Average Population --> Person Years
    :param s: pandas.Series with 'Time' level in multi-index
    :param reporting_interval: (optional) metadata from original output file on interval related to 'Time' index
    :param start_day: (optional) metadata from original output file on beginning of first aggregation interval
    :return: pandas.Series normalized to Reporting Interval as a fraction of the year
    """

    s_ix = s.index
    time_ix = s_ix.names.index('Time')
    time_levels = s_ix.levels[time_ix].values

    start_time = (start_day - 1) if start_day else 0  # metadata reported at end of first time step

    time_intervals = np.diff(np.insert(time_levels, 0, values=start_time))  # prepending simulation Start_Time

    if reporting_interval and len(time_intervals) > 1 and np.abs(time_intervals[0] - reporting_interval) > 1:
        raise Exception('Time differences between reports differ by more than integer rounding.')
    if len(time_intervals) > 2 and np.abs(time_intervals[0] - time_intervals[1]) > 1:
        raise Exception('Time differences between reports differ by more than integer rounding.')

    intervals_by_time = dict(zip(time_levels, time_intervals))
    df = s.reset_index()
    df[s.name] *= df.Time.apply(lambda x: intervals_by_time[x] / 365.0)
    return df.set_index(s_ix.names)


def convert_to_counts(rates, pops):
    """
    Convert population-normalized rates to counts
    :param rates: a pandas.Series of normalized rates
    :param pops: a pandas.Series of average population counts
    :return: a pandas.Series (same binning as rates)
    """
    rate_idx = rates.index.names
    pop_idx = pops.index.names

    # Join rates to population counts on the binning of the latter
    df = rates.reset_index().set_index(pop_idx)\
              .join(pops, how='left')\
              .reset_index().set_index(rate_idx)

    counts = (df[rates.name] * df[pops.name]).rename(rates.name)
    return counts


def age_from_birth_cohort(df):
    """
    Reinterpret 'Time' as 'Age Bin' for a birth cohort
    :param df: a pandas.DataFrame of counts and 'Time' in days
    :return: a pandas.DataFrame including an additional (or overwritten) 'Age Bin' column
    """

    df['Age Bin'] = df['Time'] / 365.0   # Time in days but Age in years
    return df


def aggregate_on_index(df, index, keep=slice(None)):
    """
    Aggregate and re-index data on specified (multi-)index (levels and) intervals
    :param df: a pandas.DataFrame with columns matching the specified (Multi)Index (level) names
    :param index: pandas.(Multi)Index of categorical values or right-bin-edges, e.g. ['early', 'late'] or [5, 15, 100]
    :param keep: optional list of columns to keep, default=all
    :return: pandas.Series or DataFrame of specified channels aggregated and indexed on the specified binning
    """

    if isinstance(index, pd.MultiIndex):
        levels = index.levels
    else:
        levels = [index]  # Only one "level" for Index. Put into list for generic pattern as for MultiIndex

    for ix in levels:
        logger.debug("%s (%s) : %s" % (ix.name, ix.dtype, ix.values))

        if ix.name not in df.columns:
            raise Exception('Cannot perform aggregation as MultiIndex level (%s) not found in DataFrame:\n%s' % (ix.name, df.head()))

        # If dtype is object, these are categorical (e.g. season='start_wet', channel='gametocytemia')
        if ix.dtype == 'object':
            # TODO: DatetimeIndex, TimedeltaIndex are implemented as int64.  Do we want to catch them separately?
            # Keep values present in reference index-level values; drop any that are not
            df = df[df[ix.name].isin(ix.values)]

        # Otherwise, the index-level values are upper-edges of aggregation bins for the corresponding channel
        elif ix.dtype in ['int64', 'float64']:
            bin_edges = np.concatenate(([-np.inf], ix.values))
            labels = ix.values  # just using upper-edges as in reference
            # TODO: more informative labels? would need to be modified also in reference to maintain useful link...
            # labels = []
            # for low, high in pairwise(bin_edges):
            #     if low == -np.inf:
            #         labels.append("<= {0}".format(high))
            #     elif high == np.inf:
            #         labels.append("> {0}".format(low))
            #     else:
            #         labels.append("{0} - {1}".format(low, high))

            df[ix.name] = pd.cut(df[ix.name], bin_edges, labels=labels)

        else:
            logger.warning('Unexpected dtype=%s for MultiIndex level (%s). No aggregation performed.', ix.dtype, ix.name)

    # Aggregate on reference MultiIndex, keeping specified channels and dropping missing data
    if keep != slice(None):
        df = df.groupby([ix.name for ix in levels]).sum()[keep].dropna()
    logger.debug('Data aggregated/joined on MultiIndex levels:\n%s', df.head(15))
    return df


def aggregate_on_month(sim, ref):
    months = list(ref['Month'].unique())
    sim = sim[sim['Month'].isin(months)]

    return sim


def ento_data(csvfilename, metadata):

    df = pd.read_csv(csvfilename)

    df = df[['date', 'gambiae_count', 'funestus_count', 'adult_house']]
    df['gambiae'] = df['gambiae_count'] / df['adult_house']
    df['funestus'] = df['funestus_count'] / df['adult_house']
    df = df.dropna()

    df['date'] = pd.to_datetime(df['date'])
    dateparser = lambda x: int(x.strftime('%m'))
    df['Month'] = df['date'].apply(lambda x: int(dateparser(x)))
    df2 = df.groupby('Month')['gambiae'].apply(np.mean).reset_index()
    df2['funestus'] = list(df.groupby('Month')['funestus'].apply(np.mean))

    # Keep only species requested
    for spec in metadata['species']:
        df1 = df2[['Month', spec]]
        df1 = df1.rename(columns={spec: 'Counts'})
        df1['Channel'] = [spec] * len(df1)
        if 'dftemp' in locals():
            dftemp = pd.concat([dftemp, df1])
        else:
            dftemp = df1.copy()

    dftemp = dftemp.sort_values(['Channel', 'Month'])
    dftemp = dftemp.set_index(['Channel', 'Month'])

    return dftemp


def multi_year_ento_data(csvfilename, metadata):

    df = pd.read_csv(csvfilename)
    if metadata['HFCA']:
        df = df[df['hf_name']==metadata['HFCA']]

    df = df[['date', 'gambiae_count', 'funestus_count', 'adult_house']]
    df['gambiae'] = df['gambiae_count'] / df['adult_house']
    df['funestus'] = df['funestus_count'] / df['adult_house']
    df = df.dropna()

    df['date'] = pd.to_datetime(df['date'])
    mask = (df['date'] <= '2017-12-31')
    df = df.loc[mask]
    dateparser = lambda x: int(x.strftime('%m'))
    df['Month'] = df['date'].apply(lambda x: int(dateparser(x)))
    df['Year'] = df['date'].apply(lambda x: int(x.strftime('%y')))
    df['Year'] = df['Year'].apply(lambda x: x % min(df['Year']))
    df['Month'] += df['Year'] * 12
    df2 = df.groupby(['Month'])['gambiae'].apply(np.mean).reset_index()
    df2['funestus'] = list(df.groupby(['Month'])['funestus'].apply(np.mean))

    # Keep only species requested
    for spec in metadata['species']:
        df1 = df2[['Month', spec]]
        df1 = df1.rename(columns={spec: 'Counts'})
        df1['Channel'] = [spec] * len(df1)
        if 'dftemp' in locals():
            dftemp = pd.concat([dftemp, df1])
        else:
            dftemp = df1.copy()

    dftemp = dftemp.sort_values(['Channel', 'Month'])
    dftemp = dftemp.set_index(['Channel', 'Month'])

    return dftemp


def multi_year_ento_data_clustered(metadata):

    df = pd.read_csv(metadata['csvfilename'])

    df = df[['month', 'gambiae', 'funestus']]

    df['date'] = pd.to_datetime(df['month'])
    mask = (df['date'] <= '2017-12-31')
    df = df.loc[mask]
    dateparser = lambda x: int(x.strftime('%m'))
    df['Month'] = df['date'].apply(lambda x: int(dateparser(x)))
    df['Year'] = df['date'].apply(lambda x: int(x.strftime('%y')))
    df['Year'] = df['Year'].apply(lambda x: x % min(df['Year']))
    df['Month'] += df['Year'] * 12
    df2 = df.groupby(['Month'])['gambiae'].apply(np.mean).reset_index()
    df2['funestus'] = list(df.groupby(['Month'])['funestus'].apply(np.mean))

    # Keep only species requested
    for spec in metadata['species']:
        df1 = df2[['Month', spec]]
        df1 = df1.rename(columns={spec: 'Counts'})
        df1['Channel'] = [spec] * len(df1)
        if 'dftemp' in locals():
            dftemp = pd.concat([dftemp, df1])
        else:
            dftemp = df1.copy()

    dftemp = dftemp.sort_values(['Channel', 'Month'])
    dftemp = dftemp.set_index(['Channel', 'Month'])

    return dftemp


