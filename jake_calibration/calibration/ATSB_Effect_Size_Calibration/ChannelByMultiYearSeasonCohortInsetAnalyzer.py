import calendar
import datetime
import logging
from abc import abstractmethod
import pandas as pd
import numpy as np
from scipy.stats import binom

from calibtool import LL_calculators
from dtk.utils.parsers.malaria_summary import summary_channel_to_pandas
from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer, thread_lock
from calibtool.analyzers.Helpers import \
    convert_annualized, convert_to_counts, age_from_birth_cohort, aggregate_on_index, aggregate_on_month

logger = logging.getLogger(__name__)


class ChannelByMultiYearSeasonCohortInsetAnalyzer(BaseCalibrationAnalyzer):
    """
    Base class implementation for similar comparisons of age-binned reference data to simulation output.
    """

    filenames = ['output/ReportMalariaFiltered.json', 'config.json']
    population_channel = 'Statistical Population'
    hbr_channel = 'Daily Bites per Human'

    site_ref_type = 'entomology_by_season'

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, site, weight=1, compare_fn=LL_calculators.euclidean_distance_pandas, **kwargs):
        super(ChannelByMultiYearSeasonCohortInsetAnalyzer, self).__init__(site, weight, compare_fn)
        self.reference = site.get_reference_data(self.site_ref_type)
        self.duration = kwargs['duration']
        self.site_name = site.name

    def apply(self, parser):
        """
        Extract data from output data and accumulate in same bins as reference.
        """

        cfg = parser.raw_data[self.filenames[1]]['parameters']
        species = cfg['Vector_Species_Names'][0]
        hfca = self.site_name
        multiplier = 1

        # Load data from simulation
        data = { x : parser.raw_data[self.filenames[0]]['Channels'][x]['Data'] for x in [self.population_channel, self.hbr_channel]}

        data = pd.DataFrame(data)
        data['Time'] = data.index
        data['Species'] = species
        data = data.rename(columns={ self.population_channel : 'Population',
                                     self.hbr_channel : 'HBR'})
        # data['Vector_per_Human'] = data['VectorPopulation'] / data['Population'] * multiplier

        data['Month'] = data['Time'].apply(lambda x: self.monthparser(x % 365))
        data['Year'] = data['Time'].apply(lambda x: int((x-1)/365))
        data['Month'] += data['Year'] * 12
        # data = data.groupby(['Month', 'Species'])['Vector_per_Human'].apply(np.mean).reset_index()

        data = data.rename(columns={'HBR': 'Counts', 'Species': 'Channel'})
        data = data.sort_values(['Channel', 'Month'])
        data = data.set_index(['Channel', 'Month'])
        channel_data_dict = {}

        for channel in self.site.metadata['species']:

            with thread_lock:  # TODO: re-code following block to ensure thread safety (Issue #758)?

                # Reset multi-index and perform transformations on index columns
                df = data.copy().reset_index()
                df = df.rename(columns={'Counts': channel})
                del df['Channel']

                # Re-bin according to reference and return single-channel Series
                rebinned = aggregate_on_index(df, self.reference.loc(axis=1)[channel].index, keep=[channel])
                channel_data_dict[channel] = rebinned[channel].rename('Counts')

        sim_data = pd.concat(channel_data_dict.values(), keys=channel_data_dict.keys(), names=['Channel'])
        sim_data = pd.DataFrame(sim_data)  # single-column DataFrame for standardized combine/compare pattern
        sim_data.sample = parser.sim_data.get('__sample_index__')
        sim_data.sim_id = parser.sim_id

        return sim_data

    @classmethod
    def plot_comparison(cls, fig, data, **kwargs):
        axs = fig.axes
        # ax = fig.gca()
        df = pd.DataFrame.from_dict(data, orient='columns')
        nrows, ncols = 1, len(df.Channel.unique())
        fmt_str = kwargs.pop('fmt', None)
        args = (fmt_str,) if fmt_str else ()
        if not axs:
            fig.set_size_inches((12, 6))  # override smaller single-panel default from SiteDataPlotter
            axs = [fig.add_subplot(nrows, ncols, iax+1) for iax in range(nrows*ncols)]
        for iax, (species, group_df) in enumerate(df.groupby('Channel')):
            ax = axs[iax]
            counts = list(group_df['Counts'])
            time = list(group_df['Month'])
            months = [calendar.month_abbr[i%12] for i in time]
            irow, icol = int(iax / ncols), (iax % ncols)
            if irow == 0:
                ax.set_title(species)
                ax.set_xlabel('Month')
            if icol==0:
                ax.set_ylabel('Human Biting Rate')
            if 'reference' in kwargs:
                kwargs2 = kwargs.copy()
                kwargs2.pop('reference')
                ax.plot(time, counts, *args, **kwargs2)
            else:
                ax.plot(time, counts, *args, **kwargs)
            ax.xaxis.set_ticks(time)
            ax.set_xticklabels(months)