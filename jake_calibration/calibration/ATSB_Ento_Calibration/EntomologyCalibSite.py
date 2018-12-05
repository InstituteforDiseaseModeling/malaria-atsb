import logging
from abc import ABCMeta

from calibtool.CalibSite import CalibSite
from malaria.study_sites.site_setup_functions import \
    config_setup_fn, vector_stats_report_fn


logger = logging.getLogger(__name__)


class EntomologyCalibSite(CalibSite):
    """
    An abstract class that implements the simulation setup for density-by-season-and-age analyses:
    - Laye, Burkina Faso
    - Dapelogo, Burkina Faso
    """

    __metaclass__ = ABCMeta

    metadata = {
        'months': {}
    }

    def get_setup_functions(self):
        return [
            config_setup_fn(duration=21915),  # 60 years (with leap years)
            vector_stats_report_fn()
        ]

    def get_reference_data(self, reference_type):
        site_ref_type = 'entomology_by_season'

        if reference_type is not site_ref_type:
            raise Exception("%s does not support %s reference_type, only %s.",
                            self.__class__.__name__, reference_type, site_ref_type)

    # def get_analyzers(self):
    #     return [ChannelBySeasonCohortAnalyzer(site=self, seasons=self.metadata['months'])]
