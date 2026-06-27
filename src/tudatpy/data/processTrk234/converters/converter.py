"""
Base converter class for processing SFDU data into structured data.
"""

from abc import ABC, abstractmethod
from pandas import DataFrame
from trk234 import SFDU
from tudatpy.estimation.observations import ObservationDataset


class Converter(ABC):
    @abstractmethod
    def extract(self, sfdu_list: list[SFDU]) -> DataFrame:
        """
        Extract data from a list of SFDU objects and return a pandas DataFrame.
        """
        pass

    @abstractmethod
    def process(
        self, merged_df: DataFrame, spacecraftName: str | None = None
    ) -> ObservationDataset:
        """
        Process a merged DataFrame (from multiple files extract outputs) into Tudat structured format.
        For observable converters, this will be an
        :class:`~tudatpy.estimation.observations.ObservationDataset` object.

        Parameters
        ----------
        merged_df : pandas.DataFrame
            Merged DataFrame containing data from multiple files.
        spacecraftName : str, optional
            The spacecraft name used for building link definitions, if None the NAIF ID of the
            spacecraft is extracted from the tracking file.

        Returns
        -------
        ObservationDataset
            A dataset containing the converted observation sets.
        """
        pass
