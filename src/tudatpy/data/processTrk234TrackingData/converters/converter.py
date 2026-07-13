"""
Base converter class for processing SFDU data into structured data.
"""

from abc import ABC, abstractmethod
from pandas import DataFrame
from trk234 import SFDU
from tudatpy.data_access.tracking import TrackingData


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
    ) -> list[TrackingData]:
        """
        Process a merged DataFrame (from multiple files extract outputs) into Tudat structured format.
        For observable converters, this will be a list of
        :class:`~tudatpy.data.TrackingData` objects.

        Parameters
        ----------
        merged_df : pandas.DataFrame
            Merged DataFrame containing data from multiple files.
        spacecraftName : str, optional
            The spacecraft name used for building link definitions, if None the NAIF ID of the
            spacecraft is extracted from the tracking file.

        Returns
        -------
        list[TrackingData]
            A list of tracking data objects.
        """
        pass
