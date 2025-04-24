"""
Base converter class for processing SFDU data into structured data.
"""

from abc import ABC, abstractmethod
from typing import List, Any, Union
from pandas import DataFrame


class Converter(ABC):
    @abstractmethod
    def extract(self, sfdu_list: List[Any]) -> DataFrame:
        """
        Extract data from a list of SFDU objects and return a pandas DataFrame.
        """
        pass

    @abstractmethod
    def process(self, merged_df, spacecraftName=None) -> Union[DataFrame, List[Any]]:
        """
        Process a merged DataFrame (from multiple files extract outputs) into Tudat structured format.
        For observable converters, this will be a list of
        tudatpy.numerical_simulation.estimation_setup.observation.single_observation_sets;
        For the ramp converter, a merged ramp DataFrame.

        Parameters
        ----------
        merged_df : DataFrame
            Merged DataFrame containing data from multiple files.
        spacecraftName : str, optional
            The spacecraft name used for building link definitions, if None the NAIF ID of the
            spacecraft is extracted from the tracking file.

        Returns
        -------
        Union[DataFrame, List[Any]]
            A DataFrame or a list of single observation sets.
        """
        pass
