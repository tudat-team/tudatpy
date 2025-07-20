import typing
from _typeshed import Incomplete
from abc import ABC, abstractmethod
from pandas import DataFrame
from typing import Any

class Converter(ABC):

    @abstractmethod
    def extract(self, sfdu_list: list[Any]) -> DataFrame:
        """
        Extract data from a list of SFDU objects and return a pandas DataFrame.
        """

    @abstractmethod
    def process(self, merged_df, spacecraftName: Incomplete | None=None) -> DataFrame | list[Any]:
        """
        Process a merged DataFrame (from multiple files extract outputs) into Tudat structured format.
        For observable converters, this will be a list of
        tudatpy.estimation.observations.single_observation_sets;
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