import numpy as np
import pandas as pd
from tudatpy.estimation import observations, observable_models_setup
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.plotting._observation_collection import plot_sky_map

class TudatObservationWrapper:
    """
    A Python wrapper for the Tudatpy C++ ObservationCollection.
    Provides convenience methods like visualization and data extraction.
    """
    def __init__(self, collection: observations.ObservationCollection):
        self.collection = collection

    @property
    def object_designation(self) -> str:
        """
        Extracts the transmitter body name from the link definitions.
        Based on expose_observations.cpp structure.
        """
        link_definitions = self.collection.link_definition_ids

        # Get the first available link definition dictionary
        if not link_definitions:
            return "Unknown"

        first_definition = next(iter(link_definitions.values()))

        # Access the transmitter LinkEndId
        transmitter_type = observable_models_setup.links.transmitter
        if transmitter_type in first_definition:
            # Returns the string name of the body (e.g., '433')
            return first_definition[transmitter_type].body_name

        return "Unknown"

    def visualize(self, projection: str = 'mollweide', figsize: tuple = (14, 7)):
        """
        Delegates to the plotting module to visualize the collection.
        """
        return plot_sky_map(self.collection, self.object_designation, projection, figsize)

class UniversalObservationAdapter:
    """
    Class that transforms structured objects (like dataframes) into
    TudatObservationWrappers.

    This can in principle be used within any MPC, Horizons, NEOCC
    and all those data APIs that do not support C++ native tudatpy functions.

    Also useful if users have external CSV that they want to easily convert to tudat-compatible objects.
    """

    def __init__(self):
        self.time_converter = time_representation.default_time_scale_converter()

    def from_astropy_table(self):
        pass
    def from_whatever_we_like(self):
        pass
    def from_df(self, df: pd.DataFrame, designator: str, center: str) -> TudatObservationWrapper:
        """
        Main entry point: Converts a DataFrame (for instance, NEOCC format) to a wrapped collection.
        """
        # 1. Coordinate Conversion (HMS/DMS -> Radians)
        radec_values = self._convert_coordinates(df)

        # 2. Time Conversion (UTC -> TDB)
        times_tdb = self._convert_times(df.iloc[:, 0])

        # 3. Setup Link Ends
        link_ends = self._setup_link_ends(designator, center)

        # 4. Create the C++ Observation Set
        single_observation_set = observations.create_single_observation_set(
            observable_models_setup.model_settings.angular_position_type,
            link_ends,
            radec_values,
            times_tdb,
            observable_models_setup.links.receiver
        )

        # 5. Build C++ Collection and wrap it
        cpp_collection = observations.ObservationCollection([single_observation_set])
        return TudatObservationWrapper(cpp_collection)

    def _convert_coordinates(self, df: pd.DataFrame) -> np.ndarray:
        """Vectorized HMS/DMS to Radians math from working script."""
        # Right Ascension conversion
        ra_h, ra_m, ra_s = df.iloc[:, 2].astype(float), df.iloc[:, 3].astype(float), df.iloc[:, 4].astype(float)
        ra_deg = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0
        ra_rad = (np.deg2rad(ra_deg) + np.pi) % (2 * np.pi) - np.pi

        # Declination conversion
        dec_d, dec_m, dec_s = df.iloc[:, 5].astype(float), df.iloc[:, 6].astype(float), df.iloc[:, 7].astype(float)
        sign = np.copysign(1, dec_d)
        dec_deg = sign * (np.abs(dec_d) + dec_m / 60.0 + dec_s / 3600.0)
        dec_rad = np.deg2rad(dec_deg)

        return np.column_stack((ra_rad, dec_rad))

    def _convert_times(self, utc_series: pd.Series) -> list:
        """Converts UTC datetime column to TDB seconds since J2000."""
        tdb_epochs = []
        for date in utc_series:
            utc_epoch = DateTime.from_python_datetime(date).to_epoch()
            tdb_epoch = self.time_converter.convert_time(
                time_representation.utc_scale,
                time_representation.tdb_scale,
                utc_epoch
            )
            tdb_epochs.append(tdb_epoch)
        return tdb_epochs

    def _setup_link_ends(self, designator: str, center: str) -> dict:
        """Determines link end IDs based on the observation center."""
        norm_center = '500' if center.lower() == "earth" else center

        if norm_center == '500':
            return {
                observable_models_setup.links.transmitter:
                    observable_models_setup.links.body_origin_link_end_id(designator),
                observable_models_setup.links.receiver:
                    observable_models_setup.links.body_origin_link_end_id("Earth")
            }
        else:
            return {
                observable_models_setup.links.transmitter:
                    observable_models_setup.links.body_origin_link_end_id(designator),
                observable_models_setup.links.receiver:
                    observable_models_setup.links.body_reference_point_link_end_id("Earth", norm_center)
            }