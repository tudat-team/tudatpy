import numpy
import pybind11_stubgen.typing_ext
from ....astro import time_representation
from .... import data
from ....dynamics import environment
from ....estimation.observable_models import observables_simulation
from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
from ....estimation.observations_setup import ancillary_settings
from ....estimation.observations_setup import observations_simulation_settings
import typing
__all__ = ['ProcessedOdfFileContents', 'create_compressed_doppler_collection', 'create_odf_observed_observation_collection', 'create_pseudo_observations_and_models', 'create_tracking_txtfile_observation_collection', 'observations_from_fdets_files', 'observations_from_ifms_files', 'observations_from_multi_station_ifms_files', 'observations_from_odf_files', 'process_odf_data_multiple_files', 'process_odf_data_single_file', 'set_existing_observations', 'set_odf_information_in_bodies', 'simulate_observations', 'single_type_observation_collection']

class ProcessedOdfFileContents:
    """No documentation found."""

    def define_antenna_id(self, spacecraft_name: str, antenna_name: str) -> None:
        """
        No documentation found.
        """

    @property
    def ground_station_names(self) -> list[str]:
        """
        No documentation found.
        """

    @property
    def ignored_ground_stations(self) -> list[str]:
        """
        No documentation found.
        """

    @property
    def ignored_odf_observable_types(self) -> list[int]:
        """
        No documentation found.
        """

    @property
    def processed_observable_types(self) -> list[model_settings.ObservableType]:
        """
        No documentation found.
        """

    @property
    def raw_odf_data(self) -> list[data.OdfRawFileContents]:
        """
        No documentation found.
        """

    @property
    def start_and_end_time(self) -> tuple[float, float]:
        """
        No documentation found.
        """

def create_compressed_doppler_collection(*args, **kwargs) -> ...:
    """No documentation found."""

def create_odf_observed_observation_collection(processed_odf_file: ProcessedOdfFileContents, observable_types_to_process: list[model_settings.ObservableType], start_and_end_times_to_process: tuple[time_representation.Time, time_representation.Time]) -> ...:
    """No documentation found."""

def create_pseudo_observations_and_models(bodies: environment.SystemOfBodies, observed_bodies: list[str], central_bodies: list[str], initial_time: time_representation.Time, final_time: time_representation.Time, time_step: time_representation.Time) -> tuple[list[model_settings.ObservationSettings], ..., ..., ...]:
    """No documentation found."""

def create_tracking_txtfile_observation_collection(raw_tracking_txtfile_contents: data.TrackingTxtFileContents, spacecraft_name: str, observable_types_to_process: list[model_settings.ObservableType]=[], earth_fixed_ground_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=..., ancillary_settings: ancillary_settings.ObservationAncilliarySimulationSettings=...) -> ...:
    """No documentation found."""

def observations_from_fdets_files(ifms_file_name: str, base_frequency: float, column_types: list[str], target_name: str, transmitting_station_name: str, receiving_station_name: str, reception_band: ancillary_settings.FrequencyBands, transmission_band: ancillary_settings.FrequencyBands, earth_fixed_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ...:
    """No documentation found."""

def observations_from_ifms_files(ifms_file_names: list[str], bodies: environment.SystemOfBodies, target_name: str, ground_station_name: str, reception_band: ancillary_settings.FrequencyBands, transmission_band: ancillary_settings.FrequencyBands, apply_troposphere_correction: bool=True, earth_fixed_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ...:
    """No documentation found."""

def observations_from_multi_station_ifms_files(ifms_file_names: list[str], bodies: environment.SystemOfBodies, target_name: str, ground_station_names: list[str], reception_band: ancillary_settings.FrequencyBands, transmission_band: ancillary_settings.FrequencyBands, apply_troposphere_correction: bool=True, earth_fixed_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ...:
    """No documentation found."""

def observations_from_odf_files(bodies: environment.SystemOfBodies, odf_file_names: list[str], target_name: str, verbose_output: bool=True, earth_fixed_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ...:
    """No documentation found."""

def process_odf_data_multiple_files(file_names: list[str], spacecraft_name: str, verbose: bool=True, earth_fixed_ground_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ProcessedOdfFileContents:
    """No documentation found."""

def process_odf_data_single_file(file_name: str, spacecraft_name: str, verbose: bool=True, earth_fixed_ground_station_positions: dict[str, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]]=...) -> ProcessedOdfFileContents:
    """No documentation found."""

def set_existing_observations(observations: dict[model_settings.ObservableType, tuple[dict[links.LinkEndType, links.LinkEndId], tuple[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], list[time_representation.Time]]]], reference_link_end: links.LinkEndType, ancilliary_settings_per_observatble: dict[model_settings.ObservableType, ancillary_settings.ObservationAncilliarySimulationSettings]={}) -> ...:
    ...

def set_odf_information_in_bodies(processed_odf_file: ProcessedOdfFileContents, bodies: environment.SystemOfBodies, body_with_ground_stations_name: str='Earth', turnaround_ratio_function: typing.Callable[[ancillary_settings.FrequencyBands, ancillary_settings.FrequencyBands], float]=...) -> None:
    """No documentation found."""

def simulate_observations(simulation_settings: list[observations_simulation_settings.ObservationSimulationSettings], observation_simulators: list[observables_simulation.ObservationSimulator], bodies: environment.SystemOfBodies) -> ...:
    """Function to simulate observations.
    
    Function to simulate observations from set observation simulators and observation simulator settings.
    Automatically iterates over all provided observation simulators, generating the full set of simulated observations.
    
    
    Parameters
    ----------
    observation_to_simulate : List[ :class:`ObservationSimulationSettings` ]
        List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.
    
    observation_simulators : List[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
        List of :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.
    
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
    
    Returns
    -------
    :class:`~tudatpy.estimation.observations.ObservationCollection`
        Object collecting all products of the observation simulation."""

def single_type_observation_collection(observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition, observations_list: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], times_list: list[time_representation.Time], reference_link_end: links.LinkEndType, ancilliary_settings: ancillary_settings.ObservationAncilliarySimulationSettings=None) -> ...:
    """No documentation found."""