import typing
from _typeshed import Incomplete
from tudatpy.astro import time_representation as time_representation
from tudatpy.dynamics.environment import PiecewiseLinearFrequencyInterpolator as PiecewiseLinearFrequencyInterpolator
from tudatpy.estimation.observations import ObservationCollection as ObservationCollection

class Trk234Processor:
    """Processor for TNF files using pytrk234.
    
    For a given set of requested observables types (e.g. [\'doppler\', \'range\']),
    this processor iterates file-by-file, uses each converter’s extract method to obtain per-file data,
    merges the outputs, and then calls each converter’s process method to produce final observation
    collection. If simulation bodies are provided, ramp data are processed and used to set the
    stations\' frequency interpolator with the set_tnf_information_in_bodies() method.
    
    Examples
    --------
    >>> from tudatpy.data import Trk234Processor
    >>>
    >>> # Define TNF file paths
    >>> tnf_files = ["mro_kernels/mromagr2012_002_1426xmmmv1.tnf"]
    >>>
    >>> # Create processor for both Doppler and range data
    >>> tnf_processor = Trk234Processor(
    ...     tnf_files,
    ...     ["doppler", "range"],
    ...     spacecraft_name="MRO"
    ... )
    >>>
    >>> # Process observations
    >>> observations = tnf_processor.process()
    >>>
    >>> # Set frequency information in the bodies assuming you have a bodies object tudatpy.dynamics.environment.SystemOfBodies
    >>> tnf_processor.set_tnf_information_in_bodies(bodies)"""
    tnf_file_paths: Incomplete
    spacecraft_name: Incomplete
    bodies: Incomplete
    converters: Incomplete
    ramp_converter: Incomplete

    def __init__(self, tnf_file_paths, requested_types, spacecraft_name: Incomplete | None=None, bodies: Incomplete | None=None) -> None:
        """
        Parameters
        ----------
        tnf_file_paths : list[str]
            List of TNF file paths to be processed.
        requested_types : list[str]
            List of requested radiometric data types, e.g., [\'doppler\', \'range\'].
            Note: "ramp" should NOT be included here.
        spacecraft_name : str, optional
            The spacecraft name for building link definitions.
        bodies : object, optional
            The simulation bodies container (used for setting ramp data).
        """

    def process(self):
        """
        Process all TNF files provided at initialization. For each file, decode the SFDU data,
        and for each requested radiometric data type, extract data via the converter's extract method.
        Ramp data is also extracted separately if a ramp converter is available.
        Then, merge the per-file outputs and process them to produce the final outputs.

        Returns
        -------
        ObservationCollection or None
            An ObservationCollection containing all radiometric observation sets,
            or None if none were processed.
        """

    def set_tnf_information_in_bodies(self, bodies):
        """
        Update stations in bodies by setting the frequency interpolators from the ramp data.
        Set the transponder turnaround ratio for the spacecraft.

        NOTE: It's not optimal to set the transponder turnaround ratio here, but it's done for now

        Parameters
        ----------
        ramp_df : pd.DataFrame
            DataFrame containing ramp data.
        spacecraft_name : str
            The spacecraft name used for setting the transponder turnaround ratio.
        bodies : object
            The simulation bodies container.
        """