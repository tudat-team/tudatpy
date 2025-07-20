import typing
from . import RadioBase as RadioBase
from _typeshed import Incomplete
from tudatpy.estimation.observable_models_setup.links import link_definition as link_definition, receiver as receiver, reflector1 as reflector1
from tudatpy.estimation.observable_models_setup.model_settings import dsn_n_way_averaged_doppler as dsn_n_way_averaged_doppler
from tudatpy.estimation.observations import single_observation_set as single_observation_set
from tudatpy.estimation.observations_setup.ancillary_settings import dsn_n_way_doppler_ancilliary_settings as dsn_n_way_doppler_ancilliary_settings

class DerivedDopplerConverter(RadioBase):

    def extract(self, sfdu_list):
        ...

    def process(self, doppler_df, spacecraftName: Incomplete | None=None):
        ...

    def get_link_delays(self, sfdu):
        """
        Returns the transmit time tag delay, spacecraft transmit delay, and receive time tag delay for a given SFDU record.
        The secondary CHDO has to be decoded before calling this function.

        Parameters
        ----------
        sfdu : trk234.SFDU
            The SFDU record to extract the time tag delays from.

        Returns
        -------
        tuple(float, float, float)
            A tuple containing the transmit time tag delay, spacecraft transmit delay, and receive time tag delay.
            If the values are not valid or not provided, the function sets the delays to 0.
        """