from tudatpy.data._compat import deprecated_dir, deprecated_getattr, warn_custom_deprecation
from tudatpy.data_input.tracking_data.tnf import (
    TnfTrackingDataProcessor as _TnfTrackingDataProcessor,
)

_ALIASES = {}

__all__ = ["Trk234Processor"]


class Trk234Processor(_TnfTrackingDataProcessor):
    """Deprecated compatibility facade for the former TNF processor."""

    def __init__(self, *args, **kwargs):
        warn_custom_deprecation(
            __name__,
            "Trk234Processor",
            "Use tudatpy.data_input.tracking_data.tnf.TnfTrackingDataProcessor instead.",
        )
        super().__init__(*args, **kwargs)

    def process_observation_collection(self):
        from tudatpy.dynamics import environment_setup
        from tudatpy.dynamics.environment_setup import ground_station
        from tudatpy.estimation import observations

        tracking_data, _ = self.process()
        body_settings = environment_setup.get_default_body_settings(["Earth"], "SSB", "J2000")
        body_settings.get("Earth").ground_station_settings = ground_station.dsn_stations()
        bodies = environment_setup.create_system_of_bodies(body_settings)
        dataset = observations.create_observation_dataset_from_tracking_data(tracking_data, bodies)
        return observations.create_observation_collection_from_dataset(dataset)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
