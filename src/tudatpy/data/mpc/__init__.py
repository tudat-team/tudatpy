import warnings

import numpy as np

from tudatpy.data._compat import deprecated_dir, deprecated_getattr, warn_custom_deprecation
from tudatpy.data_input.tracking_data.mpc import BatchMPC as _TrackingDataBatchMPC

from ._vfcc17 import get_weights_VFCC17

_ALIASES = {
    "load_bias_file": "tudatpy.data_input.tracking_data.optical_utilities.load_bias_file",
    "get_biases_EFCC18": "tudatpy.data_input.tracking_data.optical_utilities.get_biases_EFCC18",
}

__all__ = sorted(set(_ALIASES) | {"BatchMPC", "get_weights_VFCC17"})


def _table_with_tdb_epochs(table):
    from tudatpy.astro import time_representation

    result = table.copy()
    converter = time_representation.default_time_scale_converter()
    result["epoch_seconds_TDB"] = [
        converter.convert_time(
            input_scale=time_representation.utc_scale,
            output_scale=time_representation.tdb_scale,
            input_value=float(epoch),
        )
        for epoch in result["epoch_seconds_UTC"]
    ]
    return result


def _legacy_observation_dataset_from_table(
    table,
    station_body,
    apply_weights_vfcc17,
    apply_star_catalog_debias,
    debias_kwargs,
    in_degrees,
    included_satellites=None,
):
    from tudatpy.data_input.tracking_data.optical_utilities import (
        create_augmented_optical_table,
        get_biases_EFCC18,
    )
    from tudatpy.estimation import observations
    from tudatpy.estimation.observable_models_setup import links, model_settings

    if hasattr(table, "to_pandas"):
        table = table.to_pandas()
    processed_table = create_augmented_optical_table(table, in_degrees=in_degrees)
    processed_table = _table_with_tdb_epochs(processed_table)

    if apply_star_catalog_debias:
        right_ascension_correction, declination_correction = get_biases_EFCC18(
            processed_table, **debias_kwargs
        )
        processed_table = processed_table.assign(
            _corrected_ra=processed_table["RA"] - right_ascension_correction,
            _corrected_dec=processed_table["DEC"] - declination_correction,
        )
        value_columns = ["_corrected_ra", "_corrected_dec"]
    else:
        value_columns = ["RA", "DEC"]

    if apply_weights_vfcc17:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            processed_table = get_weights_VFCC17(mpc_table=processed_table, return_full_table=True)

    observation_dataset = observations.ObservationDataset()
    included_satellites = included_satellites or {}
    for (target, observatory), group in processed_table.groupby(
        ["number", "observatory"], sort=True
    ):
        group = group.sort_values("epoch_seconds_TDB", kind="stable")
        link_ends = {
            links.transmitter: links.body_origin_link_end_id(str(target)),
            links.receiver: (
                links.body_origin_link_end_id(included_satellites[str(observatory)])
                if str(observatory) in included_satellites
                else links.body_reference_point_link_end_id(station_body, str(observatory))
            ),
        }
        set_id = observation_dataset.add_observation_set(
            model_settings.angular_position_type,
            links.LinkDefinition(link_ends),
            [np.asarray(value) for value in group[value_columns].to_numpy()],
            group["epoch_seconds_TDB"].to_numpy(),
            links.receiver,
        )
        if apply_weights_vfcc17:
            weights = group["weight"].to_numpy()
            observation_dataset.set_weight_vector_for_set(set_id, np.ravel([weights, weights], "F"))
    return observation_dataset


class BatchMPC(_TrackingDataBatchMPC):
    """Deprecated compatibility facade for the pre-data-refactor MPC interface."""

    def __init__(self, *args, **kwargs):
        warn_custom_deprecation(
            __name__,
            "BatchMPC",
            "Use tudatpy.data_input.tracking_data.mpc.BatchMPC instead.",
        )
        super().__init__(*args, **kwargs)

    def get_observations(self, *args, **kwargs):
        result = super().get_observations(*args, **kwargs)
        self._table = _table_with_tdb_epochs(self._table)
        self._refresh_metadata()
        return result

    def create_observations_from_astropy_table(
        self,
        table,
        station_body="Earth",
        apply_weights_VFCC17=False,
        apply_star_catalog_debias=False,
        debias_kwargs=None,
        in_degrees=True,
    ):
        """Create the deprecated ObservationCollection representation."""
        from tudatpy.estimation import observations

        dataset = _legacy_observation_dataset_from_table(
            table,
            station_body,
            apply_weights_VFCC17,
            apply_star_catalog_debias,
            debias_kwargs or {},
            in_degrees,
        )
        return observations.create_observation_collection_from_dataset(dataset)

    def to_tudat_observation_collection(
        self,
        bodies,
        included_satellites=None,
        station_body="Earth",
        add_sbdb_gravity_model=False,
        apply_weights_VFCC17=True,
        apply_star_catalog_debias=True,
        debias_kwargs=None,
    ):
        """Convert the batch to the deprecated ObservationCollection representation."""
        from tudatpy.dynamics import environment_setup
        from tudatpy.dynamics.environment_setup import ground_station
        from tudatpy.estimation import observations

        station_body_object = bodies.get(station_body)
        existing_stations = set(station_body_object.ground_station_list)
        required_stations = set(self.observatories) - set((included_satellites or {}).keys())
        for station_settings in ground_station.optical_telescope_stations():
            if (
                station_settings.station_name in required_stations
                and station_settings.station_name not in existing_stations
            ):
                environment_setup.add_ground_station(station_body_object, station_settings)

        for target in self.MPC_objects:
            if not bodies.does_body_exist(str(target)):
                bodies.create_empty_body(str(target))
                self._bodies_created[str(target)] = "empty body"
                if add_sbdb_gravity_model:
                    from tudatpy.dynamics.environment_setup import add_gravity_field_model
                    from tudatpy.dynamics.environment_setup.gravity_field import central_sbdb

                    add_gravity_field_model(bodies, str(target), central_sbdb(str(target)))
                    self._bodies_created[str(target)] = "empty body + sbdb gravity"

        dataset = _legacy_observation_dataset_from_table(
            self.table,
            station_body,
            apply_weights_VFCC17,
            apply_star_catalog_debias,
            debias_kwargs or {},
            False,
            included_satellites,
        )
        return observations.create_observation_collection_from_dataset(dataset)


def __getattr__(name):
    return deprecated_getattr(__name__, _ALIASES, name)


def __dir__():
    return deprecated_dir(globals(), _ALIASES)
