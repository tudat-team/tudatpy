from pathlib import Path
import numpy as np
import pandas as pd
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch
from tudatpy.dynamics.environment_setup.shape import polyhedron
from tudatpy.dynamics.environment_setup import rotation_model
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from warnings import warn
from tudatpy.resource_installer import download_extra_file, DEFAULT_DEST, DEFAULT_CACHE
import tarfile
import sys
from functools import cache

# Location of DAMIT archive files under tudat/resources
DAMIT_RESOURCES = DEFAULT_DEST / "damit"


@cache
def _damit_files_path() -> Path:
    """Find the path to the DAMIT catalog"""
    # Lookup files in damit resources
    if DAMIT_RESOURCES.is_dir():
        for path in DAMIT_RESOURCES.iterdir():
            if path.is_dir() and path.name.startswith("damit-"):
                return path

    # Attempt download if not present
    print("No DAMIT files found on disk.")
    _download_latest_damit_models()

    return _damit_files_path()


def _download_latest_damit_models() -> None:
    """
    Downloads the latest DAMIT shape models from the DAMIT website
    """
    while True:
        user_query = input(
            "The DAMIT archive (~4 GB) will be downloaded from the internet. Proceed? (y/n) "
        )
        if user_query.lower() == "n":
            sys.exit()
        if user_query.lower() == "y":
            break

    DAMIT_RESOURCES.mkdir(parents=True, exist_ok=True)
    url = "https://damit.cuni.cz/projects/damit/exports/complete/latest"

    tar_path = download_extra_file(
        url=url, dest_path=DAMIT_RESOURCES / "damit.tar.gz", cache_dir=DEFAULT_CACHE, force=False
    )

    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(DAMIT_RESOURCES, filter="data")

    tar_path.unlink()  # Remove compressed file

    print(f"DAMIT files saved to {str(DAMIT_RESOURCES)}")


def _load_table() -> pd.DataFrame:
    """Loads all asteroid shape models into a pandas dataframe"""
    # Asteroid data table containing primarily the DAMIT ID's of the asteroids
    asteroid_data = pd.read_csv(
        _damit_files_path() / "tables" / "asteroids.csv", dtype={"number": "Int64"}
    ).rename(columns={"number": "mpc_number", "id": "asteroid_id"})
    # Model data containing actual models
    asteroid_models = pd.read_csv(
        _damit_files_path() / "tables" / "asteroid_models.csv",
    )

    # Append the MPC number to the table to query by MPC
    asteroid_models = pd.merge(
        asteroid_data[["mpc_number", "asteroid_id"]],
        asteroid_models,
        on="asteroid_id",
        how="inner",
    ).set_index("id", drop=False)

    # Convert units to SI
    asteroid_models[["lambda", "beta", "phi0"]] = np.deg2rad(
        asteroid_models[["lambda", "beta", "phi0"]]
    )
    asteroid_models["period"] = asteroid_models["period"] * 3600
    asteroid_models["epoch"] = asteroid_models["jd0"].transform(julian_day_to_seconds_since_epoch)
    asteroid_models[["equiv_diameter", "diameter_err"]] = (
        1e3 * asteroid_models[["equiv_diameter", "diameter_err"]]
    )
    return asteroid_models


def damit_model_query(
    mpc_number: int,
) -> pd.DataFrame:
    """
    Returns a table asteroid model data for the queried MPC number. An asteroid may have one or more models,
    or none, in which case the table is returned empty.

    Parameters
    ----------
    mpc_number : int
        MPC number to query

    Returns
    -------
    pd.DataFrame
        A dataframe containing the models that match the query
    """
    table = _load_table()
    query_result = table[table["mpc_number"] == mpc_number]
    return query_result


class DamitModel:
    """
    Interface with DAMIT database for asteroid shape models and rotational parameters
    """

    def __init__(
        self,
        mpc_number: int | None = None,
        damit_model_id: int | None = None,
    ):
        """
        Initialize a DamitModel object

        Parameters
        ----------
        mpc_number : int, optional
            If given, queries a DamitModel by small-body MPC number. The model with the highest quality flag is selected.
            If there is no highest flag, selection is arbitrary
        damit_model_id : int, optional
            If given, queries a DamitModel by the Damit model Id
        """
        # Query by mpc number and choose best model
        if mpc_number and not damit_model_id:
            best_model = damit_model_query(
                mpc_number,
            )
            if best_model.empty:
                raise LookupError(f"No shape models found for {mpc_number}")
            self._data = best_model.sort_values("quality_flag").iloc[-1]

        # Retrieve by model ID
        elif damit_model_id and not mpc_number:
            try:
                self._data = _load_table().loc[damit_model_id]
            except KeyError as e:
                raise LookupError(f"No shape model found with ID {damit_model_id}") from e

        else:
            self._data = None
            raise ValueError("Must specify either a valid MPC number or Damit Model ID")

    @property
    def _path_to_shape_model(self) -> Path:
        """Returns the path to the shape model file"""
        shape_model_file = (
            _damit_files_path()
            / "files"
            / f'asteroid_{self._data["asteroid_id"]}'
            / f"model_{self._data.name}"
            / "shape.txt"
        )
        return shape_model_file

    def polyhedron_rotational_ephemeris_settings(self, base_frame: str):
        """
        Create rotational ephemeris settings of the polyhedron shape model

        Returns
        -------

        """
        if base_frame not in ("J2000", "ECLIPJ2000"):
            raise ValueError("`base_frame` must be J2000 or ECLIPJ2000")

        # Rotation matrices
        rotation_y = lambda a: np.array(
            [[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]]
        )
        rotation_z = lambda a: np.array(
            [[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]
        )

        # Initial orientation ecliptic to body-fixed
        rotation_to_body_fixed = (
            rotation_z(self._data["phi0"]).T
            @ rotation_y(np.pi / 2 - self._data["beta"]).T
            @ rotation_z(self._data["lambda"]).T
        )
        if base_frame == "J2000":
            rotation_to_body_fixed = rotation_to_body_fixed @ j2000_to_eclipj2000()

        return rotation_model.simple(
            base_frame=base_frame,
            target_frame=f"shape_model_{self._data.name}",
            initial_orientation=rotation_to_body_fixed,
            initial_time=self._data["epoch"],
            rotation_rate=2 * np.pi / self._data["period"],
        )

    def polyhedron_shape_model(
        self,
        diameter: float | None = None,
        compute_altitude_with_sign: bool = True,
        just_compute_distance_to_vertices: bool = False,
    ):
        """
        Convert the shape model into Tudat polyhedron shape model settings

        Parameters
        ----------
        diameter : float, optional
            A volume equivalent diameter can be provided, which is used when the shape model has no size information.
            It is used as a fall-back option, meaning the DAMIT sizes are always used when possible
        compute_altitude_with_sign : bool, optional
            See ...
        just_compute_distance_to_vertices : bool, optional
            See ...

        Returns
        -------

        """
        with open(self._path_to_shape_model, "r") as f:
            number_of_vertices, number_of_facets = [
                int(num) for num in f.readline().strip().split()
            ]

        vertices_and_facets = np.loadtxt(self._path_to_shape_model, skiprows=1)
        vertices = vertices_and_facets[:number_of_vertices]
        facets = vertices_and_facets[number_of_vertices:].astype(int)

        # Scale model for size
        if self._data["calibrated_size"] == 1:
            # Shape model size is in KM3, resize to m3
            scale_factor = 1e3
        elif self._data["equiv_diameter"] or diameter:
            used_diameter = (
                self._data["equiv_diameter"] if self._data["equiv_diameter"] else diameter
            )
            required_volume = (4 / 3) * np.pi * (used_diameter / 2) ** 3
            scale_factor = np.cbrt(required_volume)  # Assumes DAMIT model is unit volume
        else:
            scale_factor = 1
            warn(
                "This asteroid shape model has no size information and no diameter was given. "
                "The size of the polyhedron shape model will be unit volume"
            )

        vertices *= scale_factor

        return polyhedron(
            vertices_coordinates=vertices,
            vertices_defining_each_facet=facets - 1,  # Convert to zero based indexing
            compute_altitude_with_sign=compute_altitude_with_sign,
            just_compute_distance_to_vertices=just_compute_distance_to_vertices,
        )

    def best_fitting_ellipsoid_shape_model(self):
        pass

    def best_fitting_ellipsoid_rotational_ephemeris_settings(self):
        pass
