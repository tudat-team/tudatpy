from pathlib import Path
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch
from tudatpy.dynamics.environment_setup.shape import polyhedron
from tudatpy.dynamics.environment_setup import rotation_model
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from tudatpy.resource_installer import download_extra_file, DEFAULT_DEST, DEFAULT_CACHE
import tarfile
import sys
from functools import cache
import matplotlib.pyplot as plt

# Location of DAMIT archive files under tudat/resources
_DAMIT_RESOURCES = DEFAULT_DEST / "damit"


@cache
def _damit_files_path() -> Path:
    """Find the path to the DAMIT catalog"""
    # Lookup files in damit resources
    if _DAMIT_RESOURCES.is_dir():
        for path in _DAMIT_RESOURCES.iterdir():
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

    _DAMIT_RESOURCES.mkdir(parents=True, exist_ok=True)
    url = "https://damit.cuni.cz/projects/damit/exports/complete/latest"

    tar_path = download_extra_file(
        url=url, dest_path=_DAMIT_RESOURCES / "damit.tar.gz", cache_dir=DEFAULT_CACHE, force=False
    )

    print("Extracting files...")
    with tarfile.open(tar_path, "r:*") as tar:
        tar.extractall(_DAMIT_RESOURCES, filter="data")

    tar_path.unlink()  # Remove compressed file

    print(f"DAMIT files saved to {str(_DAMIT_RESOURCES)}")


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
    asteroid_models[["equiv_diameter", "equiv_diameter_err"]] = (
        1e3 * asteroid_models[["equiv_diameter", "equiv_diameter_err"]]
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
        diameter: float | None = None,
    ):
        """
        Initialize a DamitModel object

        Parameters
        ----------
        mpc_number : int, optional
            If given, queries a DamitModel by small-body MPC number. The model with the highest quality flag is selected.
            If there is no highest flag, selection is arbitrary
        damit_model_id : int, optional
            If given, queries a DamitModel by the Damit model ID
        diameter : float, optional
            Volume-equivalent diameter to scale the shape model to. This is used as a fall-back option, in case the
            model has no shape information available.
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

        # Load the shape model
        self._vertices, self._facets = self._load_shape_model(self._data, diameter)

    @property
    def metadata(self) -> pd.Series:
        """Returns a read-only copy of the model metadata"""
        return self._data.copy(deep=True)

    @staticmethod
    def _load_shape_model(
        model_metadata: pd.Series, diameter: float | None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Locates and loads the polyhedron shape model"""

        # Load shape model from DAMIT data
        shape_model_file = (
            _damit_files_path()
            / "files"
            / f'asteroid_{model_metadata["asteroid_id"]}'
            / f"model_{model_metadata["id"]}"
            / "shape.txt"
        )
        if not shape_model_file.exists():
            raise FileNotFoundError("Shape model file for this asteroid model could not be found")

        with open(shape_model_file, "r") as f:
            number_of_vertices, number_of_facets = [
                int(num) for num in f.readline().strip().split()
            ]

        vertices_and_facets = np.loadtxt(shape_model_file, skiprows=1)
        vertices = vertices_and_facets[:number_of_vertices]
        facets = vertices_and_facets[number_of_vertices:].astype(int) - 1  # Apply zero-indexing

        # Scale model for size
        if model_metadata["calibrated_size"] == 1:
            # Shape model size is in KM3, resize to m3
            scale_factor = 1e3
        elif pd.notna(model_metadata["equiv_diameter"]) or diameter:
            used_diameter = (
                model_metadata["equiv_diameter"]
                if pd.notna(model_metadata["equiv_diameter"])
                else diameter
            )
            required_volume = (4 / 3) * np.pi * (used_diameter / 2) ** 3
            scale_factor = np.cbrt(required_volume)  # Assumes DAMIT model is unit volume
        else:
            raise RuntimeError(
                "This asteroid shape model has no size information and no diameter was given. "
                "It cannot be converted into a shape model"
            )

        return vertices * scale_factor, facets

    def polyhedron_rotational_ephemeris_settings(self, base_frame: str) -> rotation_model.simple:
        """
        Create rotational ephemeris settings of the polyhedron shape model

        Parameters
        ----------
        base_frame : str, optional
            Base frame used in the environment. Options: "J2000", "ECLIPJ2000"

        Returns
        -------
        RotationModelSettings
            Simple rotation model settings for this DAMIT shape model
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

    def polyhedron_shape_model(self, *args, **kwargs) -> polyhedron:
        """
        Convert the shape model into Tudat polyhedron shape model settings

        Parameters
        ----------
        args, kwargs
            Remaining inputs are passed to shape.polyhedron()

        Returns
        -------

        """
        return polyhedron(
            vertices_coordinates=self._vertices,
            vertices_defining_each_facet=self._facets,
            *args,
            **kwargs,
        )

    def plot_shape_model(
        self,
        show: bool = True,
    ) -> plt.Axes:
        """
        Plots the polyhedron using Matplotlib's 3D plotting interface

        Parameters
        ----------
        show : bool, optional
            If true, calls ``matplotlib.pyplot.show``

        Returns
        -------
        plt.Axes
            Matplotlib Axes object
        """
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        ax.add_collection3d(
            Poly3DCollection(
                self._vertices[self._facets],
                facecolors="lightgrey",
                shade=True,
                label="Shape model",
            )
        )

        ax.set_aspect("equal")
        ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
        ax.set_title(f"Shape model MPC {self._data['mpc_number']} (model ID {self._data['id']})")
        ax.legend()

        if show:
            plt.show()

        return ax

    def best_fitting_ellipsoid_shape_model(self):
        pass

    def best_fitting_ellipsoid_rotational_ephemeris_settings(self):
        pass
