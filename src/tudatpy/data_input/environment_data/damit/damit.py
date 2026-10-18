from pathlib import Path
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from tudatpy.dynamics.environment_setup.rotation_model import RotationModelSettings
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch
from tudatpy.dynamics.environment_setup.shape import (
    polyhedron,
    triaxial_ellipsoid,
    PolyhedronBodyShapeSettings,
    TriAxialEllipsoidBodyShapeSettings,
)
from tudatpy.dynamics.environment_setup import rotation_model
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from tudatpy.resource_installer import download_extra_file, DEFAULT_DEST, DEFAULT_CACHE
import tarfile
import shutil
import matplotlib.pyplot as plt

try:
    import trimesh

    TRIMESH_AVAILABLE = True
except ModuleNotFoundError:
    TRIMESH_AVAILABLE = False

# File locations
# - Full archive download stored under .tudat / resource / damit / archive
# - Individual models and asteroid tables saved under .tudat / resource / damit / models
_DAMIT_DIR = DEFAULT_DEST / "damit"
_DAMIT_DIR.mkdir(parents=True, exist_ok=True)
_DAMIT_MODELS_DIR = DEFAULT_DEST / "damit" / "models"
_DAMIT_TABLES_DIR = _DAMIT_MODELS_DIR / "tables"
_DAMIT_ARCHIVE_DIR = DEFAULT_DEST / "damit" / "archive"
_ARCHIVE_AVAILABLE = (_DAMIT_ARCHIVE_DIR / "tables").is_dir() and (
    _DAMIT_ARCHIVE_DIR / "files"
).is_dir()


def _download_tables(force: bool = False) -> None:
    """Download DAMIT asteroid and asteroid model tables"""
    # Check if they exist
    asteroid_models_csv = _DAMIT_TABLES_DIR / "asteroid_models.csv"
    asteroids_csv = _DAMIT_TABLES_DIR / "asteroids.csv"
    tables_exist = asteroids_csv.exists() and asteroid_models_csv.exists()

    # Download if tables don't exist yet, or if they should be overwritten
    if not tables_exist or force:
        print("Downloading DAMIT tables...")
        _DAMIT_TABLES_DIR.mkdir(parents=True, exist_ok=True)
        download_extra_file(
            url="https://damit.cuni.cz/projects/damit/exports/table/asteroids",
            dest_path=_DAMIT_TABLES_DIR / "asteroids.csv",
            cache_dir=DEFAULT_CACHE,
            force=True,
        )
        download_extra_file(
            url="https://damit.cuni.cz/projects/damit/exports/table/asteroid_models",
            dest_path=_DAMIT_TABLES_DIR / "asteroid_models.csv",
            cache_dir=DEFAULT_CACHE,
            force=True,
        )
        print("DAMIT tables saved.")


def _download_model(model_id: int, force: bool = False) -> None:
    """Downloads an individual model"""
    # Check if model exists on disk
    model_file = _DAMIT_MODELS_DIR / str(model_id) / "obj.txt"

    # Download if model does not exist yet or if it should be overwritten
    if not model_file.exists() or force:
        _DAMIT_MODELS_DIR.mkdir(parents=True, exist_ok=True)
        base_url = "https://damit.cuni.cz/projects/damit/generated_files/open/AsteroidModel"
        url_shape_obj = f"{base_url}/{model_id}/shape.obj"

        print(f"Downloading model with model ID {model_id}...")

        download_extra_file(
            url=url_shape_obj,
            dest_path=model_file,
            cache_dir=DEFAULT_CACHE,
            force=True,
        )
        print("Model saved.")


def download_damit_archive(force: bool = False) -> None:
    """
    This function downloads the latest complete DAMIT archive containing shape models for roughly
    10,000 asteroids. The complete archive takes up roughly 4 GB of disk space.

    Parameters
    ----------
    force : bool
        Overwrites existing archive if it exists
    """
    global _ARCHIVE_AVAILABLE
    if _ARCHIVE_AVAILABLE and not force:
        return

    while True:
        user_query = input(
            "The DAMIT archive (~4 GB) will be downloaded from the internet. Proceed? (y/n) "
        )
        if user_query.lower() == "n":
            return
        if user_query.lower() == "y":
            break

    # Download archive compressed tar.gz file
    url = "https://damit.cuni.cz/projects/damit/exports/complete/latest"
    tar_path = download_extra_file(
        url=url, dest_path=_DAMIT_DIR / "damit_latest.tar.gz", cache_dir=DEFAULT_CACHE, force=True
    )

    # Clean up existing files in damit/archive
    if _DAMIT_ARCHIVE_DIR.exists():
        shutil.rmtree(_DAMIT_ARCHIVE_DIR)

    print("Extracting files...")
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(_DAMIT_ARCHIVE_DIR, filter="data")

    # Move contents to archive/
    top_folder = next(_DAMIT_ARCHIVE_DIR.glob("damit-*"))
    for folder in ("files", "tables"):
        shutil.move(top_folder / folder, _DAMIT_ARCHIVE_DIR / folder)

    top_folder.rmdir()  # Remove empty top folder
    tar_path.unlink()  # Remove compressed gz file

    print(f"DAMIT files saved to {str(_DAMIT_ARCHIVE_DIR)}")

    _ARCHIVE_AVAILABLE = True


def _load_table() -> pd.DataFrame:
    """Loads all asteroid shape models into a pandas dataframe"""
    if _ARCHIVE_AVAILABLE:
        table_dir = _DAMIT_ARCHIVE_DIR / "tables"
    else:
        _download_tables()
        table_dir = _DAMIT_TABLES_DIR

    # Asteroid data table containing primarily the DAMIT ID's of the asteroids
    asteroid_data = pd.read_csv(table_dir / "asteroids.csv", dtype={"number": "Int64"}).rename(
        columns={"number": "mpc_number", "id": "asteroid_id"}
    )

    # Model data containing actual models
    asteroid_models = pd.read_csv(
        table_dir / "asteroid_models.csv",
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
    Returns a table with asteroid models for the queried MPC number. An asteroid may have one or more models,
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

    def __init__(
        self,
        mpc_number: int | None = None,
        damit_model_id: int | None = None,
        diameter: float | None = None,
    ):
        """
        Initialize a DamitModel object
        """
        # Query by mpc number and choose best model
        if mpc_number and not damit_model_id:
            best_model = damit_model_query(
                mpc_number,
            )
            if best_model.empty:
                raise LookupError(f"No shape models found for {mpc_number}")
            self._data = best_model.sort_values("quality_flag", ascending=False).iloc[0]

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
        # Locate model file
        model_id = model_metadata["id"]
        if _ARCHIVE_AVAILABLE:
            shape_model_file = (
                _DAMIT_ARCHIVE_DIR
                / "files"
                / f'asteroid_{model_metadata["asteroid_id"]}'
                / f"model_{model_id}"
                / "obj.txt"
            )
        else:
            _download_model(model_id)
            shape_model_file = _DAMIT_MODELS_DIR / str(model_id) / "obj.txt"

        if not shape_model_file.exists():
            raise FileNotFoundError("Shape model file for this asteroid model could not be found")

        # Load OBJ file: 'v x y z' lines are vertices, 'f i j k' lines are facets
        vertices_and_facets = pd.read_csv(
            shape_model_file, sep=r"\s+", names=["v/f", "num1", "num2", "num3"]
        )
        vertices = vertices_and_facets.loc[
            vertices_and_facets["v/f"] == "v", ["num1", "num2", "num3"]
        ].to_numpy()
        facets = (
            vertices_and_facets.loc[
                vertices_and_facets["v/f"] == "f", ["num1", "num2", "num3"]
            ].to_numpy(dtype=int)
            - 1
        )  # Apply zero-indexing

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

    def _approximate_best_fitting_ellipsoid(self):
        """Approximate the best fitting ellipsoid from the moment of inertia"""
        mesh = trimesh.Trimesh(vertices=self._vertices, faces=self._facets)
        # Principal axes of the equatorial block only
        # -> it is assumed that the shape model z-axis is a principal axis itself to preserve the correct
        # axis of rotation.

        inertia_matrix = mesh.moment_inertia
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_matrix[:2, :2])
        ixx, iyy, izz = eigenvalues[0], eigenvalues[1], inertia_matrix[2, 2]

        # Convert moments of inertia to ellipsoid axes
        mass = mesh.volume  # Model has unit density
        a = np.sqrt((5.0 / (2.0 * mass)) * (iyy + izz - ixx))
        b = np.sqrt((5.0 / (2.0 * mass)) * (ixx + izz - iyy))
        c = np.sqrt((5.0 / (2.0 * mass)) * (ixx + iyy - izz))

        # Principal axis offset required to account for the difference in the ellipsoid x-axis and DAMIT model x-axis
        principal_axis_offset = (
            np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]) + np.pi / 2
        ) % np.pi - np.pi / 2

        return [a, b, c], principal_axis_offset

    def _get_rotational_ephemeris_settings(
        self, base_frame: str, extra_rotation: float = 0.0
    ) -> RotationModelSettings:
        """Get the Simple rotation model settings, with an optional extra offset around the z-axis"""
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
            rotation_z(self._data["phi0"] + extra_rotation).T
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

    def polyhedron_rotational_ephemeris_settings(self, base_frame: str) -> RotationModelSettings:
        """
        Create rotational ephemeris settings of the polyhedron shape model

        Parameters
        ----------
        base_frame : str
            Base frame used in the environment. Options: "J2000", "ECLIPJ2000"

        Returns
        -------
        RotationModelSettings
            Simple rotation model settings for this DAMIT shape model
        """
        return self._get_rotational_ephemeris_settings(base_frame)

    def polyhedron_shape_model(self, *args, **kwargs) -> PolyhedronBodyShapeSettings:
        """
        Convert the shape model into Tudat polyhedron shape model settings

        Parameters
        ----------
        args, kwargs
            Remaining inputs are passed to shape.polyhedron()

        Returns
        -------
        PolyhedronBodyShapeSettings
            Shape settings objects for the DAMIT model
        """
        return polyhedron(
            vertices_coordinates=self._vertices,
            vertices_defining_each_facet=self._facets,
            *args,
            **kwargs,
        )

    def plot_shape_model(
        self,
        hide_axes: bool = True,
        show: bool = True,
    ) -> plt.Axes:
        """
        Plots the polyhedron using Matplotlib's 3D plotting interface

        Parameters
        ----------
        hide_axes: bool, optional
            Hides the axes and gridlines
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
            )
        )

        ax.set_aspect("equal")
        ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
        ax.set_title(f"Shape model MPC {self._data['mpc_number']} (model ID {self._data['id']})")

        if hide_axes:
            ax.set_axis_off()

        if show:
            plt.show()

        return ax

    def best_fitting_ellipsoid_shape_model(self) -> TriAxialEllipsoidBodyShapeSettings:
        """
        Create the best fitting triaxial ellipsoid shape model based on the polyhedron DAMIT shape model.

        Creates a triaxial ellipsoid shape model for an ellipsoid with equivalent moment of inertia to the DAMIT
        polyhedron shape model.

        NOTE: this functionality requires the optional dependency ``trimesh`` to be installed, at this time.
        This can be installed into the conda environment using ``conda install trimesh``.

        Returns
        -------
        TriaxialEllipsoidBodyShapeSettings
            Shape settings for the best-fitting triaxial ellipsoid shape model
        """
        if not TRIMESH_AVAILABLE:
            raise ModuleNotFoundError(
                "`trimesh` is required to approximate the best-fitting ellipsoids at this time. "
                "Install it with `conda install trimesh` in the conda environment"
            )
        semi_axes, _ = self._approximate_best_fitting_ellipsoid()
        return triaxial_ellipsoid(*semi_axes)

    def best_fitting_ellipsoid_rotational_ephemeris_settings(
        self, base_frame: str
    ) -> RotationModelSettings:
        """
        Create rotational ephemeris settings of the best fitting ellipsoid shape model

        Creates the rotational ephemeris settings of the best fitting ellipsoid with equivalent moment of inertia
        to the DAMIT polyhedron shape model.

        NOTE 1: this rotational model is *not* the same as computed by:
        :meth:~`tudatpy.data_input.environment_data.damit.DamitModel.polyhedron_rotational_ephemeris_settings`,
        because an extra offset to the initial phase of rotation is necessary for the ellipsoid model.

        NOTE 2: this functionality requires the optional dependency ``trimesh`` to be installed, at this time.
        This can be installed into the conda environment using ``conda install trimesh``.

        Parameters
        ----------
        base_frame : str
            Base frame used in the environment. Options: "J2000", "ECLIPJ2000

        Returns
        -------
        RotationModelSettings
            Simple rotation model settings for the best fitting ellipsoid
        """
        if not TRIMESH_AVAILABLE:
            raise ModuleNotFoundError(
                "`trimesh` is required to approximate the best-fitting ellipsoids at this time. "
                "Install it with `conda install trimesh` in the conda environment"
            )
        _, principal_axis_offset = self._approximate_best_fitting_ellipsoid()
        return self._get_rotational_ephemeris_settings(base_frame, principal_axis_offset)
