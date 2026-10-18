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
from functools import cache

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
    """Download DAMIT asteroid and asteroid model tables, containing model metadata such as
    size information, asteroid designation, and rotational parameters"""
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


def _download_model(model_id: int, force: bool = False) -> Path:
    """Downloads an individual model shape file containing vertex coordinates and facets"""
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

    return model_file


def download_damit_archive(tables_only: bool = True, force: bool = False) -> None:
    """
    This function downloads the latest complete DAMIT archive containing shape models for roughly
    10,000 asteroids. The complete archive takes up roughly 4 GB of disk space.

    Parameters
    ----------
    tables_only : bool, optional
        If true, only downloads the model tables, containing the metadata of all the models. If set to false,
        it downloads the complete archive, including shape models, images, and light-curve data (±4 GB).
    force : bool, optional
        Overwrites any existing files
    """
    if tables_only:
        _download_tables(force=force)
        return

    global _ARCHIVE_AVAILABLE
    if _ARCHIVE_AVAILABLE and not force:
        print("DAMIT archive already downloaded.")
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


@cache
def _load_model_table() -> pd.DataFrame:
    """Loads all asteroid shape models into a pandas dataframe"""
    # Find location of asteroid models files
    if _ARCHIVE_AVAILABLE:
        table_dir = _DAMIT_ARCHIVE_DIR / "tables"
    else:
        #  Will download the tables when run for the first time
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

    return asteroid_models


def list_damit_models(
    mpc_number: int | None = None,
) -> pd.DataFrame:
    """
    Returns a table with asteroid models for the queried MPC number. An asteroid may have one or more models,
    or none, in which case the table is returned empty. Asteroid model tables will be downloaded if run for the first time

    Parameters
    ----------
    mpc_number : int, optional
        MPC number to query. If None, all models are returned.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the models that match the query
    """
    table = _load_model_table().copy()

    if mpc_number is None:
        return table

    query_result = table[table["mpc_number"] == mpc_number]
    return query_result


class DamitModel:
    """
    Interface with DAMIT database for asteroid shape models and rotational parameters
    """

    def __init__(
        self,
        name: str,
        rotational_parameters: dict,
        vertices: np.ndarray,
        facets: np.ndarray,
    ):
        """
        Initialize a DamitModel object

        Parameters
        ----------
        name : str
            Name of the model, used for the reference frame name
        rotational_parameters : dict
            Dictionary containing rotational parameters
        vertices : np.ndarray
            Vertex coordinates
        facets : np.ndarray
            Facet coordinates
        """
        # DamitModel attributes
        self._name = name
        self._rotational_parameters = rotational_parameters
        self._vertices = vertices
        self._facets = facets

    @classmethod
    def query(
        cls,
        mpc_number: int | None = None,
        damit_model_id: int | None = None,
        diameter: float | None = None,
    ):
        """
        Queries the DAMIT database by MPC number or model ID. Models are downloaded from the internet and stored on disk
        for faster subsequent retrieval. Alternatively, the function download_damit_archive can be used to download all
        available shape models instead.

        Parameters
        ----------
        mpc_number : int, optional
            If given, queries a DamitModel by small-body MPC number. The model with the highest quality flag is selected.
            If there is no highest flag, selection is arbitrary
        damit_model_id : int, optional
            If given, queries a DamitModel by the Damit model ID
        diameter : float, optional
            Volume-equivalent diameter to scale the shape model to. This is used as a fall-back option, in case the
            model has no shape information available."""

        # Get the model metadata (rotational parameters, size info etc.)
        # Query by mpc number and choose best model
        if mpc_number and not damit_model_id:
            available_models = list_damit_models(
                mpc_number,
            )
            if available_models.empty:
                raise LookupError(f"No shape models found for {mpc_number}")
            model_metadata = available_models.sort_values("quality_flag", ascending=False).iloc[0]

        # Retrieve by model ID
        elif damit_model_id and not mpc_number:
            try:
                model_metadata = list_damit_models().loc[damit_model_id]
            except KeyError as e:
                raise LookupError(f"No shape model found with ID {damit_model_id}") from e

        else:
            model_metadata = None
            raise ValueError("Must specify either a valid MPC number or Damit Model ID")

        # Get size info
        if model_metadata["calibrated_size"] == 1:
            diameter = None
        elif pd.notna(model_metadata["equiv_diameter"]):
            diameter = model_metadata["equiv_diameter"] * 1e3  # KM -> M
        elif diameter is None:
            raise ValueError(
                "model has no size information and no diameter was given. "
                "A shape model cannot be created."
            )

        # Load model
        asteroid_id, model_id = model_metadata["asteroid_id"], model_metadata["id"]
        model_file = cls._get_shape_model_file(model_id, asteroid_id)
        vertices, facets = cls._load_shape_model(model_file, diameter)

        # Rotation parameters
        rotational_parameters = dict(
            pole_longitude=np.deg2rad(model_metadata["lambda"]),
            pole_latitude=np.deg2rad(model_metadata["beta"]),
            period=model_metadata["period"] * 3600,
            phi0=np.deg2rad(model_metadata["phi0"]),
            epoch=julian_day_to_seconds_since_epoch(model_metadata["epoch"]),
        )

        return cls(
            model_id,
            rotational_parameters,
            vertices,
            facets,
        )

    @classmethod
    def from_file(
        cls,
        spin_file: str | Path,
        shape_file: str | Path,
        model_name: str,
        diameter: float | None = None,
    ):
        """
        Loads a DAMIT model from files

        Parameters
        ----------
        spin_file : str | Path
            Path to spin.txt file
        shape_file : str | Path
            Path to shape.obj(/obj.txt) file
        model_name : str
            Name of the model, used for the reference frame name
        diameter : float | None
            Diameter in meter. If given, scales the model to this size.

        Returns
        -------
        DamitModel
        """
        with open(spin_file, "r") as f:
            pole_longitude, pole_latitude, period = [
                float(num) for num in f.readline().strip().split()
            ]
            epoch, phi0 = [float(num) for num in f.readline().strip().split()]

        rotational_parameters = dict(
            pole_longitude=np.deg2rad(pole_longitude),
            pole_latitude=np.deg2rad(pole_latitude),
            period=period * 3600,
            phi0=np.deg2rad(phi0),
            epoch=julian_day_to_seconds_since_epoch(epoch),
        )
        vertices, facets = cls._load_shape_model(shape_file, diameter)

        return cls(
            model_name,
            rotational_parameters,
            vertices,
            facets,
        )

    @staticmethod
    def _get_shape_model_file(model_id: int, asteroid_id: int) -> Path:
        """Find shape model files on disk and download if not present"""
        if _ARCHIVE_AVAILABLE:
            shape_model_file = (
                _DAMIT_ARCHIVE_DIR
                / "files"
                / f"asteroid_{asteroid_id}"
                / f"model_{model_id}"
                / "obj.txt"
            )
        else:
            shape_model_file = _download_model(model_id)

        if not shape_model_file.exists():
            raise FileNotFoundError("Shape model file for this asteroid model could not be found")

        return shape_model_file

    @staticmethod
    def _load_shape_model(
        shape_model_file, diameter: float | None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Loads the shape model file and optionally scales by diameter"""
        # Check format
        with open(shape_model_file, "r") as f:
            first_line = f.readline().strip().split()
            if not (len(first_line) == 4 and first_line[0] == "v"):
                raise ValueError(f"shape file {str(shape_model_file)} must be of the obj format")

        # Load OBJ file: 'v x y z' lines are vertices, 'f i j k' lines are facets
        vertices_and_facets = pd.read_csv(
            shape_model_file, sep=r"\s+", names=["v/f", "num1", "num2", "num3"]
        )
        vertices = vertices_and_facets.loc[
            vertices_and_facets["v/f"] == "v", ["num1", "num2", "num3"]
        ].to_numpy()
        facets = vertices_and_facets.loc[
            vertices_and_facets["v/f"] == "f", ["num1", "num2", "num3"]
        ].to_numpy(dtype=int)

        # Scale model for size
        if diameter is None:
            # Assume size is in km^3
            scale_factor = 1e3

        else:
            required_volume = (4 / 3) * np.pi * (diameter / 2) ** 3
            scale_factor = np.cbrt(required_volume)  # Assumes DAMIT model is unit volume

        return vertices * scale_factor, facets - 1  # zero-indexing

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
            rotation_z(self._rotational_parameters["phi0"] + extra_rotation).T
            @ rotation_y(np.pi / 2 - self._rotational_parameters["pole_latitude"]).T
            @ rotation_z(self._rotational_parameters["pole_longitude"]).T
        )
        if base_frame == "J2000":
            rotation_to_body_fixed = rotation_to_body_fixed @ j2000_to_eclipj2000()

        return rotation_model.simple(
            base_frame=base_frame,
            target_frame=f"shape_model_{self._name}",
            initial_orientation=rotation_to_body_fixed,
            initial_time=self._rotational_parameters["epoch"],
            rotation_rate=2 * np.pi / self._rotational_parameters["period"],
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
            self._vertices,
            self._facets,
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
