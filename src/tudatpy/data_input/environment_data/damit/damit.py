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
import matplotlib.pyplot as plt
from tudatpy.astro.polyhedron_utilities import volume, inertia_tensor_from_density
from warnings import warn

# File locations
# - Shape models .tudat / resource / damit / models
# - Damit tables .tudat / resource / damit / tables
_DAMIT_DIR = DEFAULT_DEST / "damit"
_DAMIT_MODELS_DIR = _DAMIT_DIR / "models"

_DAMIT_TABLES_DIR = _DAMIT_DIR / "tables"
_DAMIT_ASTEROID_MODELS = _DAMIT_TABLES_DIR / "asteroid_models.csv"
_DAMIT_ASTEROIDS = _DAMIT_TABLES_DIR / "asteroids.csv"


def _download_tables(force: bool = False) -> None:
    """Download DAMIT asteroid and asteroid model tables, necessary for:
    - Querying models by MPC number
    - Getting volume-equivalent diameter
    - Rotational parameters
    """
    # Check if tables already exist
    tables_exist = _DAMIT_ASTEROID_MODELS.is_file() and _DAMIT_ASTEROIDS.is_file()

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
        print(f"DAMIT tables saved at {str(_DAMIT_TABLES_DIR)}.\n")


def _download_model(model_id: int, force: bool = False) -> tuple[Path, Path]:
    """Downloads an individual asteroid model
    - shape.obj: contains vertices and facets
    - spin.txt: contains rotational parameters"""
    # Check if model exists on disk
    model_dir = _DAMIT_MODELS_DIR / str(model_id)
    obj_file = model_dir / "shape.obj"
    spin_file = model_dir / "spin.txt"
    model_exists = obj_file.is_file() and spin_file.is_file()

    # Download if model does not exist yet or if it should be overwritten
    if not model_exists or force:
        print(f"Downloading model with model ID {model_id}...")
        _DAMIT_MODELS_DIR.mkdir(parents=True, exist_ok=True)
        base_url = "https://damit.cuni.cz/projects/damit/generated_files/open/AsteroidModel"

        for file_path in (obj_file, spin_file):
            download_extra_file(
                url=f"{base_url}/{model_id}/{file_path.name}",
                dest_path=file_path,
                cache_dir=DEFAULT_CACHE,
                force=True,
            )

        print(f"Model saved at {str(model_dir)}\n")

    return obj_file, spin_file


def _load_model_table(force: bool = False) -> pd.DataFrame:
    """Loads all asteroid shape models into a pandas dataframe"""
    _download_tables(force)

    # Asteroid data table containing primarily the DAMIT ID's of the asteroids
    asteroid_data = pd.read_csv(_DAMIT_ASTEROIDS, dtype={"number": "Int64"}).rename(
        columns={"number": "mpc_number", "id": "asteroid_id"}
    )

    # Model data containing actual models
    asteroid_models = pd.read_csv(
        _DAMIT_ASTEROID_MODELS,
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
    force_refresh: bool = False,
) -> pd.DataFrame:
    """
    Returns a table with asteroid models for the queried MPC number. An asteroid may have one or more models,
    or none, in which case the table is returned empty. Asteroid model tables will be downloaded if run for the first time

    Parameters
    ----------
    mpc_number : int, optional
        MPC number to query. If None, all models are returned.
    force_refresh : bool, optional
        If true, redownloads and overwrites existing tables to the latest versions.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the models that match the query
    """
    table = _load_model_table(force_refresh)

    if mpc_number is None:
        return table

    query_result = table[table["mpc_number"] == mpc_number]
    return query_result


class DamitModel:
    """
    Initialize a DamitModel object from the shape and rotational parameters.

    Normally, use either :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query`
    or :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.from_file` to initialize a model.

    Parameters
    ----------
    name : str
        Name of the model, used for the reference frame name
    vertices : np.ndarray
        Vertex coordinates in meters
    facets : np.ndarray
        Facet coordinates, zero-indexed
    pole_longitude : float
        Longitude in ecliptic J2000 coordinates (rad)
    pole_latitude : float
        Latitude in ecliptic J2000 coordinates (rad)
    period: float
        Period in seconds. Note: changes to rotation rate due to YORP effect are ignored
    initial_rotation_angle: float
        Initial rotation angle in (rad)
    reference_epoch: float
        Rotation reference epoch in seconds since J2000
    """

    def __init__(
        self,
        name: str,
        vertices: np.ndarray,
        facets: np.ndarray,
        pole_longitude: float,
        pole_latitude: float,
        period: float,
        initial_rotation_angle: float,
        reference_epoch: float,
    ):
        # Input check
        if not vertices.shape[1] == 3 or not facets.shape[1] == 3:
            raise ValueError("Vertices and facets must have three columns")
        if facets.max() > len(vertices) - 1 or facets.min() < 0:
            raise ValueError("Facets indices out of range")

        if period <= 0:
            raise ValueError("Period must be positive")

        # DamitModel attributes
        self._name = name
        self._vertices = vertices
        self._facets = facets
        self._pole_longitude = pole_longitude
        self._pole_latitude = pole_latitude
        self._period = period
        self._initial_rotation_angle = initial_rotation_angle
        self._reference_epoch = reference_epoch

    @property
    def vertices_and_facets(self) -> tuple[np.ndarray, np.ndarray]:
        """Returns a copy of the vertices and facets of this model"""
        return self._vertices.copy(), self._facets.copy()

    @classmethod
    def query(
        cls,
        *,
        mpc_number: int | None = None,
        damit_model_id: int | None = None,
        diameter: float | None = None,
        force_refresh: bool = False,
    ):
        """
        Queries the DAMIT database by MPC number or model ID. Models are downloaded from the internet and stored on disk
        for faster subsequent retrieval.

        Parameters
        ----------
        mpc_number : int, optional
            If given, queries a DamitModel by small-body MPC number. The model with the highest quality flag is selected.
            If there is no highest flag, selection is arbitrary
        damit_model_id : int, optional
            If given, queries a DamitModel by the Damit model ID
        diameter : float, optional
            Volume-equivalent diameter to scale the shape model to. If not given, it uses the size that is available
            in the DAMIT catalog, which may be none in which case the size is garbage
        force_refresh : bool
            If true, re-downloads and overwrites models to the latest versions

        Returns
        -------
        DamitModel
        """
        # Get the model metadata (rotational parameters, size info etc.)
        # Query by mpc number and choose best model
        if mpc_number and not damit_model_id:
            available_models = list_damit_models(
                mpc_number,
                force_refresh=force_refresh,
            )
            if available_models.empty:
                raise LookupError(f"No shape models found for {mpc_number}")
            model_metadata = available_models.sort_values("quality_flag", ascending=False).iloc[0]

        # Retrieve by model ID
        elif damit_model_id and not mpc_number:
            try:
                model_metadata = list_damit_models(force_refresh=force_refresh).loc[damit_model_id]
            except KeyError as e:
                raise LookupError(f"No shape model found with ID {damit_model_id}") from e

        else:
            model_metadata = None
            raise ValueError(
                "Must specify either a valid MPC number or Damit Model ID (and not both)"
            )

        # Get diameter from metadata if size is not calibrated and diameter is not given
        if diameter is None and model_metadata["calibrated_size"] != 1:
            if pd.notna(model_metadata["equiv_diameter"]):
                diameter = model_metadata["equiv_diameter"] * 1e3  # km -> m
            else:
                warn(
                    "Queried model has no size information and no diameter was given. "
                    "The asteroid size will be non-physical"
                )

        # Load model
        # Download shape models
        model_id = model_metadata["id"]
        shape_file, spin_file = _download_model(model_id, force_refresh)

        return cls.from_file(spin_file, shape_file, str(model_id), diameter)

    @classmethod
    def from_file(
        cls,
        spin_file: str | Path,
        shape_file: str | Path,
        model_name: str,
        diameter: float | None = None,
    ):
        """
        Loads a DAMIT model from files. It is assumed the files are in the exact same format as those downloaded from
        the DAMIT database.

        Note: if no diameter is given, the model uses whatever the volume of the shape is already in (scaled for the
        correct units), which may be non-physical. To be safe, it is recommended to pass a diameter.

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

        vertices, facets = cls._load_shape_model(shape_file, diameter)

        return cls(
            name=model_name,
            vertices=vertices,
            facets=facets,
            pole_longitude=np.deg2rad(pole_longitude),
            pole_latitude=np.deg2rad(pole_latitude),
            period=period * 3600,
            initial_rotation_angle=np.deg2rad(phi0),
            reference_epoch=julian_day_to_seconds_since_epoch(epoch),
        )

    @staticmethod
    def _load_shape_model(
        shape_model_file,
        diameter: float | None,
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
        facets = (
            vertices_and_facets.loc[
                vertices_and_facets["v/f"] == "f", ["num1", "num2", "num3"]
            ].to_numpy(dtype=int)
            - 1
        )  # Apply zero-indexing

        # Scale the model size
        if diameter is not None:
            current_model_volume = volume(vertices, facets)
            if current_model_volume <= 0:
                raise RuntimeError("Broken shape model: shape model has negative volume")
            required_volume = (4 / 3) * np.pi * (diameter / 2) ** 3
            vertices *= np.cbrt(required_volume / current_model_volume)
        else:
            vertices *= 1e3  # converts km^3 to m^3 for models that are not scaled

        return vertices, facets

    def _approximate_best_fitting_ellipsoid(self) -> tuple[list, float]:
        """Approximate the best fitting ellipsoid from the moment of inertia"""
        # Principal axes of the equatorial block only
        # -> it is assumed that the shape model z-axis is a principal axis itself to preserve the correct
        # axis of rotation.

        inertia_matrix = inertia_tensor_from_density(self._vertices, self._facets, 1)
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_matrix[:2, :2])
        ixx, iyy, izz = eigenvalues[0], eigenvalues[1], inertia_matrix[2, 2]

        # Convert moments of inertia to ellipsoid axes
        mass = volume(self._vertices, self._facets)  # Model has unit density
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
            rotation_z(self._initial_rotation_angle + extra_rotation).T
            @ rotation_y(np.pi / 2 - self._pole_latitude).T
            @ rotation_z(self._pole_longitude).T
        )
        if base_frame == "J2000":
            rotation_to_body_fixed = rotation_to_body_fixed @ j2000_to_eclipj2000()

        return rotation_model.simple(
            base_frame=base_frame,
            target_frame=f"shape_model_{self._name}",
            initial_orientation=rotation_to_body_fixed,
            initial_time=self._reference_epoch,
            rotation_rate=2 * np.pi / self._period,
        )

    def polyhedron_rotational_ephemeris_settings(self, base_frame: str) -> RotationModelSettings:
        """
        Create rotational ephemeris settings of the polyhedron shape model. Rotation is assumed to be around a fixed
        axis in inertial space with a fixed rate of rotation

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
            Remaining inputs are passed to :func:`~tudatpy.dynamics.environment_setup.shape.polyhedron`

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

        Returns
        -------
        TriaxialEllipsoidBodyShapeSettings
            Shape settings for the best-fitting triaxial ellipsoid shape model
        """
        semi_axes, _ = self._approximate_best_fitting_ellipsoid()
        return triaxial_ellipsoid(*semi_axes)

    def best_fitting_ellipsoid_rotational_ephemeris_settings(
        self, base_frame: str
    ) -> RotationModelSettings:
        """
        Create rotational ephemeris settings of the best fitting ellipsoid shape model.

        Creates the rotational ephemeris settings of the best fitting ellipsoid with equivalent moment of inertia
        to the DAMIT polyhedron shape model. Rotation is assumed to be around a fixed axis in inertial space with a
        fixed rate of rotation

        NOTE: this rotational model is **not** the same as computed by:
        :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.polyhedron_rotational_ephemeris_settings`,
        because an extra offset to the initial phase of rotation is necessary for the ellipsoid model.

        Parameters
        ----------
        base_frame : str
            Base frame used in the environment. Options: "J2000", "ECLIPJ2000

        Returns
        -------
        RotationModelSettings
            Simple rotation model settings for the best fitting ellipsoid
        """
        _, principal_axis_offset = self._approximate_best_fitting_ellipsoid()
        return self._get_rotational_ephemeris_settings(base_frame, principal_axis_offset)
