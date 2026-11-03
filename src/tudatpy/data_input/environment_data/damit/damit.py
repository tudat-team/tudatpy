from pathlib import Path
import numpy as np
import pandas as pd
from tudatpy.astro.time_representation import julian_day_to_seconds_since_epoch
from tudatpy.astro.element_conversion import j2000_to_eclipj2000
from tudatpy.astro.frame_conversion import inertial_to_body_fixed_rotation_matrix
from tudatpy.resource_installer import download_extra_file, DEFAULT_DEST, DEFAULT_CACHE
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
    for table in (_DAMIT_ASTEROIDS, _DAMIT_ASTEROID_MODELS):
        download_extra_file(
            url=f"https://damit.cuni.cz/projects/damit/exports/table/{table.stem}",
            dest_path=table,
            cache_dir=DEFAULT_CACHE,
            force=force,
        )


def _download_model(model_id: int, force: bool = False) -> tuple[Path, Path]:
    """Downloads an individual asteroid model
    - shape.obj: contains vertices and facets
    - spin.txt: contains rotational parameters"""
    model_dir = _DAMIT_MODELS_DIR / str(model_id)
    obj_file = model_dir / "shape.obj"
    spin_file = model_dir / "spin.txt"
    base_url = "https://damit.cuni.cz/projects/damit/generated_files/open/AsteroidModel"

    for file_path in (obj_file, spin_file):
        download_extra_file(
            url=f"{base_url}/{model_id}/{file_path.name}",
            dest_path=file_path,
            cache_dir=DEFAULT_CACHE,
            force=force,
        )

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

    To create shape and rotation model settings from a DamitModel, use
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_polyhedron`,
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_best_fitting_ellipsoid`,
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit_polyhedron` and
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit_best_fitting_ellipsoid`.

    Parameters
    ----------
    vertices : np.ndarray
        Vertex coordinates in meters
    facets : np.ndarray
        Facet coordinates, zero-indexed
    pole_longitude : float
        Longitude in ecliptic J2000 coordinates (rad)
    pole_latitude : float
        Latitude in ecliptic J2000 coordinates (rad)
    period: float
        Period in seconds. Note: changes to rotation rate due to YORP effect are ignored. For models obtained from
        :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query`, a warning is given once per body when
        the rotation of such a model is used.
    initial_rotation_angle: float
        Initial rotation angle in (rad)
    reference_epoch: float
        Rotation reference epoch in seconds since J2000
    """

    def __init__(
        self,
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
        self._vertices = vertices
        self._facets = facets
        self._pole_longitude = pole_longitude
        self._pole_latitude = pole_latitude
        self._period = period
        self._initial_rotation_angle = initial_rotation_angle
        self._reference_epoch = reference_epoch
        self._rotation_warning = None
        self._number_of_equally_good_models = None

    def _warn_if_rotation_incomplete(self) -> None:
        """Warns if rotation is used of a model with ignored rotational data. Python reports a warning with the same
        text only once, and the text identifies the body, so this is reported once per body."""
        if self._rotation_warning is not None:
            warn(self._rotation_warning)

    @property
    def vertices(self) -> np.ndarray:
        """**read-only**

        Copy of the vertex coordinates (m)
        """
        return self._vertices.copy()

    @property
    def facets(self) -> np.ndarray:
        """**read-only**

        Copy of the facet vertex indices, zero-indexed
        """
        return self._facets.copy()

    @property
    def pole_longitude(self) -> float:
        """**read-only**

        Pole longitude in ecliptic J2000 coordinates (rad)
        """
        self._warn_if_rotation_incomplete()
        return self._pole_longitude

    @property
    def pole_latitude(self) -> float:
        """**read-only**

        Pole latitude in ecliptic J2000 coordinates (rad)
        """
        self._warn_if_rotation_incomplete()
        return self._pole_latitude

    @property
    def period(self) -> float:
        """**read-only**

        Rotation period (s)
        """
        self._warn_if_rotation_incomplete()
        return self._period

    @property
    def initial_rotation_angle(self) -> float:
        """**read-only**

        Rotation angle at the reference epoch (rad)
        """
        self._warn_if_rotation_incomplete()
        return self._initial_rotation_angle

    @property
    def reference_epoch(self) -> float:
        """**read-only**

        Rotation reference epoch in seconds since J2000
        """
        self._warn_if_rotation_incomplete()
        return self._reference_epoch

    @property
    def number_of_equally_good_models(self) -> int | None:
        """**read-only**

        Number of models (including this one) sharing the highest quality flag of the asteroid, when this model was
        selected by :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query` from an MPC number. None for
        models obtained in any other way.
        """
        return self._number_of_equally_good_models

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
            If several models share the highest flag, the one with the lowest model ID is selected, and their number is
            given by :attr:`~tudatpy.data_input.environment_data.damit.DamitModel.number_of_equally_good_models`.
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
        number_of_equally_good_models = None

        # Query by mpc number and choose best model
        if mpc_number and not damit_model_id:
            available_models = list_damit_models(
                mpc_number,
                force_refresh=force_refresh,
            )
            if available_models.empty:
                raise LookupError(f"No shape models found for {mpc_number}")

            # Highest quality flag (a missing flag ranks lowest), and lowest model ID (the table index) of those
            quality_flags = available_models["quality_flag"].fillna(-1)
            best_models = available_models[quality_flags == quality_flags.max()].sort_index()
            model_metadata = best_models.iloc[0]
            number_of_equally_good_models = len(best_models)

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

        model = cls.from_file(spin_file, shape_file, diameter)
        model._number_of_equally_good_models = number_of_equally_good_models

        # YORP is only available in the table, not in the spin file
        if pd.notna(model_metadata["yorp"]):
            model._rotation_warning = (
                f"DAMIT asteroid {model_metadata['asteroid_id']} (MPC number {model_metadata['mpc_number']}) "
                "has a YORP change in rotation rate, which is ignored: its rotation is modelled with a constant rate."
            )

        return model

    @classmethod
    def from_file(
        cls,
        spin_file: str | Path,
        shape_file: str | Path,
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
        ].to_numpy(dtype=float, copy=True)
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
            vertices *= 1e3  # converts km to m for models that are not scaled

        return vertices, facets

    def best_fitting_ellipsoid(self) -> tuple[list, float]:
        r"""
        Approximate the best fitting (in the sense of inertia) triaxial ellipsoid of the polyhedron shape model.

        The ellipsoid is solved for by solving the Eigenvalue problem:

        .. math::

            I\omega = \lambda\omega

        where :math:`I` is the inertia tensor of the polyhedron, :math:`\lambda` is an eigenvalue and :math:`\omega` is
        an eigenvector. The eigenvectors and eigenvalues represent principal axes of the ellipsoid and the moments of
        inertia about those axes, respectively. Rotation is assumed to be about the same axis as for the polyhedron
        (usually the axis of highest inertia). The x/y-axes of the ellipsoid may not necessarily align axes in which the
        polyhedron is defined, in which case an extra rotation about the z-axis is required (the principal axis offset).

        From the computed moments of inertia, the ellipsoid semi-axes are computed as:

        .. math::

            a = \sqrt{\frac{5}{2M}\left(I_B + I_C - I_A\right)}, \quad
            b = \sqrt{\frac{5}{2M}\left(I_A + I_C - I_B\right)}, \quad
            c = \sqrt{\frac{5}{2M}\left(I_A + I_B - I_C\right)}

        The principal axis offset :math:`\psi` is computed as follows, where :math:`v_A` is the eigenvector
        corresponding to eigenvalue :math:`I_A`

        .. math::

            \psi = \operatorname{atan2}\left(v_{A,y}, v_{A,x}\right)

        Returns
        -------
        tuple[list, float]
            Semi-axes [a, b, c] of the ellipsoid (m), and the principal axis offset around the z-axis (rad)
        """
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

    def initial_orientation(self, base_frame: str, extra_rotation: float = 0.0) -> np.ndarray:
        """
        Rotation matrix from the base frame to the body-fixed frame at the reference epoch.

        Parameters
        ----------
        base_frame : str
            Base frame used in the environment. Options: "J2000", "ECLIPJ2000"
        extra_rotation : float, optional
            Extra rotation angle around the body z-axis (rad), e.g. from the ellipsoid principal axis offset

        Returns
        -------
        np.ndarray
            3x3 rotation matrix from the base frame to the body-fixed frame
        """
        if base_frame not in ("J2000", "ECLIPJ2000"):
            raise ValueError("`base_frame` must be J2000 or ECLIPJ2000")

        self._warn_if_rotation_incomplete()

        # Initial orientation ecliptic to body-fixed. The prime meridian longitude is measured from the ascending node
        # of the body equator, which is 90 deg ahead of the direction from which DAMIT measures its rotation angle.
        rotation_to_body_fixed = inertial_to_body_fixed_rotation_matrix(
            pole_declination=self._pole_latitude,
            pole_right_ascension=self._pole_longitude,
            prime_meridian_longitude=self._initial_rotation_angle + extra_rotation - np.pi / 2,
        )
        if base_frame == "J2000":
            rotation_to_body_fixed = rotation_to_body_fixed @ j2000_to_eclipj2000()

        return rotation_to_body_fixed
