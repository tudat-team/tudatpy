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
    """Download the DAMIT asteroid and asteroid model tables, which give the MPC numbers, sizes and YORP rates of the
    models."""
    for table in (_DAMIT_ASTEROIDS, _DAMIT_ASTEROID_MODELS):
        download_extra_file(
            url=f"https://damit.cuni.cz/projects/damit/exports/table/{table.stem}",
            dest_path=table,
            cache_dir=DEFAULT_CACHE,
            force=force,
        )


def _download_model(model_id: int, force: bool = False) -> tuple[Path, Path]:
    """Download the shape (shape.obj) and rotation (spin.txt) files of a DAMIT model, and return their paths."""
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
    """Load the table of all DAMIT models, with the MPC number of each model's asteroid, indexed by model ID."""
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
    """List the DAMIT models of an asteroid, or of all asteroids.

    An asteroid can have several models, or none, in which case the returned table is empty. The DAMIT tables are
    downloaded on the first call, and stored for later use.

    Parameters
    ----------
    mpc_number : int, optional
        MPC number of the asteroid. If None, all models are listed.
    force_refresh : bool, optional
        If True, downloads the DAMIT tables again, replacing the stored ones. Defaults to False.

    Returns
    -------
    pandas.DataFrame
        Table of the models, indexed by DAMIT model ID, with the columns of the DAMIT model table and the MPC number
        of the asteroid (``mpc_number``).
    """
    table = _load_model_table(force_refresh)

    if mpc_number is None:
        return table

    query_result = table[table["mpc_number"] == mpc_number]
    return query_result


class DamitModel:
    """Asteroid model from DAMIT: a polyhedron shape and a rotation about a fixed axis.

    A model is normally created with :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query` or
    :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.from_file`. Shape model settings are created from it
    with :func:`~tudatpy.dynamics.environment_setup.shape.damit_polyhedron` or
    :func:`~tudatpy.dynamics.environment_setup.shape.damit_best_fitting_ellipsoid`, and rotation model settings with
    :func:`~tudatpy.dynamics.environment_setup.rotation_model.damit`.

    The body-fixed frame of the model has its z-axis along the rotation axis, which is assumed to be a principal axis,
    and its x-axis along the principal axis of minimum moment of inertia in the equatorial plane (for a constant
    density). The vertices are rotated about the z-axis to this frame, and the rotation angle is changed by the same
    angle, so that the orientation of the body is unchanged. Body-fixed longitudes therefore differ from those in the
    frame of the DAMIT shape model by this angle.

    The rotation rate is constant: a YORP change of the rotation rate is ignored. For a model from
    :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query` that has one, a warning is given when its
    rotation is used.

    Parameters
    ----------
    vertices : numpy.ndarray
        Cartesian coordinates of the vertices in the frame of the shape model (m), one row per vertex.
    facets : numpy.ndarray
        Indices of the three vertices of each facet (zero-indexed), one row per facet.
    pole_longitude : float
        Ecliptic longitude of the rotation axis, w.r.t. ECLIPJ2000 (rad).
    pole_latitude : float
        Ecliptic latitude of the rotation axis, w.r.t. ECLIPJ2000 (rad).
    period : float
        Sidereal rotation period (s).
    initial_rotation_angle : float
        Rotation angle of the frame of the shape model at the reference epoch (rad).
    reference_epoch : float
        Reference epoch of the rotation, in seconds since J2000.
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
        # Copy the input, so that later changes to the input arrays do not change the model
        vertices = np.array(vertices, dtype=float)
        facets = np.array(facets, dtype=int)

        # Input check
        if not vertices.shape[1] == 3 or not facets.shape[1] == 3:
            raise ValueError("Vertices and facets must have three columns")
        if facets.max() > len(vertices) - 1 or facets.min() < 0:
            raise ValueError("Facets indices out of range")

        if period <= 0:
            raise ValueError("Period must be positive")

        # Rotate the vertices about the z-axis by the smallest angle that puts the x-axis along the principal axis of
        # minimum moment of inertia in the equatorial plane (r' = R_z(-offset) r), and add this angle to the rotation
        # angle, so that the orientation of the body is unchanged
        _, eigenvectors = np.linalg.eigh(inertia_tensor_from_density(vertices, facets, 1)[:2, :2])
        principal_axis_offset = (
            np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]) + np.pi / 2
        ) % np.pi - np.pi / 2
        cos_offset, sin_offset = np.cos(principal_axis_offset), np.sin(principal_axis_offset)
        vertices[:, :2] = vertices[:, :2] @ np.array(
            [[cos_offset, -sin_offset], [sin_offset, cos_offset]]
        )
        initial_rotation_angle += principal_axis_offset

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

    @property
    def vertices(self) -> np.ndarray:
        """**read-only**

        Copy of the body-fixed Cartesian coordinates of the vertices (m), one row per vertex.

        :type: numpy.ndarray
        """
        return self._vertices.copy()

    @property
    def facets(self) -> np.ndarray:
        """**read-only**

        Copy of the indices of the three vertices of each facet (zero-indexed), one row per facet.

        :type: numpy.ndarray
        """
        return self._facets.copy()

    @property
    def pole_longitude(self) -> float:
        """**read-only**

        Ecliptic longitude of the rotation axis, w.r.t. ECLIPJ2000 (rad).

        :type: float
        """
        return self._pole_longitude

    @property
    def pole_latitude(self) -> float:
        """**read-only**

        Ecliptic latitude of the rotation axis, w.r.t. ECLIPJ2000 (rad).

        :type: float
        """
        return self._pole_latitude

    @property
    def period(self) -> float:
        """**read-only**

        Sidereal rotation period (s).

        :type: float
        """
        return self._period

    @property
    def initial_rotation_angle(self) -> float:
        """**read-only**

        Rotation angle of the body-fixed frame at the reference epoch (rad).

        :type: float
        """
        return self._initial_rotation_angle

    @property
    def reference_epoch(self) -> float:
        """**read-only**

        Reference epoch of the rotation, in seconds since J2000.

        :type: float
        """
        return self._reference_epoch

    @property
    def number_of_equally_good_models(self) -> int | None:
        """**read-only**

        Number of models of the asteroid (including this one) that share the highest quality flag, if this model was
        selected by :meth:`~tudatpy.data_input.environment_data.damit.DamitModel.query` from an MPC number, and None
        otherwise.

        :type: int | None
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
        """Load a DAMIT model by MPC number or by DAMIT model ID.

        Exactly one of ``mpc_number`` and ``damit_model_id`` must be given. The DAMIT tables and model files are
        downloaded on first use, and stored for later use.

        Parameters
        ----------
        mpc_number : int, optional
            MPC number of the asteroid. Of its models, the one with the highest quality flag is selected. If several
            models share the highest flag, the one with the lowest model ID is selected, and their number is given by
            :attr:`~tudatpy.data_input.environment_data.damit.DamitModel.number_of_equally_good_models`.
        damit_model_id : int, optional
            DAMIT model ID.
        diameter : float, optional
            Volume-equivalent diameter (m) to which the shape is scaled. If not given, the size from the DAMIT tables
            is used. If DAMIT has no size information for the model, a warning is given, and the size is not physical.
        force_refresh : bool, optional
            If True, downloads the DAMIT tables and model files again, replacing the stored ones. Defaults to False.

        Returns
        -------
        DamitModel
            Selected DAMIT model.
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
            raise ValueError(
                "Must specify either a valid MPC number or DAMIT model ID (and not both)"
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
        model_id = int(model_metadata["id"])
        shape_file, spin_file = _download_model(model_id, force_refresh)

        model = cls.from_file(spin_file, shape_file, diameter)
        model._number_of_equally_good_models = number_of_equally_good_models

        # YORP is only available in the table, not in the spin file
        if pd.notna(model_metadata["yorp"]):
            model._rotation_warning = (
                f"DAMIT asteroid {model_metadata['asteroid_id']} (MPC number {model_metadata['mpc_number']}) "
                "has a YORP change in rotation rate, which is ignored: its rotation is modeled with a constant rate."
            )

        return model

    @classmethod
    def from_file(
        cls,
        spin_file: str | Path,
        shape_file: str | Path,
        diameter: float | None = None,
    ):
        """Load a DAMIT model from files in the format of those downloaded from the DAMIT database.

        Parameters
        ----------
        spin_file : str or pathlib.Path
            Path to the rotation file (spin.txt).
        shape_file : str or pathlib.Path
            Path to the shape file, in OBJ format (shape.obj).
        diameter : float, optional
            Volume-equivalent diameter (m) to which the shape is scaled. If not given, the vertex coordinates are taken
            to be in km, which gives a non-physical size for a model that DAMIT gives with unit volume.

        Returns
        -------
        DamitModel
            DAMIT model from the files.
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
        """Load the vertices (m) and zero-indexed facets of an OBJ shape file, scaled to the given volume-equivalent
        diameter (m), or converted from km if no diameter is given."""
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
            if not diameter > 0:
                raise ValueError("Diameter must be positive")
            current_model_volume = volume(vertices, facets)
            if current_model_volume <= 0:
                raise RuntimeError("Broken shape model: shape model has negative volume")
            required_volume = (4 / 3) * np.pi * (diameter / 2) ** 3
            vertices *= np.cbrt(required_volume / current_model_volume)
        else:
            vertices *= 1e3  # converts km to m for models that are not scaled

        return vertices, facets

    def best_fitting_ellipsoid(self) -> list[float]:
        r"""Compute the tri-axial ellipsoid with the same moments of inertia as the polyhedron shape.

        The axes of the ellipsoid are those of the body-fixed frame, which are principal axes of the polyhedron. With
        the moments of inertia :math:`I_A \le I_B` about the x- and y-axes, :math:`I_C` about the z-axis, and the mass
        :math:`M` of the polyhedron, the semi-axes are:

        .. math::

            a = \sqrt{\frac{5}{2M}\left(I_B + I_C - I_A\right)}, \quad
            b = \sqrt{\frac{5}{2M}\left(I_A + I_C - I_B\right)}, \quad
            c = \sqrt{\frac{5}{2M}\left(I_A + I_B - I_C\right)}

        The inertia tensor is computed for a constant density, about the origin of the model, which is assumed to be
        its center of mass.

        Returns
        -------
        list[float]
            Semi-axes :math:`[a, b, c]` of the ellipsoid along the x-, y- and z-axes (m).
        """
        ixx, iyy, izz = np.diag(inertia_tensor_from_density(self._vertices, self._facets, 1))

        # Convert moments of inertia to ellipsoid axes
        mass = volume(self._vertices, self._facets)  # Model has unit density
        a = np.sqrt((5.0 / (2.0 * mass)) * (iyy + izz - ixx))
        b = np.sqrt((5.0 / (2.0 * mass)) * (ixx + izz - iyy))
        c = np.sqrt((5.0 / (2.0 * mass)) * (ixx + iyy - izz))

        return [a, b, c]

    def initial_orientation(self, base_frame: str) -> np.ndarray:
        """Compute the rotation matrix from the base frame to the body-fixed frame at the reference epoch.

        If DAMIT gives a YORP change in rotation rate for the model, a warning is given, as this change is ignored.

        Parameters
        ----------
        base_frame : str
            Base frame orientation of the environment: ``"J2000"`` or ``"ECLIPJ2000"``.

        Returns
        -------
        numpy.ndarray
            Rotation matrix from the base frame to the body-fixed frame.
        """
        if base_frame not in ("J2000", "ECLIPJ2000"):
            raise ValueError("`base_frame` must be J2000 or ECLIPJ2000")

        if self._rotation_warning is not None:
            warn(self._rotation_warning)

        # Initial orientation ecliptic to body-fixed. The prime meridian longitude is measured from the ascending node
        # of the body equator, which is 90 deg ahead of the direction from which DAMIT measures its rotation angle.
        rotation_to_body_fixed = inertial_to_body_fixed_rotation_matrix(
            pole_declination=self._pole_latitude,
            pole_right_ascension=self._pole_longitude,
            prime_meridian_longitude=self._initial_rotation_angle - np.pi / 2,
        )
        if base_frame == "J2000":
            rotation_to_body_fixed = rotation_to_body_fixed @ j2000_to_eclipj2000()

        return rotation_to_body_fixed
