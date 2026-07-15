import pandas as pd
import numpy as np
from astroquery.mpc import MPC
import copy
import importlib
from tudatpy.data_input.tracking_data.optical_utilities import (
    create_augmented_optical_table,
    optical_table_to_tracking_data,
    standardize_optical_dataframe,
)

unpackers = importlib.import_module("tudatpy.data_input.tracking_data.80_column.unpackers")
OBS_TYPES_TO_DROP = unpackers.OBS_TYPES_TO_DROP


def read_mpc_data(
    MPCcodes: list[str | int],
    id_types: list[str | None] | None = None,
    drop_misc_observations: bool = True,
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
):
    """Retrieve MPC observations and return tracking-data containers.

    This function is a convenience interface around :class:`BatchMPC`, which
    uses the ``astroquery`` MPC interface to retrieve optical observations for
    asteroids and comets from the Minor Planet Center. The retrieved
    :attr:`BatchMPC.table` is standardized and passed to
    :func:`~tudatpy.data_input.tracking_data.optical_utilities.optical_table_to_tracking_data`,
    which creates Tudat optical tracking data through the common optical-data
    conversion pipeline. During this conversion the data can be augmented with
    default optical weights based on :cite:t:`veres2017`, star-catalog bias
    corrections based on :cite:t:`eggl2020`, and ancillary metadata.

    Parameters
    ----------
    MPCcodes : list[str | int]
        MPC object identifiers for which observations are retrieved.
    id_types : list[str | None] | None, default None
        Optional identifier type for each MPC object identifier.
    drop_misc_observations : bool, default True
        Whether to drop MPC records with Note 2 flags for replaced discovery
        observations (``x``, ``X``), roving observers (``V``, ``v``, ``W``,
        ``w``), radar observations (``R``, ``r``, ``Q``, ``q``), offset
        observations of natural satellites (``O``), and satellite/space-based
        records (``S``, ``s``, ``T``, ``t``). This filtering is based on the
        MPC Note 2 flag, not on the observatory code.
    custom_name : str | None, default None
        Optional body name used in the generated tracking data.
    add_weights : bool, default False
        Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
    add_star_catalog_corrections : bool, default False
        Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
    add_ancillary_data : bool, default False
        Whether to attach available optical ancillary metadata.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data objects and supplementary data objects.
    """
    batch = BatchMPC()
    batch.get_observations(
        MPCcodes,
        id_types=id_types,
        drop_misc_observations=drop_misc_observations,
        custom_name=custom_name,
    )
    return optical_table_to_tracking_data(
        batch.table,
        add_weights,
        add_star_catalog_corrections,
        add_ancillary_data,
    )


class BatchMPC:
    """Interface between MPC optical observations and Tudat tracking data.

    This class wraps the MPC interface of ``astroquery`` and provides
    Tudat-specific processing for asteroid and comet observations, including
    optional observation weights based on :cite:t:`veres2017` and
    star-catalog bias corrections based on :cite:t:`eggl2020`.

    Notes
    ----------
    Currently, transformations between reference frames are not implemented.
    As such observations are only returned in the J2000 Frame.

    Examples
    ----------
    Basic sequence of usage:

    Initialise and retrieve data:

    >>> MPCcodes = [1, 4] # Ceres and Vesta
    >>> batch = BatchMPC()
    >>> batch.get_observations(MPCcodes)

    Transform to Tudat tracking data:
    >>> tracking_data, supplementary_data = batch.to_tracking_dataset()

    For simple loading, users can call read_mpc_data(...) without explicitly
    constructing a BatchMPC instance.

    """

    def __init__(self) -> None:
        """Create an empty MPC batch."""
        self._table: pd.DataFrame = pd.DataFrame()
        self._observatories: list[str] = []
        self._space_telescopes: list[str] = []
        self._bands: list[str] = []
        self._MPC_codes: list[str] = []
        self._size: int = 0

        self._epoch_start: float = 0.0
        self._epoch_end: float = 0.0

        self._bodies_created = {}

        self._EFCC18_applied = False
        self._custom_weights_set = False

    def __copy__(self):
        new = BatchMPC()

        new._table = copy.deepcopy(self._table)
        new._refresh_metadata()

        return new

    def copy(self) -> "BatchMPC":
        """Create a copy of this MPC batch.

        Parameters
        ----------
        None

        Returns
        -------
        BatchMPC
            Copy of batch.
        """
        return copy.copy(self)

    # getters to make everything read-only
    @property
    def table(self) -> pd.DataFrame:
        """**read-only**

        Pandas dataframe with observation data.

        :type: pd.DataFrame
        """
        return self._table

    @property
    def observatories(self) -> list[str]:
        """**read-only**

        List of observatories in batch.

        :type: list[str]
        """
        return self._observatories

    @property
    def space_telescopes(self) -> list[str]:
        """**read-only**

        List of satellite observatories in batch.

        :type: list[str]
        """
        return self._space_telescopes

    @property
    def MPC_objects(self) -> list[str]:
        """**read-only**

        List of MPC objects.

        :type: list[str]
        """
        return self._MPC_codes

    @property
    def size(self) -> int:
        """**read-only**

        Number of observations in batch.

        :type: int
        """
        return self._size

    @property
    def bands(self) -> list[str]:
        """**read-only**

        List of bands in batch.

        :type: list[str]
        """
        return self._bands

    @property
    def epoch_start(self) -> float:
        """**read-only**

        Epoch of oldest observation in batch.

        :type: float
        """
        return self._epoch_start

    @property
    def epoch_end(self) -> float:
        """**read-only**

        Epoch of latest observation in batch.

        :type: float
        """
        return self._epoch_end

    @property
    def bodies_created(self) -> dict:
        """**read-only**

        Dictionary with generated-body bookkeeping.

        :type: dict
        """
        return self._bodies_created

    def __len__(self):
        return self._size

    def __add__(self, other):
        temp = BatchMPC()

        temp._table = (
            pd.concat([self._table, other._table])
            .sort_values("epoch")  # this is expressed in Julian Days
            .drop_duplicates()
        )

        temp._refresh_metadata()

        return temp

    # helper functions
    def _refresh_metadata(self) -> None:
        """Internal. Update batch metadata."""
        self._table.drop_duplicates()

        if self._table.empty:
            self._observatories = []
            self._space_telescopes = []
            self._bands = []
            self._MPC_codes = []
            self._size = 0
            self._epoch_start = 0.0
            self._epoch_end = 0.0
            return

        self._observatories = list(self._table.observatory.unique())
        self._space_telescopes = []
        if "band" in self._table.columns:
            self._bands = list(self._table.band.unique())

        # if user gives custom name, set that as body name, else MPC code
        if "custom_name" in self._table.columns and self._table["custom_name"].notna().any():
            self._MPC_codes = list(self._table["custom_name"].unique())
        else:
            self._MPC_codes = list(self._table.number.unique())
        self._size = len(self._table)

        self._epoch_start = self._table.epoch.min()
        self._epoch_end = self._table.epoch.max()

    def _add_custom_name_column(self, table: pd.DataFrame, custom_name) -> pd.DataFrame:
        augmented_table = table.copy()
        if custom_name is not None or "custom_name" not in augmented_table.columns:
            augmented_table["custom_name"] = [custom_name] * len(augmented_table)
        return augmented_table

    ###########################################################################################
    # MPC Astroquery Data Retrieval: get_observations
    ###########################################################################################
    def get_observations(
        self,
        MPCcodes: list[str | int],
        id_types: list[str | None] | None = None,
        drop_misc_observations: bool = True,
        custom_name: str | None = None,
    ) -> None:
        """Retrieve all observations for a set of MPC listed objects.

        This method uses astroquery to retrieve the observations from the MPC.
        An internet connection is required, observations are cached for faster subsequent retrieval.
        Removes duplicate observations, and can filter selected MPC Note 2
        observation types before the data are converted to the common augmented
        optical table used by Tudat.

        Parameters
        ----------
        MPCcodes : list[str | int]
            List of integer or str MPC object codes for minor planets or comets.
        id_types : list[str | None] | None, default None
            A list of identification types ('asteroid_number', 'comet_number', 'comet_designation') corresponding to each MPCcode.
            If an element is None, the type is considered unknown. If the entire list is None,
            all types are considered unknown.
        drop_misc_observations : bool, default True
            Whether to drop MPC records with Note 2 flags for replaced
            discovery observations (``x``, ``X``), roving observers (``V``,
            ``v``, ``W``, ``w``), radar observations (``R``, ``r``, ``Q``,
            ``q``), offset observations of natural satellites (``O``), and
            satellite/space-based records (``S``, ``s``, ``T``, ``t``). This
            filtering is based on the MPC Note 2 flag, not on the observatory
            code.

        Raises
        ------
        ValueError
            If a valid identifier (a numbered designation, or a provisional
            designation) cannot be determined for one of the provided codes,
            indicating that the code is likely invalid.
        """

        if not isinstance(MPCcodes, list):
            raise ValueError("MPCcodes parameter must be a list of integers/strings")

        # 1. If id_types is not provided, create a list of Nones
        if id_types is None:
            id_types = [None] * len(MPCcodes)

        # 2. Ensure the lists have the same length for a 1-to-1 mapping
        if len(MPCcodes) != len(id_types):
            raise ValueError("MPCcodes and id_types must have the same number of elements.")

        for code, id_type in zip(MPCcodes, id_types):
            if not (isinstance(code, int) or isinstance(code, str)):
                raise ValueError("All codes in the MPCcodes parameter must be integers or strings")

            # 3. Conditionally call the function based on whether id_type exists
            if id_type is not None:
                obs = MPC.get_observations(code, id_type=id_type).to_pandas()
            else:
                obs = MPC.get_observations(code).to_pandas()

            obs = standardize_optical_dataframe(obs)
            obs = self._add_custom_name_column(obs, custom_name)

            # convert deg to rad, wrapping as tudat wants [-pi, pi]
            obs = obs.assign(
                RA=lambda x: (np.radians(x.RA) + np.pi) % (2 * np.pi) - np.pi,
                DEC=lambda x: np.radians(x.DEC),
            )

            identifier = None
            if drop_misc_observations:
                obs = obs.query("note2 not in @OBS_TYPES_TO_DROP")

                # Check for Comets/Interstellars (Astroquery returns 'comet_type' or 'comettype')
                # If we have a number and a type, combine them (e.g., 3 + I = 3I)
                type_col = None
                if "comet_type" in obs.columns:
                    type_col = "comet_type"
                elif "comettype" in obs.columns:
                    type_col = "comettype"

                if type_col and pd.notna(obs[type_col].iloc[0]):  # checks first digit is not NA
                    # It is a comet or interstellar object
                    number_part = str(obs["number"].iloc[0])
                    type_part = str(obs[type_col].iloc[0])
                    identifier = f"{number_part}{type_part}"  # Result: "3I"

                elif "number" in obs.columns:
                    pd.set_option("future.no_silent_downcasting", True)
                    valid_numbers = (
                        obs["number"].dropna().astype(str).replace("<NA>", np.nan).dropna()
                    )

                    if not valid_numbers.empty:
                        potential_id = valid_numbers.iloc[0]
                    else:
                        # fallback to designation if no number has been assigned yet
                        valid_designations = (
                            obs["desig"].dropna().astype(str).replace("<NA>", np.nan).dropna()
                        )
                        potential_id = (
                            valid_designations.iloc[0] if not valid_designations.empty else None
                        )

                    if potential_id is not None:
                        # We allow alphanumeric strings now (to support packed numbers like D4341)
                        # We only pad if it is a short digit string (e.g. '1' -> '00001')
                        # Packed strings are always 5 chars long, so they won't be affected by zfill(5)
                        if len(potential_id) < 5:
                            potential_id = potential_id.zfill(5)

                        try:
                            # Try to unpack it. This handles '00001' and 'D4341'.
                            identifier = unpackers.unpack_permanent_minor_planet(potential_id)
                        except Exception:
                            # If unpacking fails (e.g. it was already unpacked or invalid),
                            # we keep the potential_id as is.
                            identifier = potential_id

            if identifier is None and "desig" in obs.columns and pd.notna(obs["desig"].iloc[0]):
                identifier = str(obs["desig"].iloc[0])

            if identifier is None:
                raise ValueError(
                    f"Could not find a valid identifier (number or designation) "
                    f"for object code {code}. The provided code is likely invalid."
                )

            # Assign the identifier to the 'number' column for the entire DataFrame.
            obs.loc[:, "number"] = identifier
            obs = create_augmented_optical_table(
                obs, in_degrees=False, frame="J2000", custom_name=custom_name
            )
            self._table = pd.concat([self._table, obs])

        self._refresh_metadata()

    def to_tracking_dataset(
        self,
        add_weights: bool | None = False,
        add_star_catalog_corrections: bool | None = False,
        add_ancillary_data: bool | None = False,
    ):
        """Convert the current MPC batch to tracking-data containers.

        Parameters
        ----------
        add_weights : bool, default False
            Whether to assign the default optical weighing scheme of :cite:p:`veres2017`.
        add_star_catalog_corrections : bool, default False
            Whether to attach star-catalog bias corrections following :cite:p:`eggl2020`.
        add_ancillary_data : bool, default False
            Whether to attach available optical ancillary metadata.

        Returns
        -------
        tuple[list[TrackingData], list[TrackingSupplementaryData]]
            Tracking data objects and supplementary data objects.
        """
        return optical_table_to_tracking_data(
            self._table,
            add_weights=add_weights,
            add_star_catalog_corrections=add_star_catalog_corrections,
            add_ancillary_data=add_ancillary_data,
        )
