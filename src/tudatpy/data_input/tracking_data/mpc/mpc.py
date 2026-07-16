import pandas as pd
import numpy as np
from astroquery.mpc import MPC
import copy
from astropy.table import Table
from tudatpy.data_input.tracking_data.optical_utilities import (
    create_augmented_optical_table,
    filter_augmented_optical_table,
    optical_table_to_tracking_data,
    standardize_optical_dataframe,
)
from tudatpy.data_input.tracking_data.optical_utilities.optical_utilities import (
    _datetime_to_utc_seconds,
)
from tudatpy.data_input.tracking_data.obs_80_cols import unpackers
from tudatpy.data_input.tracking_data.obs_80_cols.parsers import parse_80cols_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    empty_radar_table,
    filter_radar_data,
    radar_data_from_table,
    radar_data_to_tracking_data,
    radar_station_id_list,
    radar_station_observation_counts,
    radar_target_body_list,
)

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
    asteroids and comets from the Minor Planet Center. The retrieved data are
    stored in a :class:`BatchMPC` instance and converted
    through :meth:`BatchMPC.to_tracking_dataset`. Optical observations can be
    augmented with default optical weights based on :cite:t:`veres2017`,
    star-catalog bias corrections based on :cite:t:`eggl2020`, and ancillary
    metadata.

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
    return batch.to_tracking_dataset(
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
        self._radar_table: pd.DataFrame = empty_radar_table()
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
        """Create a shallow-copy protocol implementation for BatchMPC."""
        new = BatchMPC()

        new._table = copy.deepcopy(self._table)
        new._radar_table = copy.deepcopy(self._radar_table)
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

    def filter(
        self,
        bands: list[str] | None = None,
        catalogs: list[str] | None = None,
        observation_types: list[str] | None = None,
        observatories: list[str] | None = None,
        observatories_exclude: list[str] | None = None,
        epoch_start=None,
        epoch_end=None,
        in_place: bool = True,
    ) -> "None | BatchMPC":
        """Filter observations in the batch.

        This method filters the augmented MPC optical table stored in
        :attr:`table` and the canonical radar table stored in
        :attr:`radar_table`. Epoch filters are interpreted as UTC seconds since
        J2000, or as datetime-like objects that can be converted to UTC seconds
        since J2000.

        Parameters
        ----------
        bands : list[str] | None, default None
            Observation bands to keep.
        catalogs : list[str] | None, default None
            Star-catalog codes to keep.
        observation_types : list[str] | None, default None
            MPC Note 2 observation types to keep. Radar filters may use MPC
            radar flags (``"R"``, ``"V"``) or the canonical radar observable
            type strings.
        observatories : list[str] | None, default None
            Observatory codes to keep.
        observatories_exclude : list[str] | None, default None
            Observatory codes to remove.
        epoch_start : float | datetime.datetime | None, default None
            Lower epoch bound. See
            :func:`~tudatpy.data_input.tracking_data.optical_utilities.filter_augmented_optical_table`
            for accepted epoch formats.
        epoch_end : float | datetime.datetime | None, default None
            Upper epoch bound. See
            :func:`~tudatpy.data_input.tracking_data.optical_utilities.filter_augmented_optical_table`
            for accepted epoch formats.
        in_place : bool, default True
            Whether to modify this batch. If ``False``, a filtered copy is
            returned.

        Returns
        -------
        None | BatchMPC
            ``None`` if ``in_place`` is ``True``; otherwise a filtered copy.
        """
        for name, value in {
            "bands": bands,
            "catalogs": catalogs,
            "observation_types": observation_types,
            "observatories": observatories,
            "observatories_exclude": observatories_exclude,
        }.items():
            if value is not None and not isinstance(value, list):
                raise ValueError(f"{name} must be a list or None")

        if observatories is not None and observatories_exclude is not None:
            raise ValueError("Include or exclude observatories, not both at the same time.")

        batch = self if in_place else self.copy()
        table = batch._table

        if bands is not None and "band" in table.columns:
            table = table.loc[table["band"].isin(bands)]
        if catalogs is not None and "catalog" in table.columns:
            table = table.loc[table["catalog"].isin(catalogs)]
        if observation_types is not None and "note2" in table.columns:
            table = table.loc[table["note2"].isin(observation_types)]

        if not table.empty:
            table = filter_augmented_optical_table(
                table,
                epoch_start=epoch_start,
                epoch_end=epoch_end,
                observatories=observatories,
                observatories_exclude=observatories_exclude,
            )

        if bands is not None or catalogs is not None:
            radar_table = empty_radar_table()
        else:
            radar_table = filter_radar_data(
                batch._radar_table,
                epoch_start=(
                    None if epoch_start is None else _datetime_to_utc_seconds(epoch_start)
                ),
                epoch_end=None if epoch_end is None else _datetime_to_utc_seconds(epoch_end),
                observable_type=self._radar_observable_types_from_mpc_filter(observation_types),
                station_ids=observatories,
                exclude_station_ids=observatories_exclude,
            )

        batch._table = table
        batch._radar_table = radar_table
        batch._refresh_metadata()

        if in_place:
            return None
        return batch

    def summary(self) -> None:
        """Print a short summary of the current batch."""
        print()
        print("   Batch Summary:")
        print(f"1. Batch includes {len(self._MPC_codes)} minor planets:")
        print("  ", self.MPC_objects)
        print(f"2. Batch includes {self.size} observations.")

        if self._table.empty and self._radar_table.empty:
            print("3. The batch contains no observations.")
            print()
            return

        print(f"3. The observations range from {self.epoch_start} to {self.epoch_end}")
        if not self._table.empty:
            print(f"   In Julian Days: {self._table.epoch.min()} to {self._table.epoch.max()}")
        print(
            f"4. The batch contains observations from {len(self.observatories)} "
            + f"observatories, including {len(self.space_telescopes)} space telescopes"
        )
        print()

    def observatories_table(
        self,
        only_in_batch: bool = True,
        only_space_telescopes: bool = False,
        exclude_space_telescopes: bool = False,
        include_positions: bool = False,
    ) -> pd.DataFrame:
        """Return a table with observatory counts for this batch.

        Parameters
        ----------
        only_in_batch : bool, default True
            Whether to return only observatories present in this batch. Since
            the new MPC interface does not maintain a separate global MPC
            observatory catalog, this argument is retained for compatibility
            and has no effect when ``True``.
        only_space_telescopes : bool, default False
            Whether to return only observatories marked as space telescopes in
            the batch metadata.
        exclude_space_telescopes : bool, default False
            Whether to remove observatories marked as space telescopes in the
            batch metadata.
        include_positions : bool, default False
            Retained for compatibility. The current table does not include
            observatory positions.

        Returns
        -------
        pandas.DataFrame
            Dataframe with columns ``Code``, ``Name`` and ``count``.
        """
        if self._table.empty and self._radar_table.empty:
            return pd.DataFrame(columns=["Code", "Name", "count"])

        frames = []
        if not self._table.empty:
            frames.append(
                self._table.groupby("observatory")
                .size()
                .rename("count")
                .reset_index()
                .rename(columns={"observatory": "Code"})
            )
        if not self._radar_table.empty:
            frames.append(
                radar_station_observation_counts(self._radar_table).rename(
                    columns={"observatory": "Code"}
                )
            )
        table = (
            pd.concat(frames, ignore_index=True)
            .groupby("Code", as_index=False)["count"]
            .sum()
            .assign(Name=lambda x: x["Code"])
        )

        if only_space_telescopes:
            table = table.loc[table["Code"].isin(self.space_telescopes)]
        if exclude_space_telescopes:
            table = table.loc[~table["Code"].isin(self.space_telescopes)]

        columns = ["Code", "Name", "count"]
        return table.loc[:, columns].reset_index(drop=True)

    # getters to make everything read-only
    @property
    def table(self) -> pd.DataFrame:
        """**read-only**

        Pandas dataframe with observation data.

        :type: pd.DataFrame
        """
        return self._table

    @property
    def radar_table(self) -> pd.DataFrame:
        """**read-only**

        Pandas dataframe with canonical radar tracking data.

        :type: pd.DataFrame
        """
        return self._radar_table.copy()

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
        """Return the total number of optical and radar observations."""
        return self._size

    def __add__(self, other):
        """Return a new batch containing observations from both operands."""
        temp = BatchMPC()

        temp._table = pd.concat([self._table, other._table], ignore_index=True, sort=False)
        if "epoch" in temp._table.columns:
            temp._table = temp._table.sort_values("epoch")  # this is expressed in Julian Days
        temp._table = temp._table.drop_duplicates()
        temp._radar_table = pd.concat(
            [self._radar_table, other._radar_table],
            ignore_index=True,
            sort=False,
        ).drop_duplicates()

        temp._refresh_metadata()

        return temp

    # helper functions
    def _refresh_metadata(self) -> None:
        """Internal. Update batch metadata."""
        self._table = self._table.drop_duplicates()
        self._radar_table = self._radar_table.drop_duplicates()

        if self._table.empty and self._radar_table.empty:
            self._observatories = []
            self._space_telescopes = []
            self._bands = []
            self._MPC_codes = []
            self._size = 0
            self._epoch_start = 0.0
            self._epoch_end = 0.0
            return

        optical_observatories = (
            self._table["observatory"].dropna().astype(str).tolist()
            if "observatory" in self._table
            else []
        )
        radar_observatories = radar_station_id_list(self._radar_table)
        self._observatories = list(dict.fromkeys(optical_observatories + radar_observatories))
        self._space_telescopes = []
        if "band" in self._table.columns:
            self._bands = list(self._table.band.unique())

        # if user gives custom name, set that as body name, else MPC code
        if "custom_name" in self._table.columns and self._table["custom_name"].notna().any():
            optical_targets = self._table["custom_name"].dropna().astype(str).tolist()
        else:
            optical_targets = (
                self._table["number"].dropna().astype(str).tolist()
                if "number" in self._table
                else []
            )
        radar_targets = radar_target_body_list(self._radar_table)
        self._MPC_codes = list(dict.fromkeys(optical_targets + radar_targets))
        self._size = len(self._table) + len(self._radar_table)

        epoch_sources = []
        if "epoch_seconds_UTC" in self._table.columns and not self._table.empty:
            epoch_sources.extend(
                [
                    self._table.epoch_seconds_UTC.min(),
                    self._table.epoch_seconds_UTC.max(),
                ]
            )
        elif "epoch" in self._table.columns and not self._table.empty:
            epoch_sources.extend([self._table.epoch.min(), self._table.epoch.max()])
        if not self._radar_table.empty:
            epoch_sources.extend(
                [
                    self._radar_table.epoch_seconds_UTC.min(),
                    self._radar_table.epoch_seconds_UTC.max(),
                ]
            )
        self._epoch_start = min(epoch_sources) if epoch_sources else 0.0
        self._epoch_end = max(epoch_sources) if epoch_sources else 0.0

    def _add_custom_name_column(self, table: pd.DataFrame, custom_name) -> pd.DataFrame:
        """Return a copy of an optical table with a ``custom_name`` column."""
        augmented_table = table.copy()
        if custom_name is not None or "custom_name" not in augmented_table.columns:
            augmented_table["custom_name"] = [custom_name] * len(augmented_table)
        return augmented_table

    @staticmethod
    def _radar_observable_types_from_mpc_filter(observation_types: list[str] | None):
        """Map MPC Note 2 radar flags to canonical radar observable names."""
        if observation_types is None:
            return None
        type_map = {
            "R": RANGE_OBSERVABLE,
            "r": RANGE_OBSERVABLE,
            RANGE_OBSERVABLE: RANGE_OBSERVABLE,
            "V": DOPPLER_OBSERVABLE,
            "v": DOPPLER_OBSERVABLE,
            DOPPLER_OBSERVABLE: DOPPLER_OBSERVABLE,
        }
        return [type_map.get(str(value), str(value)) for value in observation_types]

    def _add_parsed_mpc80_table(self, parsed_table: Table, custom_name: str | None = None) -> None:
        """Internal. Add a parsed MPC80 table, including radar metadata if present."""
        if len(parsed_table) > 0:
            # Optical and space-astrometry rows are converted to the shared
            # augmented optical representation. Radar rows are held separately
            # as canonical radar DataFrames.
            optical_table = create_augmented_optical_table(
                parsed_table,
                in_degrees=False,
                frame="J2000",
                custom_name=custom_name,
            )
            self._table = pd.concat([self._table, optical_table], ignore_index=True, sort=False)

        radar_table = radar_data_from_table(parsed_table)
        if custom_name is not None and not radar_table.empty:
            radar_table["target_body"] = str(custom_name)
        if not radar_table.empty:
            if self._radar_table.empty:
                self._radar_table = radar_table.copy()
            else:
                self._radar_table = pd.concat(
                    [self._radar_table, radar_table],
                    ignore_index=True,
                    sort=False,
                )

    @staticmethod
    def _fetch_mpc80_records(code: str | int, id_type: str | None = None) -> list[str]:
        """Retrieve raw MPC 80-column records through astroquery."""
        query_kwargs = {"get_mpcformat": True}
        if id_type is not None:
            query_kwargs["id_type"] = id_type
        raw_observations = MPC.get_observations(code, **query_kwargs)
        if "obs" not in raw_observations.colnames:
            raise RuntimeError("MPC get_mpcformat query did not return an 'obs' column.")
        return [str(row) for row in raw_observations["obs"]]

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

            if not drop_misc_observations:
                parsed_table = parse_80cols_data(self._fetch_mpc80_records(code, id_type))
                self._add_parsed_mpc80_table(parsed_table, custom_name=custom_name)
                continue

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
        if self._table.empty:
            optical_tracking_data, optical_supplementary_data = [], []
        else:
            optical_tracking_data, optical_supplementary_data = optical_table_to_tracking_data(
                self._table,
                add_weights=add_weights,
                add_star_catalog_corrections=add_star_catalog_corrections,
                add_ancillary_data=add_ancillary_data,
            )
        radar_tracking_data, radar_supplementary_data = radar_data_to_tracking_data(
            self._radar_table
        )
        return (
            optical_tracking_data + radar_tracking_data,
            optical_supplementary_data + radar_supplementary_data,
        )
