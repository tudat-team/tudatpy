import pandas as pd
import numpy as np
from astroquery.mpc import MPC
import copy
from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data.optical_utilities import (
    create_augmented_optical_table,
    filter_augmented_optical_table,
    optical_table_to_tracking_data,
    standardize_optical_dataframe,
)
from tudatpy.dynamics import environment_setup
from tudatpy.data_input.tracking_data.obs_80_cols import unpackers
from tudatpy.data_input.tracking_data.obs_80_cols.parsers import parse_80cols_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    empty_radar_table,
    filter_radar_data,
    radar_data_from_table,
    radar_data_to_tracking_data,
)

from tudatpy.data_input.tracking_data.optical_utilities.optical_utilities import (
    _resolve_optical_target_names,
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
    use_mpc80_format: bool = False,
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
    use_mpc80_format : bool, default False
        Retrieve raw MPC 80-column records instead of Astroquery's parsed table.
        This is required for radar and space-based observations.
        Drop-misc filtering is ignored in this mode.

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
        use_mpc80_format=use_mpc80_format,
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

    _SPACECRAFT_POSITION_COLUMNS = (
        "spacecraft_position_x",
        "spacecraft_position_y",
        "spacecraft_position_z",
    )

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
            MPC Note 2 observation types to keep. ``"R"`` keeps all radar
            observations; ``"NWayRange"`` or ``"DopplerMeasuredFrequency"``
            keeps one radar observable. Radar observations are removed when
            ``bands`` or ``catalogs`` is given, since they have neither.
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
                epoch_start=epoch_start,
                epoch_end=epoch_end,
                observable_type=self._radar_observable_types(observation_types),
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
        """Return MPC observatory names and observation counts for this batch.

        Observatory metadata is retrieved only when this method is called;
        constructing or loading a batch does not require the station catalog.

        Parameters
        ----------
        only_in_batch : bool, default True
            Whether to return only observatories present in this batch.
        only_space_telescopes : bool, default False
            Whether to return only observatories marked as space telescopes in
            the MPC observatory catalog.
        exclude_space_telescopes : bool, default False
            Whether to remove observatories marked as space telescopes in the
            MPC observatory catalog.
        include_positions : bool, default False
            Retained for compatibility. The current table does not include
            observatory positions.

        Returns
        -------
        pandas.DataFrame
            Dataframe with columns ``Code``, ``Name`` and ``count``.
        """
        if self._table.empty and self._radar_table.empty and only_in_batch:
            return pd.DataFrame(columns=["Code", "Name", "count"])

        catalog = MPC.get_observatory_codes().to_pandas()
        catalog["Code"] = catalog["Code"].astype(str).str.strip().str.zfill(3)
        space_telescopes = catalog.loc[catalog["Longitude"].isna(), "Code"]

        radar = self._radar_table
        counts = pd.concat(
            [
                self._table.get("observatory", pd.Series(dtype=str)),
                radar["transmitter"],
                radar["receiver"].where(radar["receiver"] != radar["transmitter"]),
            ]
        ).dropna().astype(str).value_counts().rename("count")

        if only_in_batch:
            table = (
                counts.rename_axis("Code")
                .reset_index()
                .merge(catalog[["Code", "Name"]], on="Code", how="left")
            )
            table["Name"] = table["Name"].fillna(table["Code"])
        else:
            table = catalog[["Code", "Name"]].copy()
            table["count"] = table["Code"].map(counts).fillna(0).astype(int)

        if only_space_telescopes:
            table = table.loc[table["Code"].isin(space_telescopes)]
        if exclude_space_telescopes:
            table = table.loc[~table["Code"].isin(space_telescopes)]

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

        Canonical radar table.

        :type: pd.DataFrame
        """
        return self._radar_table

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
        optical = self._table
        radar = self._radar_table

        self._space_telescopes = []
        self._bands = list(optical["band"].unique()) if "band" in optical.columns else []
        self._observatories = pd.unique(
            pd.concat(
                [
                    optical.get("observatory", pd.Series(dtype=str)),
                    radar["transmitter"],
                    radar["receiver"],
                ]
            )
            .dropna()
            .astype(str)
        ).tolist()
        optical_targets = (
            list(_resolve_optical_target_names(optical).values()) if not optical.empty else []
        )
        self._MPC_codes = list(
            dict.fromkeys(optical_targets + radar["target_body"].astype(str).tolist())
        )
        self._size = len(optical) + len(radar)

        epoch_columns = [
            column
            for column in [
                optical.get("epoch_seconds_UTC", pd.Series(dtype=float)),
                radar["epoch_seconds_UTC"],
            ]
            if len(column)
        ]
        epochs = pd.concat(epoch_columns) if epoch_columns else pd.Series(dtype=float)
        self._epoch_start, self._epoch_end = (
            (float(epochs.min()), float(epochs.max())) if len(epochs) else (0.0, 0.0)
        )

    def _add_custom_name_column(self, table: pd.DataFrame, custom_name) -> pd.DataFrame:
        """Return a copy of an optical table with a ``custom_name`` column."""
        augmented_table = table.copy()
        if custom_name is not None or "custom_name" not in augmented_table.columns:
            augmented_table["custom_name"] = [custom_name] * len(augmented_table)
        return augmented_table

    @staticmethod
    def _radar_observable_types(observation_types: list[str] | None):
        if observation_types is None:
            return None
        radar_types = [RANGE_OBSERVABLE, DOPPLER_OBSERVABLE]
        if "R" in observation_types:
            return radar_types
        return [value for value in observation_types if value in radar_types]

    def _add_parsed_mpc80_table(self, parsed_table, custom_name: str | None = None) -> None:
        if len(parsed_table) > 0:
            optical_table = create_augmented_optical_table(
                parsed_table,
                in_degrees=False,
                frame="J2000",
                custom_name=custom_name,
            )
            self._table = pd.concat([self._table, optical_table], ignore_index=True, sort=False)

        radar_table = radar_data_from_table(parsed_table)
        if custom_name is not None:
            radar_table = radar_table.assign(target_body=str(custom_name))
        if not radar_table.empty:
            self._radar_table = (
                radar_table.copy()
                if self._radar_table.empty
                else pd.concat([self._radar_table, radar_table], ignore_index=True)
            )

    @staticmethod
    def _fetch_mpc80_records(code: str | int, id_type: str | None = None) -> list[str]:
        query_kwargs = {"get_mpcformat": True}
        if id_type is not None:
            query_kwargs["id_type"] = id_type
        return [str(row) for row in MPC.get_observations(code, **query_kwargs)["obs"]]

    def _spacecraft_codes_with_positions(self) -> set[str]:
        """Return observatory codes with complete MPC spacecraft positions."""
        position_columns = list(self._SPACECRAFT_POSITION_COLUMNS)
        if "observatory" not in self._table or not set(position_columns).issubset(
            self._table.columns
        ):
            return set()

        return set(
            self._table.loc[
                self._table[position_columns].notna().all(axis=1),
                "observatory",
            ]
            .astype(str)
            .str.strip()
            .str.zfill(3)
        )

    def _resolve_spacecraft_observatory_code(self, satellite_name: str) -> str:
        """Resolve a spacecraft observatory name or code present in this batch."""
        satellite_name = str(satellite_name).strip()
        candidate_code = satellite_name.zfill(3)
        available_codes = self._spacecraft_codes_with_positions()

        if candidate_code in available_codes:
            return candidate_code

        observatory_catalog = MPC.get_observatory_codes().to_pandas()
        codes = observatory_catalog["Code"].astype(str).str.strip().str.zfill(3)
        names = observatory_catalog["Name"].astype(str).str.strip().str.lower()
        matching_rows = observatory_catalog.loc[
            (codes == candidate_code) | (names == satellite_name.lower())
        ]
        if not matching_rows.empty:
            resolved_code = str(matching_rows.iloc[0]["Code"]).strip().zfill(3)
            if resolved_code in available_codes:
                return resolved_code

        raise ValueError(
            f"Satellite '{satellite_name}' has no spacecraft positions in this batch. "
            f"Available observatory codes are: {sorted(available_codes)}."
        )

    def get_satellite_state_history(
        self,
        satellite_name: str,
        add_finite_difference_velocity: bool = True,
    ) -> dict[float, np.ndarray]:
        """Return the MPC-derived state history for a satellite observatory.

        The keys are TDB seconds since J2000. Each value is a Cartesian state
        in the Earth-centred, equatorial J2000 frame. MPC records provide only
        positions, so velocities are estimated with finite differences unless
        ``add_finite_difference_velocity`` is ``False``.

        Parameters
        ----------
        satellite_name : str
            MPC observatory code, such as ``"C51"``, or its exact name in the
            MPC observatory catalogue.
        add_finite_difference_velocity : bool, default True
            Whether to estimate velocities from the tabulated positions. If
            ``False``, zero velocities are appended.

        Returns
        -------
        dict[float, numpy.ndarray]
            TDB epochs mapped to six-element Cartesian states.
        """
        spacecraft_code = self._resolve_spacecraft_observatory_code(satellite_name)
        position_columns = list(self._SPACECRAFT_POSITION_COLUMNS)
        missing_columns = [
            column for column in position_columns if column not in self._table.columns
        ]
        if missing_columns:
            raise ValueError(
                "This batch does not contain MPC spacecraft-position columns. "
                "Load the observations from raw MPC 80-column records."
            )

        normalized_codes = self._table["observatory"].astype(str).str.strip().str.zfill(3)
        epoch_column = (
            "epoch_seconds_TDB"
            if "epoch_seconds_TDB" in self._table.columns
            else "epoch_seconds_UTC"
        )
        spacecraft_table = self._table.loc[
            (normalized_codes == spacecraft_code)
            & self._table[position_columns].notna().all(axis=1),
            [epoch_column, *position_columns],
        ]
        if spacecraft_table.empty:
            raise ValueError(
                f"No MPC spacecraft positions are available for satellite "
                f"'{satellite_name}' (observatory code '{spacecraft_code}')."
            )

        spacecraft_table = (
            spacecraft_table.groupby(epoch_column, as_index=True)[position_columns]
            .mean()
            .sort_index()
        )
        epochs = spacecraft_table.index.to_numpy(dtype=float)
        if epoch_column == "epoch_seconds_UTC":
            converter = time_representation.default_time_scale_converter()
            epochs = np.asarray(
                [
                    converter.convert_time(
                        input_scale=time_representation.utc_scale,
                        output_scale=time_representation.tdb_scale,
                        input_value=float(epoch),
                    )
                    for epoch in epochs
                ]
            )

        positions = spacecraft_table.to_numpy(dtype=float)
        if add_finite_difference_velocity and len(epochs) > 1:
            velocities = np.gradient(
                positions,
                epochs,
                axis=0,
                edge_order=2 if len(epochs) > 2 else 1,
            )
        else:
            velocities = np.zeros_like(positions)

        states = np.hstack((positions, velocities))
        return {float(epoch): state for epoch, state in zip(epochs, states)}

    @staticmethod
    def satellite_state_history_to_ephemeris_settings(
        state_history: dict[float, np.ndarray],
        frame_origin: str = "Earth",
        frame_orientation: str = "J2000",
    ):
        """Convert a TDB satellite state history to tabulated ephemeris settings."""
        return environment_setup.ephemeris.tabulated(
            body_state_history=state_history,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
        )

    def get_satellite_ephemeris_settings(
        self,
        satellite_name: str,
        frame_origin: str = "Earth",
        frame_orientation: str = "J2000",
        add_finite_difference_velocity: bool = True,
    ):
        """Return tabulated ephemeris settings from MPC satellite positions."""
        return self.satellite_state_history_to_ephemeris_settings(
            self.get_satellite_state_history(
                satellite_name,
                add_finite_difference_velocity=add_finite_difference_velocity,
            ),
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
        )

    ###########################################################################################
    # MPC Astroquery Data Retrieval: get_observations
    ###########################################################################################
    def get_observations(
        self,
        MPCcodes: list[str | int],
        id_types: list[str | None] | None = None,
        drop_misc_observations: bool = True,
        custom_name: str | None = None,
        use_mpc80_format: bool = False,
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
        use_mpc80_format : bool, default False
            Retrieve raw MPC 80-column records instead of Astroquery's parsed
            table. Required for radar and space-based observations.
            Drop-misc filtering is ignored in this mode.

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

            if use_mpc80_format:
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
