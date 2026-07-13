import pandas as pd
import numpy as np
from astroquery.mpc import MPC
import astropy
import copy
import importlib
from tudatpy.astro import time_representation
from tudatpy.data_access.tracking import TrackingData
from tudatpy.data_access.tracking import optical_utilities

_eighty_column = importlib.import_module("tudatpy.data_access.tracking.80_column")
parse_80cols_file = _eighty_column.parse_80cols_file
unpackers = importlib.import_module("tudatpy.data_access.tracking.80_column.unpackers")
OBS_TYPES_TO_DROP = unpackers.OBS_TYPES_TO_DROP


REQUIRED_OPTICAL_COLUMNS = ["number", "epoch", "RA", "DEC", "band", "observatory"]
ANCILLARY_STRING_COLUMNS = ["band", "catalog", "note2", "custom_name", "mag", "discovery"]


def standardize_optical_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure IDs are strings and observatory codes are zero-padded."""
    df = df.copy()

    if "observatory" in df.columns:
        df["observatory"] = df["observatory"].astype(str).str.strip().str.zfill(3)

    if "number" in df.columns:
        df["number"] = df["number"].astype(str).str.strip()

    return df


def validate_optical_table(
    table: astropy.table.QTable | astropy.table.Table | pd.DataFrame, frame: str
) -> None:
    """Validate the common optical astrometry table input."""
    if frame != "J2000":
        raise NotImplementedError("Only observations in J2000 are supported currently")

    if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
        colnames = table.colnames
        nrows = len(table)
    elif isinstance(table, pd.DataFrame):
        colnames = table.columns
        nrows = len(table)
    else:
        raise TypeError(f"Unsupported table type: {type(table).__name__}")

    if not set(REQUIRED_OPTICAL_COLUMNS).issubset(set(colnames)):
        raise ValueError(
            f"Table must include a set of mandatory columns: {REQUIRED_OPTICAL_COLUMNS}"
        )

    if nrows == 0:
        raise ValueError("Table contains zero rows: no valid observations were parsed.")


def create_augmented_optical_table(
    table: astropy.table.QTable | astropy.table.Table | pd.DataFrame,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
) -> pd.DataFrame:
    """Create the shared augmented optical table used by all optical sources."""
    validate_optical_table(table, frame)

    if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
        augmented_table = table.to_pandas()
    else:
        augmented_table = table.copy()

    augmented_table = standardize_optical_dataframe(augmented_table)

    if custom_name is not None or "custom_name" not in augmented_table.columns:
        augmented_table["custom_name"] = [custom_name] * len(augmented_table)

    if in_degrees:
        augmented_table = augmented_table.assign(
            RA=lambda x: (np.radians(x.RA) + np.pi) % (2 * np.pi) - np.pi,
            DEC=lambda x: np.radians(x.DEC),
        )

    if "mag" not in augmented_table.columns and "magnitude" in augmented_table.columns:
        augmented_table["mag"] = augmented_table["magnitude"]

    for column in ["catalog", "note2", "mag", "discovery"]:
        if column not in augmented_table.columns:
            augmented_table[column] = None

    augmented_table["epoch_seconds_UTC"] = [
        time_representation.julian_day_to_seconds_since_epoch(jd)
        for jd in list(augmented_table["epoch"])
    ]

    return augmented_table


def _datetime_to_utc_seconds(epoch) -> float:
    if hasattr(epoch, "to_epoch"):
        return float(epoch.to_epoch())
    if hasattr(epoch, "to_float"):
        return float(epoch.to_float())
    if hasattr(epoch, "year") and hasattr(epoch, "month") and hasattr(epoch, "day"):
        return float(time_representation.DateTime.from_python_datetime(epoch).to_epoch())
    return float(epoch)


def filter_augmented_optical_table(
    table: pd.DataFrame,
    epoch_start=None,
    epoch_end=None,
    observatories: list[str | int] | None = None,
    observatories_exclude: list[str | int] | None = None,
) -> pd.DataFrame:
    """Filter an augmented optical astrometry table.

    Epoch filters are interpreted as UTC seconds since J2000, or as Python
    datetimes that are converted to UTC seconds since J2000.
    """

    if "epoch_seconds_UTC" not in table.columns:
        raise ValueError("Expected an augmented optical table with an epoch_seconds_UTC column.")

    filtered = table.copy()
    if epoch_start is not None:
        filtered = filtered.loc[
            filtered["epoch_seconds_UTC"] >= _datetime_to_utc_seconds(epoch_start)
        ]
    if epoch_end is not None:
        filtered = filtered.loc[
            filtered["epoch_seconds_UTC"] <= _datetime_to_utc_seconds(epoch_end)
        ]

    if observatories is not None:
        included = {str(observatory).strip().zfill(3) for observatory in observatories}
        filtered = filtered.loc[filtered["observatory"].isin(included)]
    if observatories_exclude is not None:
        excluded = {str(observatory).strip().zfill(3) for observatory in observatories_exclude}
        filtered = filtered.loc[~filtered["observatory"].isin(excluded)]

    return filtered.reset_index(drop=True)


def optical_table_to_tracking_data(
    table: pd.DataFrame,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Convert an augmented optical table to TrackingData and supplementary data lists."""
    table = create_augmented_optical_table(table, in_degrees=False)
    if add_weights and not weighing_scheme:
        weighing_scheme = "VFCC17"

    if add_star_catalog_corrections:
        RA_corr, DEC_corr = optical_utilities.get_biases_EFCC18(mpc_table=table)
        table = table.assign(_RA_corr=RA_corr, _DEC_corr=DEC_corr)

    if add_ancillary_data:
        table = table.copy()
        for column in ANCILLARY_STRING_COLUMNS:
            table[column] = table[column].fillna("").astype(str)

    tracking_data_objects = []
    for (target, observatory), group in table.groupby(["number", "observatory"]):
        observable_type, reference_link_end_type = "AngularPosition", "receiver"

        link_ends = [
            ((str(target), ""), "transmitter"),
            (("Earth", str(observatory)), reference_link_end_type),
        ]
        observations = [np.array([ra, dec]) for ra, dec in zip(group["RA"], group["DEC"])]

        tracking_data_object = TrackingData(
            observable_type=observable_type,
            link_ends=link_ends,
            observations=observations,
            epochs=group["epoch_seconds_UTC"],
            reference_link_end=reference_link_end_type,
            time_scale="UTC",
            weighing_scheme=weighing_scheme,
        )

        if add_star_catalog_corrections:
            corrections_list = [
                np.array([ra_c, dec_c])
                for ra_c, dec_c in zip(group["_RA_corr"], group["_DEC_corr"])
            ]
            tracking_data_object.set_observation_corrections(corrections_list)

        if add_ancillary_data:
            for column in ANCILLARY_STRING_COLUMNS:
                tracking_data_object.add_ancillary_settings(
                    column,
                    group[column].fillna("").astype(str).tolist(),
                )

        tracking_data_objects.append(tracking_data_object)

    return tracking_data_objects, list()


def read_pandas_optical_data(
    table: pd.DataFrame,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Read optical astrometry from a pandas table into TrackingData objects."""
    augmented_table = create_augmented_optical_table(table, in_degrees, frame, custom_name)
    return optical_table_to_tracking_data(
        augmented_table,
        add_weights,
        add_star_catalog_corrections,
        add_ancillary_data,
        weighing_scheme,
    )


def read_astropy_optical_data(
    table: astropy.table.QTable | astropy.table.Table,
    in_degrees: bool = True,
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Read optical astrometry from an astropy table into TrackingData objects."""
    augmented_table = create_augmented_optical_table(table, in_degrees, frame, custom_name)
    return optical_table_to_tracking_data(
        augmented_table,
        add_weights,
        add_star_catalog_corrections,
        add_ancillary_data,
        weighing_scheme,
    )


def read_80_column_data(
    filename: str | list[str],
    frame: str = "J2000",
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Read MPC 80-column files into TrackingData objects."""
    return read_astropy_optical_data(
        parse_80cols_file(filename),
        in_degrees=False,
        frame=frame,
        custom_name=custom_name,
        add_weights=add_weights,
        add_star_catalog_corrections=add_star_catalog_corrections,
        add_ancillary_data=add_ancillary_data,
        weighing_scheme=weighing_scheme,
    )


def read_mpc_data(
    MPCcodes: list[str | int],
    id_types: list[str | None] | None = None,
    drop_misc_observations: bool = True,
    custom_name: str | None = None,
    add_weights: bool | None = False,
    add_star_catalog_corrections: bool | None = False,
    add_ancillary_data: bool | None = False,
    weighing_scheme: str = "",
):
    """Retrieve MPC observations and return TrackingData and supplementary data lists."""
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
        weighing_scheme,
    )


class BatchMPC:
    """This class provides an interface between observations
    in the Minor Planet Centre database and Tudat.

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

        # for manual additions of table (from_pandas, from_astropy)
        self._req_cols = ["number", "epoch", "RA", "DEC", "band", "observatory"]

    def __copy__(self):
        new = BatchMPC()

        new._table = copy.deepcopy(self._table)
        new._refresh_metadata()

        return new

    def copy(self) -> "BatchMPC":
        """Create a copy of the batch, equivalent to copy.copy(batchMPC())

        Returns
        -------
        BatchMPC
            Copy of batch.
        """
        return copy.copy(self)

    # getters to make everything read-only
    @property
    def table(self) -> pd.DataFrame:
        """Pandas dataframe with observation data."""
        return self._table

    @property
    def observatories(self) -> list[str]:
        """List of observatories in batch."""
        return self._observatories

    @property
    def space_telescopes(self) -> list[str]:
        """List of satellite_observatories in batch."""
        return self._space_telescopes

    @property
    def MPC_objects(self) -> list[str]:
        """List of MPC objects."""
        return self._MPC_codes

    @property
    def size(self) -> int:
        """Number of observations in batch."""
        return self._size

    @property
    def bands(self) -> list[str]:
        """List of bands in batch."""
        return self._bands

    @property
    def epoch_start(self) -> float:
        """Epoch of oldest observation in batch in seconds since J2000 TDB."""
        return self._epoch_start

    @property
    def epoch_end(self) -> float:
        """Epoch of latest observation in batch in seconds since J2000 TDB."""
        return self._epoch_end

    @property
    def bodies_created(self) -> dict:
        """Dictionary with the bodies created by to_tudat and details."""
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

    def _standardize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Internal helper to ensure IDs are strings and observatories are zero-padded."""
        return standardize_optical_dataframe(df)

    def _add_custom_name_column(self, table: pd.DataFrame, custom_name) -> pd.DataFrame:
        augmented_table = table.copy()
        if custom_name is not None or "custom_name" not in augmented_table.columns:
            augmented_table["custom_name"] = [custom_name] * len(augmented_table)
        return augmented_table

    def _add_table(self, table: pd.DataFrame, custom_name: str | None, in_degrees: bool = True):
        """Internal. Formats a table of observations, used in from_astropy and in from_pandas."""
        obs = create_augmented_optical_table(
            table=table, in_degrees=in_degrees, frame="J2000", custom_name=custom_name
        )

        # convert object mpc code to string
        self._table = pd.concat([self._table, obs])
        self._refresh_metadata()

    def _validate_table(
        self, table: astropy.table.QTable | astropy.table.Table | pd.DataFrame, frame: str
    ) -> None:
        """Internal helper to validate the frame and required columns of a table.

        Parameters
        ----------
        table : astropy.table.QTable | astropy.table.Table | pd.DataFrame
            The table to validate.
        frame : str
            The reference frame to check.
        """
        validate_optical_table(table, frame)

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
        Removes duplicate and irrelevant observation data by default (see `drop_misc_observations`).

        Parameters
        ----------
        MPCcodes : list[str | int]
            List of integer or str MPC object codes for minor planets or comets.
        id_types : list[str | None] | None, default None
            A list of identification types ('asteroid_number', 'comet_number', 'comet_designation') corresponding to each MPCcode.
            If an element is None, the type is considered unknown. If the entire list is None,
            all types are considered unknown.
        drop_misc_observations : bool, default True
            Drops observations made by method: radar and offset (natural satellites).
            Drops observations made by roaming observers.
            Drops duplicate listings to denote first observation.

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

            obs = self._standardize_dataframe(obs)
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

    ################################################################################
    # Allow external observations to be added: from_file
    ################################################################################
    def from_file(
        self,
        filename: str,
        in_degrees: bool = False,
        frame: str = "J2000",
        custom_name: str | None = None,
    ) -> None:
        """
        Loads observations from a local MPC 80-column text file.

        This method serves as a high-level convenience function that orchestrates the
        parsing of a raw 80-column file and loading the data into the batch. It uses
        the `tudatpy.data_access.tracking.80_column.parse_80cols_file` function internally.

        The parser returns an Astropy Table with RA/DEC values in radians, so this
        method subsequently calls `from_astropy` to ingest the data. The `in_degrees`
        parameter should therefore be `False` (the default).

        Note
        ----
        If you wish to perform intermediate operations on the parsed data using
        pandas, you can call the parser directly and then use the `from_pandas`
        method:

        .. code-block:: python

            parse_80cols_file = importlib.import_module(
                "tudatpy.data_access.tracking.80_column"
            ).parse_80cols_file

            # 1. Parse the file to an Astropy Table
            astropy_table = parse_80cols_file("my_obs.txt")

            # 2. Convert to a pandas DataFrame for manipulation
            pandas_df = astropy_table.to_pandas()
            # ... perform custom pandas operations on pandas_df ...

            # 3. Load the processed DataFrame into the batch
            batch = BatchMPC()
            batch.from_pandas(pandas_df, in_degrees=False)


        Parameters
        ----------
        filename : str
            The path to the MPC 80-column formatted text file.
        in_degrees : bool, optional
            Specifies the unit of RA/DEC in the file. Since the internal parser
            handles the conversion to radians, this should be left as `False`.
            Defaults to False.
        frame : str, optional
            The reference frame of the observations. Currently, only "J2000" is
            supported. Defaults to "J2000".
        """
        # Use the dedicated parser submodule to parse the external file.
        # This function returns an astropy.Table with RA/DEC in radians.
        astropy_table = parse_80cols_file(filename)

        # Use the from_astropy method to validate and add the data.
        # in_degrees is False because the parser has already converted to radians.
        self.from_astropy(
            astropy_table, in_degrees=in_degrees, frame=frame, custom_name=custom_name
        )

    #########################################################
    # FROM ASTROPY OBJECT INTO Tudat BatchMPC: from_astropy
    #########################################################
    def from_astropy(
        self,
        table: astropy.table.QTable | astropy.table.Table,
        in_degrees: bool = True,
        frame: str = "J2000",
        custom_name: str | None = None,
    ) -> None:
        """Loads observations from an Astropy Table into the BatchMPC object.

        This method provides a convenient way to import observation data that
        is stored in an Astropy Table or QTable. It serves
        as a wrapper that validates the input before converting it to a pandas DataFrame
        for internal processing.

        Args:
            table (astropy.table.Table):
                The Astropy Table or QTable containing the observation data. It must
                include the following columns: 'number', 'epoch', 'RA', 'DEC',
                'band', 'observatory'.
            in_degrees (bool, optional):
                If True, 'RA' and 'DEC' columns are assumed to be in degrees.
                If False, they are assumed to be in radians. Defaults to True.
            frame (str, optional):
                The reference frame of the observations. Currently, only "J2000" is
                supported. Defaults to "J2000".

        Raises:
            ValueError: If the input `table` is not an Astropy Table/QTable or if it
                is missing any of the required columns.
            NotImplementedError: If a `frame` other than 'J2000' is provided.
        """
        if not isinstance(table, (astropy.table.QTable, astropy.table.Table)):
            raise ValueError("Table must be of type astropy.table.QTable or astropy.table.Table")

        self._validate_table(table, frame)
        self._add_table(table=table.to_pandas(), in_degrees=in_degrees, custom_name=custom_name)

    ######################################################
    # FROM PANDAS OBJECT INTO Tudat BatchMPC: from_pandas
    ######################################################
    def from_pandas(
        self,
        table: pd.DataFrame,
        in_degrees: bool = True,
        frame: str = "J2000",
        custom_name: str | None = None,
    ) -> None:
        """
        Loads observations from a pandas DataFrame into the BatchMPC object.

        The DataFrame must contain the following columns:
        - 'number': The MPC object code (e.g., '433' for Eros).
        - 'epoch': The observation epoch in Julian Day format.
        - 'RA': Right Ascension.
        - 'DEC': Declination.
        - 'band': The observation band.
        - 'observatory': The MPC observatory code.

        Parameters
        ----------
        table : pd.DataFrame
            The pandas DataFrame containing the observation data.
        in_degrees : bool, optional
            If True, RA and DEC columns in the DataFrame are assumed to be in degrees.
            If False, they are assumed to be in radians. Defaults to True.
        frame : str, optional
            The reference frame of the observations. Currently, only "J2000" is supported.
            Defaults to "J2000".
        """
        if not isinstance(table, pd.DataFrame):
            raise ValueError("Table must be of type pandas.DataFrame")

        self._validate_table(table, frame)
        self._add_table(table=table, in_degrees=in_degrees, custom_name=custom_name)

    def to_tracking_dataset(
        self,
        add_weights: bool | None = False,
        add_star_catalog_corrections: bool | None = False,
        add_ancillary_data: bool | None = False,
        weighing_scheme: str = "",
    ):
        return optical_table_to_tracking_data(
            self._table,
            add_weights=add_weights,
            add_star_catalog_corrections=add_star_catalog_corrections,
            add_ancillary_data=add_ancillary_data,
            weighing_scheme=weighing_scheme,
        )
