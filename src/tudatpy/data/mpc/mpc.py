from tudatpy.dynamics import environment_setup
from astropy.table import Table
import pandas as pd
import numpy as np
import datetime
from astroquery.mpc import MPC
import astropy
import copy
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.data.mpc.parser_80col import parse_80cols_file
from tudatpy.data.mpc.parser_80col import unpackers

# do not remove this line, even if it looks liek an unused import line
from tudatpy.data.mpc.parser_80col.unpackers import OBS_TYPES_TO_DROP


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
    >>> batch.filter(epoch_start=datetime.datetime(2016, 1, 1))

    Transform to Tudat format:
    >>> ...
    >>> bodies = environment_setup.create_system_of_bodies(body_settings)
    >>> observation_collection = batch.to_tudat(bodies=bodies, included_satellites=None)

    """

    def __init__(self) -> None:
        """Create an empty MPC batch."""
        self._table: pd.DataFrame = pd.DataFrame()
        self._observatories: list[str] = []
        self._space_telescopes: list[str] = []
        self._bands: list[str] = []
        self._MPC_codes: list[str] = []
        self._size: int = 0

        self._observatory_info: pd.DataFrame | None = None
        self._MPC_space_telescopes: list[str] = []

        self._get_station_info()

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
        self._space_telescopes = [x for x in self._observatories if x in self._MPC_space_telescopes]
        if "bands" in self._table.columns:
            self._bands = list(self._table.band.unique())

        # if user gives custom name, set that as body name, else MPC code
        if "custom_name" in self._table.columns and self._table["custom_name"].notna().any():
            self._MPC_codes = list(self._table["custom_name"].unique())
        else:
            self._MPC_codes = list(self._table.number.unique())
        self._size = len(self._table)

        self._epoch_start = self._table.epoch_seconds_TDB.min()
        self._epoch_end = self._table.epoch_seconds_TDB.max()

    def _get_station_info(self) -> None:
        """Internal. Retrieve data on MPC listed observatories."""
        try:
            temp = MPC.get_observatory_codes().to_pandas()
            temp["Code"] = temp["Code"].astype(str).str.strip().str.zfill(3)
            # This query checks if Longitude is Nan: non-terrestrial telescopes
            sats = list(temp.query("Longitude != Longitude").Code.values)
            self._observatory_info = temp
            self._MPC_space_telescopes = sats
        except Exception as e:
            print("An error occured while retrieving observatory data")
            print(e)

    def _standardize_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Internal helper to ensure IDs are strings and observatories are zero-padded."""
        # Work on a copy to avoid SettingWithCopyWarnings
        df = df.copy()

        # Standardize Observatory Codes (e.g., 673 -> "673", 89 -> "089", C51 -> "C51")
        if "observatory" in df.columns:
            df["observatory"] = df["observatory"].astype(str).str.strip().str.zfill(3)

        # Standardize MPC Numbers (e.g., 433 -> "433")
        if "number" in df.columns:
            df["number"] = df["number"].astype(str).str.strip()

        return df

    def _add_time_columns(self, table: pd.DataFrame) -> pd.DataFrame:
        """
        Internal helper to add standardized time columns to an observation table.

        This function takes a DataFrame with an 'epoch' column (in Julian Days, UTC)
        and adds 'epoch_seconds_TDB' and 'epoch_seconds_UTC' columns.

        Parameters
        ----------
        table : pd.DataFrame
            The input DataFrame, which must contain an 'epoch' column.

        Returns
        -------
        pd.DataFrame
            The DataFrame with the new time columns added.
        """

        # Create DateTime objects from the Julian Day 'epoch' column
        dt_objects = [DateTime.from_julian_day(jd) for jd in table["epoch"]]

        # Get the default time scale converter
        time_scale_converter = time_representation.default_time_scale_converter()

        augmented_table = table.copy()
        # Add 'epoch_seconds_UTC' column by converting DateTime Objects to epoch
        augmented_table["epoch_seconds_UTC"] = [dt_obj.epoch() for dt_obj in dt_objects]

        # Add 'epoch_seconds_TDB' column by converting from UTC to TDB
        augmented_table["epoch_seconds_TDB"] = [
            time_scale_converter.convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=t_utc,
            )
            for t_utc in augmented_table["epoch_seconds_UTC"]
        ]

        return augmented_table

    def _add_custom_name_column(self, table: pd.DataFrame, custom_name) -> pd.DataFrame:
        augmented_table = table.copy()
        augmented_table["custom_name"] = [custom_name] * len(augmented_table)
        return augmented_table

    ################################################################################
    # DATA RETRIEVAL OPTION 1): from_file: allows external observations to be added
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
        the `tudatpy.data.mpc.parser_80col.parse_80cols_file` function internally.

        The parser returns an Astropy Table with RA/DEC values in radians, so this
        method subsequently calls `from_astropy` to ingest the data. The `in_degrees`
        parameter should therefore be `False` (the default).

        Note
        ----
        If you wish to perform intermediate operations on the parsed data using
        pandas, you can call the parser directly and then use the `from_pandas`
        method:

        .. code-block:: python

            from tudatpy.data.mpc.parser_80col import parse_80cols_file

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

    ###########################################################################################
    # DATA RETRIEVAL OPTION 2): get_observations: retrieves data from mpc through astroquery.
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

            # convert JD to J2000 and UTC, convert deg to rad
            obs = self._add_time_columns(obs)
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
            self._table = pd.concat([self._table, obs])

        self._refresh_metadata()

    def _add_table(self, table: pd.DataFrame, custom_name: str | None, in_degrees: bool = True):
        """Internal. Formats a manually entered table of observations, used in from_astropy and in from_pandas."""
        obs = table
        obs = self._add_time_columns(obs)
        obs = self._add_custom_name_column(obs, custom_name)
        if in_degrees:
            obs = obs.assign(RA=lambda x: (np.radians(x.RA) + np.pi) % (2 * np.pi) - np.pi).assign(
                DEC=lambda x: np.radians(x.DEC)
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
        if frame != "J2000":
            raise NotImplementedError("Only observations in J2000 are supported currently")

        # Get column names depending on table type
        if isinstance(table, (astropy.table.QTable, astropy.table.Table)):
            colnames = table.colnames
            nrows = len(table)
        elif isinstance(table, pd.DataFrame):
            colnames = table.columns
            nrows = len(table)
        else:
            raise TypeError(f"Unsupported table type: {type(table).__name__}")

        if not set(self._req_cols).issubset(set(colnames)):
            raise ValueError(f"Table must include a set of mandatory columns: {self._req_cols}")

        if nrows == 0:
            raise ValueError("Table contains zero rows: no valid observations were parsed.")

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

        This method provides a convenient way to import observation data that has been
        processed or filtered and is stored in an Astropy Table or QTable. It serves
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
