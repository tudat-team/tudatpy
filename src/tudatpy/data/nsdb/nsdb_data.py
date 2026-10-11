from io import StringIO

import bs4
import pandas as pd
import pyparsing as pp
import requests

import tudatpy.dynamics.environment as env
import tudatpy.estimation.observable_models_setup as obs_model_setup
import tudatpy.estimation.observations as obs

# Pyparsing: newlines are significant in NSDB data, only skip spaces/tabs
pp.ParserElement.set_default_whitespace_chars(" \t")

CATALOGUES = {
    "Mars": "https://nsdb.imcce.fr/obspos/OBS_COLL/M/",
    "Jupiter": "https://nsdb.imcce.fr/obspos/OBS_COLL/J/",
    "Saturn": "https://nsdb.imcce.fr/obspos/OBS_COLL/S/",
    "Uranus": "https://nsdb.imcce.fr/obspos/OBS_COLL/U/",
    "Neptune": "https://nsdb.imcce.fr/obspos/OBS_COLL/N/",
    "Pluto": "https://nsdb.imcce.fr/obspos/OBS_COLL/P/",
    "Asteroids": "https://nsdb.imcce.fr/obspos/OBS_COLL/abin/",
}

_CATALOGUE_ABBREVIATIONS = {
    "M": "Mars",
    "J": "Jupiter",
    "S": "Saturn",
    "U": "Uranus",
    "N": "Neptune",
    "P": "Pluto",
    "abin": "Asteroids",
}

DEFAULT_STATION_INFO_SOURCE = ""  # Observatories.txt


class NSDBQuery:
    """
    Class to query the NSDB database.
    """

    def __init__(self):
        self.parsed_datasets = {}  # dataset name -> {contents, informations, format}
        self.data_files = {}  # dataset name -> data file content

    def get_dataset(self, datasets: list[str]):
        self._get_metadata(datasets, raise_on_missing=True)
        return self._get_observations(datasets, self.parsed_datasets, raise_on_missing=True)

    def _get_observations(self, datasets: list[str], metadata: pd.DataFrame, raise_on_missing: bool = False) -> pd.DataFrame:
        missing_datasets = []
        for dataset in datasets:
            # Identify primary XYzzzz where X is primary
            primary = dataset[0].upper()
            target = _CATALOGUE_ABBREVIATIONS.get(primary, None)
            if target is None:
                raise ValueError(f"Unknown primary {primary} in dataset {dataset}.")
            catalogue_url = CATALOGUES[target]
            # Get the data file for the dataset
            data_file = requests.get(f"{catalogue_url}/{dataset}.txt")
            try:
                data_file.raise_for_status()  # Raise an error if the request failed
            except requests.exceptions.HTTPError as e:
                if data_file.status_code == 404:
                    missing_datasets.append(dataset)
                    continue
                else:
                    raise e
            # parse data into a pandas DataFrame, where columns are separated by some number of whitespaces, and rows by newline
            df = pd.read_csv(
                StringIO(data_file.text),
                sep=r"\s+",
                header=None,
                comment="#",
            )

            df.columns = metadata[dataset]["format"].values()
            self.data_files[dataset] = df
            # TODO: Preprocess the dataframe for units, time, frame, etc.
            self._process_data(dataset, metadata[dataset])
            
        return self.data_files

    def _process_data(self, dataset: str, metadata: pd.DataFrame):
        """
        Process the data file for a given dataset into a pandas DataFrame.

        This method applies:
         - Time scale conversion
         - Epoch of equinox precession
         - Reference frame correction (Astrometric vs Apparent)
         - Reduction correction (Only when Reduction is known, e.g. UCAC2 -> ICRF)

        Args:
            dataset (str): The name of the dataset.
            metadata (pd.DataFrame): The metadata for the dataset.

        Returns:
            pd.DataFrame: A DataFrame containing the processed data.
        """
        # TODO: Implement data processing logic
        pass

    def _ensure_observatory(self, observatory_code: str, bodies: env.SystemOfBodies) -> bool:
        """
        Ensure that the observatory information is available for a given observatory code.

        Args:
            observatory_code (str): The code of the observatory.
            bodies (env.SystemOfBodies): The system of bodies containing the observatory information.
        Returns:
            bool: True if the observatory information is available, False otherwise.
        """
        try:
            bodies.get("Earth").get_ground_station(observatory_code)
        except RuntimeError:
            return False # Fails if the observatory is not found
        else:
            return True # Succeeds if the observatory is found

    def _get_metadata(
        self, datasets: list[str], raise_on_missing: bool = False
    ) -> pd.DataFrame:
        """
        Get the observations for the given datasets from the NSDB database.

        Returns:
            pd.DataFrame: A DataFrame containing the observations with columns:
                - 'dataset': The dataset name.
                - 'primary': The primary body (e.g., 'Neptune').
                - 'target': The target body (e.g., 'Triton').
                - 'type': The type of observation (e.g., 'absolute', 'relative').
                - 'observatory': The observatory code and name.
        """

        html_pages = []
        missing_datasets = []
        for dataset in datasets:
            # Identify primary XYzzzz where X is primary
            primary = dataset[0].upper()
            target = _CATALOGUE_ABBREVIATIONS.get(primary, None)
            if target is None:
                raise ValueError(f"Unknown primary {primary} in dataset {dataset}.")
            catalogue_url = CATALOGUES[target]

            # Get the HTML page for the dataset
            html = requests.get(f"{catalogue_url}/{dataset}.html")
            try:
                html.raise_for_status()  # Raise an error if the request failed
            # Catch 404
            except requests.exceptions.HTTPError as e:
                if html.status_code == 404:
                    missing_datasets.append(dataset)
                    continue
                else:
                    raise e
            html_pages.append(html)

            # Pass the HTML page to the parser
            self._parse_html(html, dataset)

            

        if missing_datasets and raise_on_missing:
            raise ValueError(
                f"The following datasets were not found in the NSDB database: {missing_datasets}"
            )
        return self.parsed_datasets

    def _parse_html(self, html, dataset_name: str = ""):
        """
        Parse an NSDB dataset HTML page into sections.

        The HTML contains a <pre> block with sections in the order:
            Contents., [Reference.], Informations., [Comments.], Format.
        Some sections (Reference., Comments.) may be absent.
        """
        soup = bs4.BeautifulSoup(html.content, "html.parser")
        description = soup.find("pre")
        lines = description.get_text().splitlines()

        # Known section headers in their canonical order
        section_headers = [
            "Contents.",
            "Reference.",
            "Informations.",
            "Comments.",
            "Format.",
        ]

        # Find which headers are present and their line indices
        header_positions = {}
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped in section_headers:
                header_positions[stripped] = i

        if not header_positions:
            raise ValueError("No section headers found in the HTML.")

        # Sort by line number to handle any ordering
        sorted_headers = sorted(header_positions.items(), key=lambda x: x[1])

        # Extract text between consecutive headers
        sections = {}
        for j, (header, start_line) in enumerate(sorted_headers):
            if j + 1 < len(sorted_headers):
                end_line = sorted_headers[j + 1][1]
            else:
                end_line = len(lines)
            content_lines = lines[start_line + 1 : end_line]
            sections[header] = "\n".join(content_lines).strip()

        contents_section = sections.get("Contents.", "")
        informations_section = sections.get("Informations.", "")
        format_section = sections.get("Format.", "")

        # Validate required sections
        if not contents_section:
            raise ValueError("Contents section is empty.")
        if not informations_section:
            raise ValueError("Informations section is empty.")
        if not format_section:
            raise ValueError("Format section is empty.")

        contents = self._parse_contents_section(contents_section)
        informations = self._parse_contents_section(informations_section)
        format_data = self._parse_format(format_section)

        # Store per-dataset parsed data
        self.parsed_datasets[dataset_name] = {
            "contents": contents,
            "informations": informations,
            "format": format_data,
        }

    @staticmethod
    def _parse_contents_section(text):
        """
        Parse a key-value section (Contents or Informations) into a dict using pyparsing.

        Handles multi-line values where continuation lines are indented
        and do not start a new key-value pair::

            key: first line
              continuation line
              another continuation line

        Single-line values are returned as strings.
        Multi-line values are returned as lists of strings, preserving line order.

        Raises ``ValueError`` on duplicate keys.
        """
        NL = pp.LineEnd().suppress()

        # A metadata key: word characters, spaces, dots (NOT newlines)
        key = pp.Regex(r"[\w .]+")

        # A field declaration: ``key: value``
        field_decl = (
            key("key")
            + pp.Suppress(pp.Literal(":"))
            + pp.Regex(r"[^\n]*")("value")
            + NL
        )

        section_headers = frozenset(
            [
                "Contents.",
                "Reference.",
                "Informations.",
                "Comments.",
                "Format.",
            ]
        )

        lines = text.splitlines()
        result = {}
        current_key = None
        current_values = []

        def _finalize():
            nonlocal current_key, current_values
            if current_key is not None:
                if len(current_values) == 1:
                    result[current_key] = current_values[0]
                else:
                    result[current_key] = current_values
                current_key = None
                current_values = []

        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped in section_headers:
                _finalize()
                break

            # Use pyparsing to recognise field declarations
            try:
                parsed = field_decl.parse_string(stripped + "\n")
                _finalize()
                k = parsed.key.strip()
                v = parsed.value.strip()
                if k in result:
                    raise ValueError(f"Duplicate key '{k}' found in section.")
                current_key = k
                current_values = [v]
            except pp.ParseException:
                # Continuation line – append to the current field
                if current_key is not None:
                    current_values.append(stripped)

        _finalize()
        return result

    @staticmethod
    def _parse_format(text):
        """
        Parse the Format section into a dict mapping column numbers to descriptions.

        Grammar: zero or more numbered items of the form ``N. description``.
        Separator lines and column-header lines are silently skipped.
        """
        NL = pp.LineEnd().suppress()
        number = pp.Regex(r"\d+")
        description = pp.Regex(r"[^\n]*")
        numbered_item = pp.Group(
            number("number")
            + pp.Suppress(pp.Literal("."))
            + description("description")
            + NL
        )

        result = {}
        for parsed, _, _ in numbered_item.scan_string(text):
            result[parsed[0].number] = parsed[0].description.strip()
        return result

    def table():
        pass

    def observatories():
        pass

    def objects():
        pass

    def size():
        pass

    def epoch_start():
        pass

    def epoch_end():
        pass

    def __add__(self, other):
        pass

    def _refresh_metadata(self):
        pass

    def _get_station_info(self):
        pass

    def set_weights(self, weights):
        pass

    def to_observation_collection(self):
        pass

    def to_tudat(self):
        pass

    def plot_prefit_residuals(self):
        pass

    def get_observation_models(self):
        pass

    @staticmethod
    def get_available_datasets(planet: str):
        """
        Get the available datasets for a given planet from the NSDB database.
        :param planet: Name of the planet (e.g., "Neptune").
        :return: A pandas DataFrame containing the available datasets.
        """

        # Get table of available datasets for the given planet
        planet = planet.capitalize()
        planet_page = requests.get(CATALOGUES[planet])
        planet_page.raise_for_status()  # Raise an error if the request failed
        soup = bs4.BeautifulSoup(planet_page.content, "html.parser")
        table = soup.find("table")
        if table is None:
            raise ValueError(
                f"No table found for planet {planet}. Ensure it is available on NSDB."
            )
        else:
            # Dump the table rows
            rows = table.find_all("tr")
            for row in rows:
                cols = row.find_all("td")
                cols = [ele.text.strip() for ele in cols]
                # print(cols)  # Print the columns of each row

        # Cast the table to a pandas DataFrame
        data = []
        for row in rows:  # Skip the header row
            cols = row.find_all("td")
            cols = [ele.text.strip() for ele in cols]
            data.append(cols)
        df = pd.DataFrame(data)

        # Drop empty columns and formatting rows
        df = df.drop(columns=[0])  # Drop ICO/Type column (Empty)
        df = df.drop(columns=[4])  # Drop Description column (Empty)
        df = df.drop([0, 1, 2, len(df) - 1])  # Drop formatting rows (0, 1, 2, -1)
        df.columns = ["Name", "Last modified", "Size"]

        dataset_names = df["Name"].apply(lambda x: str(x).split(".")[0]).unique()
        # print(f"Found {len(dataset_names)} unique datasets for planet {planet}.")

        # Build new frame with columns:
        # Name, descriptor file size, descriptor file last modified, data file size, data file last modified
        # Where txt is the data file and html is the descriptor file
        dataset_info = []
        for dataset in dataset_names:
            descriptor_file = df[df["Name"] == f"{dataset}.html"]
            data_file = df[df["Name"] == f"{dataset}.txt"]
            if not descriptor_file.empty and not data_file.empty:
                dataset_info.append(
                    {
                        "Name": dataset,
                        "HTML size": descriptor_file["Size"].values[0],
                        "HTML last modified": descriptor_file["Last modified"].values[
                            0
                        ],
                        "Data file size": data_file["Size"].values[0],
                        "Data file last modified": data_file["Last modified"].values[0],
                    }
                )

        dataset_df = pd.DataFrame(dataset_info)
        # print(dataset_df.head())
        return dataset_df

