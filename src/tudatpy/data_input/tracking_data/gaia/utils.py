from pathlib import Path
import gzip
import itertools
import tempfile
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pv
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from tqdm import tqdm
from tudatpy.resource_installer import (
    DEFAULT_CACHE,
    DEFAULT_DEST,
    _verify_md5,
    download_extra_file,
)

_DATA_RELEASE = "gfpr"  # Gaia Focused Product Release
_DEFAULT_GAIA_RESOURCES = DEFAULT_DEST / "gaia"


def download_gaia_archive(
    astrometry: bool = True,
    asteroid_orbits: bool = True,
    destination_path: Path | str | None = None,
    force: bool = False,
) -> None:
    """
    This function downloads and stores a local copy of the Gaia archives ( https://gea.esac.esa.int/archive/)

    This function downloads and stores a local copy of the Gaia Focused Product Release (FPR), containing the observations,
    observation metadata and orbit data for all asteroids. This allows for offline and very fast
    retrieval of the Gaia observations. The data is stored as Apache Parquet files, taking up roughly 6GB in total.
    Archives whose .parquet file already exists are skipped, unless ``force`` is set. The downloaded .csv.gz files
    are removed once the .parquet file has been written. While downloading and converting, roughly 20GB of free disk
    space is needed.

    Parameters
    ----------
    astrometry : bool, optional
        Downloads the observations and metadata (``sso_observation``). Stored as 'sso_observation.parquet' (~6GB for FPR)
    asteroid_orbits : bool, optional
        Downloads asteroid orbit and covariance data (``sso_source``). Stored as 'sso_source.parquet' (~60MB for FPR)
    destination_path : Path | str, optional
        Define path to which output .parquet files are written. If None, it is stored under
        ``~/.tudat/resource/gaia/``
    force : bool, optional
        If true, downloads the archives again and overwrites existing files.
    """
    if destination_path is None:
        destination_path = _DEFAULT_GAIA_RESOURCES
    destination_path = Path(destination_path).expanduser()

    destination_path.mkdir(parents=True, exist_ok=True)

    archives = []
    if astrometry:
        # Observation and metadata
        archives.append("sso_observation")
    if asteroid_orbits:
        # Asteroid orbits and covariances
        archives.append("sso_source")

    for catalog_name in archives:
        parquet_path = destination_path / f"{catalog_name}.parquet"
        if parquet_path.exists() and not force:
            print(f"{parquet_path} already exists, skipping download")
            continue

        print(f"Downloading {catalog_name}...")
        csv_files = _download_gaia_data(
            catalog_name=catalog_name,
            destination_path=destination_path,
            force=force,
        )
        _generate_parquet_from_download(csv_files=csv_files, parquet_path=parquet_path)


def _download_gaia_data(
    catalog_name: str,
    destination_path: Path,
    max_tries: int = 3,
    force: bool = False,
) -> list[Path]:
    """Download all .csv.gz files of a Gaia catalog from https://cdn.gea.esac.esa.int/ and return their paths"""
    base_url = f"https://cdn.gea.esac.esa.int/Gaia/{_DATA_RELEASE}/Solar_system/{catalog_name}/"

    # The MD5 file lists every file in the catalog, so it also serves as the list of files to download
    md5_file = destination_path / "_MD5SUM.txt"
    download_extra_file(
        url=base_url + md5_file.name,
        dest_path=md5_file,
        cache_dir=DEFAULT_CACHE,
        force=True,
    )
    md5sums = {}
    with open(md5_file) as f:
        for line in f:
            if line.strip():
                md5, md5_file_name = line.split()
                md5sums[md5_file_name] = md5
    md5_file.unlink()  # Remove MD5 file

    file_names = sorted(name for name in md5sums if name.endswith(".csv.gz"))
    if not file_names:
        raise RuntimeError(
            f"No files listed for {catalog_name} in Gaia data release {_DATA_RELEASE}"
        )

    csv_files = []
    for file_name in file_names:
        url = base_url + file_name
        path_to_file = destination_path / file_name
        csv_files.append(path_to_file)

        # Keep intact files left over from an earlier, interrupted run
        if path_to_file.exists() and not force and _verify_md5(path_to_file, md5sums[file_name]):
            continue

        for i in range(max_tries):

            download_extra_file(
                url=url, dest_path=path_to_file, cache_dir=DEFAULT_CACHE, force=True
            )

            if _verify_md5(path_to_file, md5sums[file_name]):
                break

            if i < max_tries - 1:
                print(f"Download corrupted, retrying... {i+2}/{max_tries}")
            else:
                raise RuntimeError(
                    f"Checksum of {file_name} does not match after {max_tries} attempts. Aborting download"
                )

    return csv_files


def _generate_parquet_from_download(
    csv_files: list[Path],
    parquet_path: Path,
    rows_per_chunk: int = 2_000_000,
) -> None:
    """
    Combine the downloaded .csv.gz files into one Apache .parquet file, sorted by MPC number, and clean up compressed
    files afterwards.

    To prevent loading all data into memory at once, each file is first sorted and written to
    a temporary .parquet file, after which the temporary files are merged one range of MPC numbers at a time.

    Parameters
    ----------
    csv_files : list[Path]
        Paths to the downloaded .csv.gz files
    parquet_path : Path
        Path of the generated parquet file
    rows_per_chunk : int, optional
        Approximate number of rows held in memory while merging, by default 2 000 000
    """
    with tempfile.TemporaryDirectory(dir=parquet_path.parent) as temp_dir:
        temp_dir = Path(temp_dir)

        # Convert each file to a temporary .parquet file sorted by MPC number, one file at a time
        sorted_files = []
        number_mp_counts = []
        column_types = None
        for i, csv_file in enumerate(
            tqdm(csv_files, desc=f"Sorting {parquet_path.stem}", unit="file")
        ):
            table = pv.read_csv(
                csv_file,
                read_options=pv.ReadOptions(skip_rows=_count_comment_lines(csv_file)),
                convert_options=pv.ConvertOptions(
                    column_types=column_types, strings_can_be_null=True
                ),
            )
            if column_types is None:
                # Read all other files with the column types of the first, so they share one schema
                column_types = table.schema
                array_columns = [
                    name for name in table.column_names if _is_literal_array(table[name])
                ]

            # Arrays stored as text, e.g. '[1.0,2.0]', are converted to lists
            for name in array_columns:
                table = table.set_column(
                    table.schema.get_field_index(name), name, _literal_array_to_list(table[name])
                )

            table = table.sort_by("number_mp")
            number_mp_counts.append(table["number_mp"].to_pandas().value_counts())
            sorted_files.append(temp_dir / f"sorted_{i}.parquet")
            pq.write_table(
                table, sorted_files[-1], row_group_size=32_768, **_parquet_options(table.schema)
            )
            del table

        # Merge the sorted files, one range of MPC numbers at a time
        counts = pd.concat(number_mp_counts).groupby(level=0).sum()
        dataset = ds.dataset(sorted_files)
        number_mp = ds.field("number_mp")
        merged_file = temp_dir / parquet_path.name
        chunks = counts.groupby((counts.cumsum() - 1) // rows_per_chunk)
        with pq.ParquetWriter(
            merged_file, dataset.schema, **_parquet_options(dataset.schema)
        ) as writer:
            for _, chunk_counts in tqdm(
                chunks, total=chunks.ngroups, desc=f"Merging {parquet_path.stem}", unit="chunk"
            ):
                chunk = dataset.to_table(
                    filter=(number_mp >= chunk_counts.index[0])
                    & (number_mp <= chunk_counts.index[-1])
                )
                writer.write_table(chunk.sort_by("number_mp"), row_group_size=250_000)

        # Only move the file into place once complete, so an interrupted run is not mistaken for a finished one
        merged_file.replace(parquet_path)

    # Clean up gz files
    for csv_file in csv_files:
        csv_file.unlink()

    print(f"{parquet_path.name} saved at {parquet_path.parent}")


def _count_comment_lines(csv_file: Path) -> int:
    """Count the lines starting with '#' at the top of a .csv.gz file (the ECSV header)"""
    with gzip.open(csv_file, "rt") as file:
        return sum(1 for _ in itertools.takewhile(lambda line: line.startswith("#"), file))


def _is_literal_array(column: pa.ChunkedArray) -> bool:
    """Check if a column contains arrays stored as text, i.e. all its values start with '['"""
    if not pa.types.is_string(column.type):
        return False
    values = column.drop_null()
    return len(values) > 0 and pc.all(pc.starts_with(values, "[")).as_py()


def _literal_array_to_list(column: pa.ChunkedArray) -> pa.ChunkedArray:
    """Convert arrays stored as text, e.g. '[1.0,2.0]', to lists of floats"""
    return pc.cast(pc.split_pattern(pc.utf8_trim(column, "[]"), ","), pa.list_(pa.float64()))


def _parquet_options(schema: pa.Schema) -> dict:
    """
    Settings for writing .parquet files. Storing floats with the byte-stream-split encoding and zstd compression makes
    the files roughly 40% smaller than the default settings, without loss of precision.
    """
    float_columns = [field.name for field in schema if pa.types.is_floating(field.type)]
    return dict(
        compression="zstd",
        use_byte_stream_split=float_columns,
        use_dictionary=[name for name in schema.names if name not in float_columns],
    )
