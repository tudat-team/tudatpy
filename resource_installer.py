import argparse
import sys
import requests
from pathlib import Path
from resource_data_registry import DATA_REGISTRY
import zenodo_get

ZENODO_DOI = "10.5281/zenodo.1234567" # sample DOI

BUNDLED_FILES = {
    "spice_kernels/inpop19a_TCB_m100_p100_asc": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TCB_m100_p100_asc.tar.gz",
    "spice_kernels/inpop19a_TDB_m100_p100_asc": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TDB_m100_p100_asc.tar.gz",
    "spice_kernels/inpop19a_TDB_m100_p100_spice": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TDB_m100_p100_spice.tar.gz"
}

# Downloads a file from a URL
def download_file(name: str, url: str, dest_path: Path) -> None:
    print(f"[DOWNLOAD] Fetching {url} -> {dest_path}")
    # To be implemented
    # Assumes dest_path is valid (checks & creation of missing dirs happen before)
    # Ideally it should work with URLs from: GitHub, External
    # If one-size fits all is not possible, consider different functions for each data provider

def download_file_from_zenodo(name: str, dest_path: Path) -> None:
    print(f"[DOWNLOAD] Fetching {name} from {ZENODO_DOI} -> {dest_path}")
    # To be implemented
    # Assumes dest_path is valid (checks & creation of missing dirs happen before)
    # Consider using zenodo_get library for easy download
    zenodo_get.download(
        record_or_doi=ZENODO_DOI,
        output_dir=Path._str
    )

def download_files_from_tarball(filenames: list[str], tarball_url: str, dest_path: Path) -> None:
    # Install filenames in dest_path (directory, assume valid)
    # Files are inside a tarball
    # Files in tarball that are not in the filenames list should be discarded
    # if filenames list is empty, download and install all files in tarball
    return 0


# Installs all files; overwrite any exisiting files on disc
def install_all_files(dest_path: Path, from_origin_source: bool = False) -> None:
    # Install all files that can be downloaded individually
    prefixes_to_exclude = tuple(f"{k}/" for k in BUNDLED_FILES.keys())
    filtered_registry = {
        key: value 
        for key, value in DATA_REGISTRY.items() 
        if not key.startswith(prefixes_to_exclude)
    }
    install_files(filtered_registry, dest_path, from_origin_source, True)

    # Install files that come bundled
    if from_origin_source:
        for dir, url in BUNDLED_FILES:
            download_files_from_tarball([], url, dest_path)
    #TODO: else (from Zenodo)


# Installs a given subset of files
def install_files(files: dict, dest_path: Path, from_origin_source: bool = False, force: bool = False) -> None:
    return 0
    # TODO: implement




def find_in_registry(search_string: str) -> dict:
    """
    Searches the DATA_REGISTRY for any keys containing search_string
    and returns a new dictionary with all matching keys and their values.
    """
    return {
        key: value 
        for key, value in DATA_REGISTRY.items() 
        if search_string in key
    }

def main():

    # TODO 1: Implement CLI arguments
    usr_search_str = "earth_deformation"
    from_origin_source = True
    force = True
 
    # TODO 2: retrieve installation directory from env variable
    dest_path = "~/.tudat_resources"

    # TODO 3. if usr_search_stry is empty, update all files
    #    Function to check to check for missing files locally (conda)
    
    # User did not provide a list of files to update; update all
    if not usr_search_str:
        files_to_update = DATA_REGISTRY

    files_to_update = find_in_registry(usr_search_str)

    if not files_to_update:
        raise RuntimeError(f"'{usr_search_str}' is not part of tudat-resources.")

    # 3. Create file/dir install path if it does not exist in local directory

    install_files(files_to_update, from_origin_source, force)



if __name__ == "__main__":
    main()