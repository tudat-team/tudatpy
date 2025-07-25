import json
import getpass
import os
from collections import defaultdict
import numpy as np
from tudatpy.dynamics import environment, environment_setup, propagation_setup
from datetime import datetime, timedelta
import requests


class SpaceTrackQuery:
    def __init__(self, username=None, password=None):
        if not username:
            username = input('space-track username: ')
        if not password:
            password = getpass.getpass('space-track password: ')

        self.username = username
        self.password = password
        self.spacetrack_url = 'https://www.space-track.org'
        self.download_tle = self.DownloadTle(self)  # Pass parent to subclass

    class DownloadTle:
        def __init__(self, parent):
            self.parent = parent  # Store reference to parent SpaceTrackQuery

            # Optionally prompt again if parent credentials aren't available
            if not self.parent.username:
                self.parent.username = input('space-track username: ')
            if not self.parent.password:
                self.parent.password = getpass.getpass('space-track password: ')

            self.username = self.parent.username
            self.password = self.parent.password
            self.spacetrack_url = self.parent.spacetrack_url

        def latest_on_orbit(self):

            """
            This function retrieves the newest propagable element set for all on-orbit objects (according to Space-Track.org)
            :return: None
            This is a static function that creates a .json file as output
            """
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                response = s.get(os.path.join(self.spacetrack_url,'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/decay_date/null-val/epoch/%3Enow-30/orderby/norad_cat_id/format/json'))
                if response.status_code == 200:
                    new_data = response.json()
                    json_name = "latest_on_orbit.json"
                    try:
                        with open(json_name, "r") as json_file:
                            existing_data = json.load(json_file)
                    except (FileNotFoundError, json.decoder.JSONDecodeError):
                        existing_data = []

                    existing_data.append(new_data)

                    with open(json_name, "w") as json_file:
                        json.dump(existing_data, json_file, indent=4)
                        print(f"Data appended to {json_name} file.")
                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

        def descending_epoch(self, N = None):

            """
            This function retrieves the general perturbations (GP) class
            (newest SGP4 keplerian element set for each man-made earth-orbiting object tracked by the 18th Space Defense Squadron)

            :param N: Number of first N objects (optional) to download. Default is None (full catalog is downloaded).
            :return: None

            This is a static function that creates a .json file as output
            """
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                if N is not None and float(N):
                    json_name = f"gp_descending_{N}.json"
                    response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc/limit/{N}/format/json'))
                else:
                    json_name = f"gp_descending.json"
                    response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/orderby/epoch desc/format/json'))

                if response.status_code == 200:
                    new_data = response.json()
                    try:
                        with open(json_name, "r") as json_file:
                            existing_data = json.load(json_file)
                    except (FileNotFoundError, json.decoder.JSONDecodeError):
                        existing_data = []

                    existing_data.append(new_data)

                    with open(json_name, "w") as json_file:
                        json.dump(existing_data, json_file, indent=4)
                        print(f"Data appended to {json_name} file.")
                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

        def single_norad_id(self, norad_id, N = 1):

            """
            This function retrieves the general perturbations (GP or GP history) class
            (newest SGP4 keplerian element set for each man-made earth-orbiting object tracked by the 18th Space Defense Squadron)

            :param N: Number of N TLEs (optional) to download. Default is 1.
            :return: None

            This is a static function that creates a .json file as output
            """
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                if N is not None and float(N):
                    json_name = f"single_{norad_id}_{N}.json"
                    if N == 1:
                        response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/NORAD_CAT_ID/{norad_id}/orderby/epoch desc/limit/{N}/format/json'))
                    else:
                        response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp_history/OBJECT_TYPE/PAYLOAD/NORAD_CAT_ID/{norad_id}/orderby/epoch desc/limit/{N}/format/json'))

                if response.status_code == 200:
                    new_data = response.json()
                    try:
                        with open(json_name, "r") as json_file:
                            existing_data = json.load(json_file)
                    except (FileNotFoundError, json.decoder.JSONDecodeError):
                        existing_data = []

                    existing_data.append(new_data)

                    with open(json_name, "w") as json_file:
                        json.dump(existing_data, json_file, indent=4)

                    # SOMETIMES, Duplicate TLE for a given EPOCH for a given OBJECT are found in the gp_history.
                    # The following function filters the .json created above, keeping only the latest TLEs
                    # among the duplicated ones (latest creation date)

                    updated_N = self.parent.TleUtils.filter_tles_keep_latest_creation_from_json(self.parent, json_name)
                    parts_json_name = json_name.split('.')[0].split('_')
                    parts_json_name[-1] = f'{updated_N}'
                    new_name = '_'.join(parts_json_name) + '.json'
                    os.system(f'mv temp_filtered_tles.json {new_name}')

                    if new_name != json_name:
                        os.system(f'rm {json_name}')
                    print(f"Data appended to {new_name} file. (Name was: {json_name}.)")

                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

        def filtered_by_oe_dict(self, filter_oe_dict, limit=100, output_file='filtered_results.json'):

            base_url = "basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/"
            session = requests.Session()
            session.post(self.spacetrack_url + "/ajaxauth/login", json={'identity': self.username, 'password': self.password})

            all_results = None

            for oe, bounds in filter_oe_dict.items():
                min_val = bounds[0] if len(bounds) > 0 else None
                max_val = bounds[1] if len(bounds) > 1 else None

                results_min, results_max = [], []

                if min_val is not None:
                    url_min = f"{base_url}{oe}/>{min_val}/orderby/epoch desc/limit/{limit}/format/json/"
                    resp_min = session.get(os.path.join(self.spacetrack_url, url_min))
                    if resp_min.status_code == 200:
                        results_min = resp_min.json()

                if max_val is not None:
                    url_max = f"{base_url}{oe}/<{max_val}/orderby/epoch desc/limit/{limit}/format/json"
                    resp_max = session.get(os.path.join(self.spacetrack_url, url_max))
                    if resp_max.status_code == 200:
                        results_max = resp_max.json()

                if min_val is not None and max_val is not None:
                    ids_min = set(obj['NORAD_CAT_ID'] for obj in results_min)
                    ids_max = set(obj['NORAD_CAT_ID'] for obj in results_max)
                    common_ids = ids_min.intersection(ids_max)
                    filtered = [obj for obj in results_min if obj['NORAD_CAT_ID'] in common_ids]
                elif min_val is not None:
                    filtered = results_min
                else:
                    filtered = results_max

                if all_results is None:
                    all_results = filtered
                else:
                    current_ids = set(obj['NORAD_CAT_ID'] for obj in filtered)
                    all_results = [obj for obj in all_results if obj['NORAD_CAT_ID'] in current_ids]

            if not all_results:
                print('No object found satisfying all user-defined filters. Exiting...\n')
                exit()

            with open(output_file, 'w') as f_out:
                json.dump([all_results], f_out, indent=4)


            return all_results



    class TleUtils:
        def __init__(self, parent):
            self.parent = parent  # Store reference to parent SpaceTrackQuery

            # Optionally prompt again if parent credentials aren't available
            if not self.parent.username:
                self.parent.username = input('space-track username: ')
            if not self.parent.password:
                self.parent.password = getpass.getpass('space-track password: ')

            self.username = self.parent.username
            self.password = self.parent.password
            self.spacetrack_url = self.parent.spacetrack_url

        def get_norad_id_name_map(self, limit=None):
            """
            This function gives an norda_cat_id - object_name mapping of the full satcat catalogue.
            Users only need to run it once.

            :param limit [optional], limits the numbers of considered satellites.
            :return: Static, creates a .json file containing the mapping dictionary.
            """
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                # Download satcat data
                query = f"{self.spacetrack_url}/basicspacedata/query/class/satcat/orderby/NORAD_CAT_ID/format/json"
                if limit:
                    query = f"{self.spacetrack_url}/basicspacedata/query/class/satcat/orderby/NORAD_CAT_ID/limit/{limit}/format/json"

                response = s.get(query)
                if response.status_code != 200:
                    raise RuntimeError(f"Failed to fetch data: {response.status_code}")

                data = response.json()
                mapping = {int(obj['NORAD_CAT_ID']): obj['OBJECT_NAME'] for obj in data if obj['OBJECT_NAME']}
                # Create the dictionary mapping
                json_name = 'norad_id_to_name.json'
                with open(json_name, "w") as f:
                    json.dump(mapping, f, indent=4)
                    print(f"Saved {len(mapping)} entries to {json_name}.")

        def filter_tles_keep_latest_creation_from_json(self, json_filename):
            """
            Reads TLE data from a JSON file, filters duplicates by EPOCH keeping
            only the latest CREATION_DATE for each unique EPOCH,
            and returns the filtered list of TLE dictionaries.

            Used within the single_norad_id function.

            :param json_filename: path to JSON file
            :return: Static, it modifies the input .json file
            """
            with open(json_filename, "r") as f:
                objects_list = json.load(f)

            filtered = {}
            for i,object in enumerate(objects_list[0]):
                epoch = object['EPOCH']
                creation_date = datetime.fromisoformat(object['CREATION_DATE'])

                if epoch not in filtered:
                    filtered[epoch] = object
                else:
                    print(f'Number of TLEs in the JSON file for OBJECT {object["NORAD_CAT_ID"]}: {len(objects_list[0])}')
                    print(f'Warning: Found Duplicate TLE for EPOCH {epoch} for OBJECT {object["NORAD_CAT_ID"]}. Keeping only the latest one...')
                    print(f'Updating JSON file accordingly...')
                    existing_creation_date = datetime.fromisoformat(filtered[epoch]['CREATION_DATE'])
                    if creation_date > existing_creation_date:
                        filtered[epoch] = object

            with open("temp_filtered_tles.json", "w") as f_out:
                json.dump([list(filtered.values())], f_out, indent=4)

            print(f'Number of Filtered TLEs in the JSON file for OBJECT {object["NORAD_CAT_ID"]}: {len(filtered.values())}')

            return len(filtered.values())

        def get_tle_dict_from_json(self, json_filename):

            """
            THIS HAS TO BE FIXED BECAUSE IT ONLY WORKS WELL WITH THE GP_DESCENDING_5.json file, and not with the SINGLE_5000_4.json
            """

            tle_dict = defaultdict(list)
            with open(json_filename, "r") as f:
                objects_list = json.load(f)

            for i,object in enumerate(objects_list[0]):
                norad_cat_id = object['NORAD_CAT_ID']
                tle_line_1 = object['TLE_LINE1']
                tle_line_2 = object['TLE_LINE2']

                tle_dict[norad_cat_id].append((tle_line_1, tle_line_2))

            return tle_dict

        def plot_earth(self, ax, radius=6378, color='lightblue', alpha=0.5, resolution=50):
            """Plot Earth as a sphere in the 3D axes."""
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            x = radius * np.outer(np.cos(u), np.sin(v))
            y = radius * np.outer(np.sin(u), np.sin(v))
            z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, rstride=1, cstride=1, color=color, alpha=alpha, edgecolor='none')

        def tle_to_TleEphemeris_object(self, tle_line_1, tle_line_2):
            object_tle = environment.Tle(tle_line_1, tle_line_2)
            ephemeris_object = environment.TleEphemeris("Earth", "J2000", object_tle, False)
            return ephemeris_object

        def get_tle_reference_epoch(self, tle_line_1):

            """
            Extracts and converts the epoch from TLE line 1 to a UTC datetime object.

            Parameters:
                tle_line1 (str): The first line of a TLE string.

            Returns:
                datetime: The epoch in UTC.
            """

            year_str = tle_line_1[18:20]
            day_str = tle_line_1[20:32]

            # Convert to full year
            year = int(year_str)
            year += 2000 if year < 57 else 1900

            # Convert fractional day to datetime
            day_of_year = float(day_str)
            epoch = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)

            return epoch
