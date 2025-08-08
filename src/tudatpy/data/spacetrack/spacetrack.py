import json
import getpass
import os
from collections import defaultdict
import numpy as np
from tudatpy.dynamics import environment, environment_setup, propagation_setup
from datetime import datetime, timedelta
import requests
import math
from tudatpy.astro import time_representation
import warnings

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
        self.omm_utils = self.OMMUtils(self)

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
            :return: dictionary corresponding to the downloaded json file
            """
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                response = s.get(os.path.join(self.spacetrack_url,'basicspacedata/query/class/gp/OBJECT_TYPE/PAYLOAD/decay_date/null-val/epoch/%3Enow-30/orderby/norad_cat_id/format/json'))
                json_name = "latest_on_orbit.json"

                if response.status_code == 200:
                    queried_data = response.json()
                    if os.path.exists(json_name):
                        os.remove(json_name)
                    json_data = queried_data
                    with open(json_name, "w") as json_file:
                        json.dump(json_data, json_file, indent=4)
                        print(f"Data downloaded to {json_name} file.")

                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

            return json_data

        def descending_epoch(self, N = None):

            """
            This function retrieves the general perturbations (GP) class
            (newest SGP4 keplerian element set for each man-made earth-orbiting object tracked by the 18th Space Defense Squadron)

            :param N: Number of first N objects (optional) to download. Default is None (full catalog is downloaded).
            :return: dictionary corresponding to the downloaded json data
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
                    queried_data = response.json()
                    if os.path.exists(json_name):
                        os.remove(json_name)
                    json_data = queried_data
                    with open(json_name, "w") as json_file:
                        json.dump(json_data, json_file, indent=4)

                    if N != 1:
                        filtered_values = self.parent.OMMUtils.filter_tles_keep_latest_creation_from_json(self, json_name)
                        parts_json_name = json_name.split('.')[0].split('_')
                        parts_json_name[-1] = f'{len(filtered_values)}'
                        new_name = '_'.join(parts_json_name) + '.json'
                        os.system(f'mv temp_filtered_tles.json {new_name}')

                        if new_name != json_name:
                            os.system(f'rm {json_name}')
                        print(f"Data downloaded to {new_name} file.")

                        print(f"Data downloaded to {json_name} file.")

                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

            return json_data

        def single_norad_id(self, norad_id, N = 1):

            """
            This function retrieves the general perturbations (GP or GP history) class
            (newest SGP4 keplerian element set for each man-made earth-orbiting object tracked by the 18th Space Defense Squadron)

            :param N: Number of N TLEs (optional) to download. Default is 1.
            :return: dictionary corresponding to the downloaded json data
            """
            filtered_dict = None
            with requests.Session() as s:
                s.post(self.spacetrack_url+"/ajaxauth/login", json={'identity':self.username, 'password':self.password})
                if N is not None and float(N):
                    json_name = f"single_{norad_id}_{N}.json"
                    if N == 1:
                        response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp/NORAD_CAT_ID/{norad_id}/orderby/epoch desc/limit/{N}/format/json'))
                    else:
                        response = s.get(os.path.join(self.spacetrack_url,f'basicspacedata/query/class/gp_history/NORAD_CAT_ID/{norad_id}/orderby/epoch desc/limit/{N}/format/json'))

                if response.status_code == 200:
                    queried_data = response.json()
                    if os.path.exists(json_name):
                        os.remove(json_name)
                    json_data = queried_data
                    with open(json_name, "w") as json_file:
                        json.dump(json_data, json_file, indent=4)

                    # SOMETIMES, Duplicate TLE for a given EPOCH for a given OBJECT are found in the gp_history.
                    # The following function filters the .json created above, keeping only the latest TLEs
                    # among the duplicated ones (latest creation date)

                    if N != 1:
                        filtered_values = self.parent.OMMUtils.filter_tles_keep_latest_creation_from_json(self, json_name)
                        parts_json_name = json_name.split('.')[0].split('_')
                        parts_json_name[-1] = f'{len(filtered_values)}'
                        new_name = '_'.join(parts_json_name) + '.json'
                        os.system(f'mv temp_filtered_tles.json {new_name}')

                        if new_name != json_name:
                            os.system(f'rm {json_name}')
                        print(f"Data downloaded to {new_name} file.")

                else:
                    print("Failed to retrieve data from the API. Status code:", response.status_code)

            if filtered_dict:
                return filtered_dict
            else:
                return json_data

        def filtered_by_oe_dict(self, filter_oe_dict, limit=100, output_file='filtered_results.json'):

            """
            :param filter_oe_dict: dictionary containing the orbital elements as keys and the list of min and max bounds as a value
            :param limit: optional, retrieves the first N = limit objects
            :param output_file: optional, name of output file
            :return: dictionary corresponding to the downloaded filtered data
            """
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



    class OMMUtils:
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
            self.supported_orbital_regimes = ['LEO_REGIME', 'MEO_REGIME', 'GEO_REGIME', 'GSO_REGIME', 'HEO_REGIME', 'OTHER']

        def get_norad_id_name_map(self, limit=None):
            """
            This function gives a norad_cat_id - object_name mapping of the full satcat catalogue.
            Users only need to run it once.

            :param limit [optional], limits the numbers of considered satellites.
            :return: Returns norad_id to name mapping.
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
                map = {int(obj['NORAD_CAT_ID']): obj['OBJECT_NAME'] for obj in data if obj['OBJECT_NAME']}
                # Create the dictionary mapping
                json_name = 'norad_id_to_name.json'
                with open(json_name, "w") as f:
                    json.dump(map, f, indent=4)
                    print(f"Saved {len(map)} entries to {json_name}.")

            return map

        def filter_tles_keep_latest_creation_from_json(self, json_filename):
            """
            Reads TLE data from a JSON file, filters duplicates by EPOCH keeping
            only the latest CREATION_DATE for each unique EPOCH,
            and returns the filtered list of TLE dictionaries.

            Used within the single_norad_id function.

            :param json_filename: path to JSON file
            :return: Dictionary corresponding to the newly created json
            """
            with open(json_filename, "r") as f:
                objects_list = json.load(f)

            filtered = {}
            for i, object_ in enumerate(objects_list):
                epoch = object_['EPOCH']

                creation_date = datetime.fromisoformat(object_['CREATION_DATE'])

                if epoch not in filtered:
                    filtered[epoch] = object_
                else:
                    print(f'Number of TLEs in the JSON file for OBJECT {object_["NORAD_CAT_ID"]}: {len(objects_list[0])}')
                    print(f'Found Duplicate TLE for EPOCH {epoch} for OBJECT {object_["NORAD_CAT_ID"]}. Keeping only the latest one...')
                    print(f'Updating JSON file accordingly...')
                    existing_creation_date = datetime.fromisoformat(filtered[epoch]['CREATION_DATE'])
                    if creation_date > existing_creation_date:
                        filtered[epoch] = object_

            with open("temp_filtered_tles.json", "w") as f_out:
                json.dump([list(filtered.values())], f_out, indent=4)

            print(f'Number of Filtered TLEs in the JSON file for OBJECT {object_["NORAD_CAT_ID"]}: {len(filtered.values())}')

            return filtered.values()

        def get_tles(self, json_dict):
            tle_dict = defaultdict(list)

            if type(json_dict) is list:
                for json in json_dict:
                    norad_cat_id = json['NORAD_CAT_ID']
                    tle_line_1 = json['TLE_LINE1']
                    tle_line_2 = json['TLE_LINE2']
                    tle_dict[norad_cat_id] = (tle_line_1, tle_line_2)

            return tle_dict

        def get_tudat_keplerian_element_set(self, json_dict):

            # retrieve orbital elements from json_dict
            a = float(json_dict[0]['SEMIMAJOR_AXIS'])*1e3 # [meters]
            e =float(json_dict[0]['ECCENTRICITY'])
            i = float(json_dict[0]['INCLINATION']) * math.pi/180 # [rad]
            omega = float(json_dict[0]['ARG_OF_PERICENTER']) * math.pi/180 # [rad]
            raan = float(json_dict[0]['RA_OF_ASC_NODE']) * math.pi/180 # [rad]
            mo = float(json_dict[0]['MEAN_ANOMALY']) * math.pi/180 # [rad]

            # Compute true anomaly from mean anomaly via Kepler's equation
            true_anomaly = self.mean_to_true_anomaly(self,mo,e)
            
            return a,e,i,omega,raan,true_anomaly

        def tle_to_sgp4_ephemeris_object(self, tle_line_1, tle_line_2, simulation_start_epoch = None, simulation_end_epoch = None, timestep_global = 5, frame_origin = 'Earth', frame_orientation = 'J2000'):
            # Checking start and end of simulation
            reference_epoch_tle = time_representation.datetime_to_tudat(self.get_tle_reference_epoch(tle_line_1)).epoch()
            if not simulation_start_epoch and not simulation_end_epoch:
                simulation_start_epoch = reference_epoch_tle
                simulation_start_epoch_datetime = time_representation.DateTime.to_python_datetime(time_representation.DateTime.from_epoch(simulation_start_epoch))
                simulation_end_epoch = time_representation.datetime_to_tudat(simulation_start_epoch_datetime + timedelta(hours=5)).epoch()
                warnings.warn(
                    "No simulation start nor end epoch provided.\n"
                    "Starting simulation at TLE reference epoch.\n"
                    "Ending simulation at TLE reference epoch + 5 hours.",
                    UserWarning
                )
            sgp4_ephemeris =  environment_setup.ephemeris.sgp4(
                tle_line_1,
                tle_line_2,
                frame_origin = 'Earth',
                frame_orientation = 'J2000')

            return sgp4_ephemeris

        def plot_earth(self, ax, radius=6378, color='lightblue', alpha=0.5, resolution=50):
            """Plot Earth as a sphere in the 3D axes."""
            u = np.linspace(0, 2 * np.pi, resolution)
            v = np.linspace(0, np.pi, resolution)
            x = radius * np.outer(np.cos(u), np.sin(v))
            y = radius * np.outer(np.sin(u), np.sin(v))
            z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z, rstride=1, cstride=1, color=color, alpha=alpha, edgecolor='none')

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

            # add n = day_of_year days to 1st of january of that year, then subtract one.
            epoch = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)

            return epoch

        def mean_to_true_anomaly(self, mo, e, tol=1e-8, max_iter=100):
            """
            Convert mean anomaly M to true anomaly ν for elliptical orbits (solving Kepler's equation).

            Parameters:
                M : float
                    Mean anomaly (radians)
                e : float
                    Eccentricity (0 <= e < 1)
                tol : float
                    Tolerance for Newton-Raphson iteration
                max_iter : int
                    Maximum number of iterations

            Returns:
                ν : float
                    True anomaly (radians)
            """
            # Solve Kepler's Equation: M = E - e*sin(E)
            E = mo if e < 0.8 else np.pi  # Initial guess
            for _ in range(max_iter):
                f = E - e * np.sin(E) - mo
                f_prime = 1 - e * np.cos(E)
                E_new = E - f / f_prime
                if abs(E_new - E) < tol:
                    break
                E = E_new
            else:
                raise RuntimeError("Kepler's equation did not converge.")

            # Convert eccentric anomaly to true anomaly
            ν = 2 * np.arctan2(
                np.sqrt(1 + e) * np.sin(E / 2),
                np.sqrt(1 - e) * np.cos(E / 2)
            )
            return ν

        def get_orbital_regime(self, single_object_json_dict):
            """
            Classify the orbital regime of a space object based on its orbital elements.

            This function evaluates the periapsis, apoapsis, eccentricity, and inclination of a space object
            (assumed to orbit the Earth) and returns the corresponding orbital regime along with the regime's
            threshold definitions.

            The supported regimes are:
                - LEO_REGIME: Low Earth Orbit (100–2000 km)
                - MEO_REGIME: Medium Earth Orbit (2000–35786 km)
                - GEO_REGIME: Geostationary Orbit (~35786 km with low eccentricity and inclination)
                - GSO_REGIME: Geosynchronous Orbit (~35786 km but not meeting GEO constraints)
                - HEO_REGIME: Highly Elliptical Orbit (e.g., Molniya-type)
                - "Other": Any orbit that does not fall into the above categories

            Parameters
            ----------
            single_object_json_dict : dict
                A dictionary containing orbital elements for a single object.
                Required keys:
                    - 'CENTER_NAME': must be 'EARTH'
                    - 'PERIAPSIS' (km above surface)
                    - 'APOAPSIS' (km above surface)
                    - 'ECCENTRICITY' (unitless)
                    - 'INCLINATION' (degrees)

            Returns
            -------
            orbital_regime : str
                The name of the identified orbital regime (e.g., "LEO_REGIME", "MEO_REGIME", etc.)

            regime_bounds : dict
                The threshold values used to classify the orbit, including altitude and (if applicable)
                eccentricity and inclination ranges.

            Raises
            ------
            ValueError
                If the object is not orbiting Earth or required keys are missing or invalid.
            """

            REGIME_THRESHOLDS = {
                "LEO_REGIME": {
                    "rp": (100, 2000),
                    "ra": (100, 2000),
                },
                "MEO_REGIME": {
                    "rp": (2000, 35786),
                    "ra": (2000, 35786),
                },
                "GEO_REGIME": {
                    "rp": (35586, 35986),  # 35786 ± 200 km
                    "ra": (35586, 35986),
                    "ecc": (0.0, 0.01),
                    "inc": (0.0, 1.0),     # degrees
                },
                "GSO_REGIME": {
                    "rp": (35586, 35986),  # Same altitude range as GEO
                    "ra": (35586, 35986),
                    # No constraints on ecc/inc — anything not satisfying GEO will fall here
                },
                "HEO_REGIME": {
                    "rp": (100, 10000),
                    "ra": (35000, 50000),
                }
            }

            # Extract orbital elements from JSON
            central_body_name = single_object_json_dict["CENTER_NAME"]
            if central_body_name != 'EARTH':
                raise ValueError(f'Central body name must be "EARTH" (Central body: {central_body_name} not supported).')

            rp_object = float(single_object_json_dict.get('PERIAPSIS', None))  # km above surface
            ra_object = float(single_object_json_dict.get('APOAPSIS', None))   # km above surface
            ecc = float(single_object_json_dict.get('ECCENTRICITY', None))
            inc = float(single_object_json_dict.get('INCLINATION', None))

            # Initialize result
            orbital_regime = "OTHER"

            # Loop through defined regimes
            for regime, thresholds in REGIME_THRESHOLDS.items():
                rp_min, rp_max = thresholds["rp"]
                ra_min, ra_max = thresholds["ra"]

                if not (rp_min <= rp_object <= rp_max and ra_min <= ra_object <= ra_max):
                    continue

                # Check eccentricity and inclination if required
                ecc_range = thresholds.get("ecc", None)
                inc_range = thresholds.get("inc", None)

                if ecc_range and (ecc is None or not ecc_range[0] <= ecc <= ecc_range[1]):
                    continue

                if inc_range and (inc is None or not inc_range[0] <= inc <= inc_range[1]):
                    continue

                # All checks passed
                orbital_regime = regime
                break

            return orbital_regime, REGIME_THRESHOLDS[orbital_regime]