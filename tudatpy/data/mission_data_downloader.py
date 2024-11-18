#!/usr/bin/env python
# coding: utf-8

# ## Class for Loading PDS files

import sys
import numpy as np
import requests
from bs4 import BeautifulSoup
import os
from urllib.request import urlretrieve
from datetime import datetime, timedelta
import glob
import re
from tudatpy.interface import spice
from dateutil.relativedelta import relativedelta

class LoadPDS:

    def __init__(self):  
        self.check = False
        self.supported_patterns = {
            "juice":
                {"ck": r"^(?P<mission>juice)_(?P<instrument>(sc|mga|sa))_(?P<data_type>(meas|crema_\w+_baseline))?_?(?P<reference_file>[a-zA-Z]+)?_?(?P<prod_ID>[a-zA-Z]+)?_?(?P<start_date_file>(\d{6}|\d{8}))?_?(?P<end_date_file>\d{6})?_?(?P<sclk>(s|t|f)\d{6})?_?(?P<version>v\d{2})(?P<extension>\.bc)$",
                "spk": r"^(?P<mission>juice)_(?P<data_type>(orbc|crema))_(?P<reference_file_prod_ID>(\d{6}|\w+_plan))?_?(?P<start_date_file>\d{6})?_?(?P<end_date_file>\d{6})?_?(?P<version>v\d{2})(?P<extension>\.bsp)$",
                "fk": r"^(?P<mission>juice)_(?P<data_type>(dsk_surfaces|events_\w+|ops|roi|sci|stations_topo))_(?P<version>v\d{2})(?P<extension>\.tf)$"
                },
            
            "mro":
                {"ck": r"^(?P<mission>mro)_(?P<instrument>(sc|hga|sa))_(?P<mission_phase>(cru|ab|psp|rel|psp))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<version>v\d{2})?(?P<extension>\.bc)$",
                "spk": r"^(?P<mission>mro)_(?P<instrument>(sc|hga|sa|struct))_(?P<mission_phase>(cru|ab|psp|rel|psp))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<version>v\d{2})?(?P<extension>\.bsp)$",
                "odf":
                r'^(?P<mission>mro)(?P<dataset>magr)(?P<date_file>\d{4}_\d{3}_\d{4})(?P<uplink>x|n)(?P<station>nn|mm|[0-9]{2})(?P<downlink>[123m])(?P<version>v\d+)(?P<extension>\.odf)$',
                "tro": r"^(?P<mission>mro)(?P<dataset>magr)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.tro)$",
                "ion": r"^(?P<mission>mro)(?P<dataset>magr)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.ion)$"
                },
            "lro":
                {"ck": r"^(?P<mission>lro)(?P<instrument>(sc|hg|sa|dv|))_(?P<start_date_file>\d{7})_(?P<end_date_file>\d{7})_(?P<version>v\d{2})(?P<extension>\.bc)$"},
            
            "cassini":
                {"ck": r"^(?P<start_date_file>\d{5})_(?P<end_date_file>\d{5})(?P<type>(p|r))(?P<version>[a-z])(?P<freeform>\w+)?(?P<extension>\.bc)$"},
            
            "insight":
                {"ck": r"^(?P<mission>insight)_(?P<purpose>(cru_rec|edl_rec|surf_ops|))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})_(?P<version>v\d{2})(?P<extension>\.bc)"
                },
            
            "mex":
             {"ck": r"^(?P<mission>MEX)?_?(?P<data>ATNM)?_?(?P<purpose>(MEASURED|T6|SA))?_?(?P<start_date_file>(\d{4}|\d{6}|P\d{12}))?_?(?P<end_date_file>\d{6})?_?(?P<sclk>S\d{6})?_?(?P<version>(V\d+|\d+))?(?P<extension>\.BC)?$",
             "spk": "^(?P<data>(ORMM|ORMF))_(?P<SPK_type>T19)_(?P<start_date_file>\d{6})_?(?P<end_date_file>\d{6})_(?P<version>\d{5})(?P<extension>\.BSP)?$", 
              "ifms": r'^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)',
             "dp2": r'^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)',
             "dpx": r'^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)',
             "dps": r'^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)'}
        }

#########################################################################################################
                        
    def generic_strptime(self, date):

        """ 
        Given a list of supported time formats, try to convert it with strptime.
        This function is used throughout the code. 
        
        For instance: 
          juice and mro orientation files support %y%m%d,
          lro supports %Y%j,
          other mro kernels might support %Y_%j,
          etc...
        """
        self.supported_time_formats = {
            "YYMMDD": (
                r"%y%m%d"
            ),
            "YYYYjjj": (
                r"%Y%j"
            ),   
            "YYYY_jjj": (
                r"%Y_%j"
            ),   
            "YYjjjhhmm": (
                r"%y%j%H%M"
            ),
            "YYYY_JJJ_hhmm":(
            r"%Y_%j_%H%M"
            ),
            "YYYY_MM_DD":(
            r"%Y_%m_%d"
            ),
            "YYYY-MM-DD":(
            r"%Y-%m-%d"
            ),
            "YYYYMMDD_":(
            r"%Y%m%d_"
            ) #added this for mex_sa_date_trick (formats %Y_%j_%H%M and %Y%m%d might be confused, so I add an artificial underscore.)
            }

        for key, str_format in self.supported_time_formats.items():
            try: 
                self.date = datetime.strptime(date,str_format)
                #print(f'converted {self.date}, with format {str_format}')
                if self.date is not None:
                    return (self.date)             
            except:
                #print('didnt work strptime')
                continue      

    def generic_strftime(self, date, format_key):

        """
        Convert a datetime object to a string based on a supported format.
        :param date_obj: datetime object to format.
        :param format_key: Key from supported_time_formats for desired output format.
        :return: Formatted date string.
        """

        self.supported_time_formats = {
            "YYMMDD": (
                r"%y%m%d"
            ),
            "YYYYjjj": (
                r"%Y%j"
            ),   
            "YYYY_jjj": (
                r"%Y_%j"
            ),   
            "YYjjjhhmm": (
                r"%y%j%H%M"
            ),
            "YYYY_JJJ_hhmm":(
            r"%Y_%j_%H%M"
            ),
            "YYYY_MM_DD":(
            r"%Y_%m_%d"
            ),
            "YYYY-MM-DD":(
            r"%Y-%m-%d"
            ),
            "YYYYMMDD_":(
            r"%Y%m%d_"
            ) #added this for mex_sa_date_trick (formats %Y_%j_%H%M and %Y%m%d might be confused, so I add an artificial underscore.)
            }
        
        if format_key in self.supported_time_formats:
            # Get the format string and use strftime to format the datetime
            str_format = self.supported_time_formats[format_key]
            self.date = date.strftime(str_format)
            return self.date
        else:
            print("Format key not supported")   
#########################################################################################################

    def get_kernels(self, input_mission, url, wanted_files, local_folder):

        """
        Download specific kernel files from a given URL if they do not already exist locally.
    
        Parameters:
        - input_mission (str): Name of the mission (for display purposes).
        - url (str): Base URL where the files are hosted.
        - wanted_files (list): List of filenames to download.
        - local_folder (str): Path to the local directory where files should be stored.
        - start_date (str, optional): Start date for filtering files (not used here but placeholder for future use).
        - end_date (str, optional): End date for filtering files (not used here but placeholder for future use).
        """
        # Print information about the download process
        data_type = url.split('/')[-2]
        ext = self.get_extension_for_data_type(data_type)
        
        # Ensure local folder exists
        os.makedirs(local_folder+data_type+'/', exist_ok=True)
        
        # Create full paths for each file to load locally
        self.files_to_load = [local_folder + data_type+ '/' + wanted_file for wanted_file in wanted_files]
        
        # Download each file if not already present
        for wanted_file, local_file in zip(wanted_files, self.files_to_load):
            if not os.path.exists(local_file):
                print(f'Downloading File: {local_file}')
                urlretrieve(url + wanted_file, local_file)
            else:
                print(f'File: {local_file} already exists in {local_folder} and will not be downloaded.')
                
        return self.files_to_load

#########################################################################################################
    
    def add_custom_mission_pattern(self, input_mission, custom_pattern):
        
        """
        THE IDEA OF THIS FUNCTION IS TO ALLOW FOR USERS TO DEFINE THEIR OWN CUSTOMISED PATTERN.
        ONCE THE PATTERN IS DEFINED, this can be added to the supported ones via:
        
        object = LoadPDS()
        object.add_custom_mission_pattern('custom_mission', 'regex_pattern')
        
        """
        custom_dict = {}
        custom_dict[input_mission] = custom_pattern
        self.supported_patterns.update(custom_dict)
        
        return self.supported_patterns

#########################################################################################################

    def is_date_in_intervals(self, date, intervals):

        """
        
        This functions checks whether a given (wanted) date 
        falls within the time interval specified in the filename. 
        It is used to download the ck spice kernels, containing an interval of dates in their name.  
        
        """
        self.in_intervals = False
        for interval in intervals:
            if ( ( date.date() >= interval[0].date() ) and ( date.date() <= interval[1].date() ) ):
                self.in_intervals = True
        return (self.in_intervals)

#########################################################################################################
    
    def match_type_extension(self, data_type, filename):
        
        extension = (filename.split(".")[-1]).lower()
        supported_radio_science_types = ['odf', 'ifms', 'dp2', 'dpx', 'dps']
        supported_spice_kernels_types = ['ck', 'dsk', 'spk', 'fk', 'mk', 'ik', 'lsk', 'pck', 'sclk']
        supported_ancillary_types = ['ion', 'tro']
        all_supported_types = supported_radio_science_types + supported_spice_kernels_types + supported_ancillary_types
    
        all_supported_tuples = [
            ('ck', 'bc'),
            ('dsk', 'bds'),
            ('spk', 'bsp'),
            ('fk', 'tpc'),
            ('mk', 'tm'),
            ('ik', 'ti'),
            ('lsk', 'tls'),
            ('pck', 'bpc'),
            ('sclk', 'tsc'),
            ('odf', 'odf'),
            ('ifms', 'tab'),
            ('dp2', 'tab'),
            ('dps', 'tab'),
            ('dpx', 'tab'),
            ('ion', 'ion'),
            ('tro', 'tro')
        ]
        
        # Check if the data_type is supported
        for type_, ext in all_supported_tuples:
            if data_type == type_:
                return ext == extension  # Return True if the extension matches, False if it does not.
                
#########################################################################################################

    def get_extension_for_data_type(self, data_type):
        # Mapping each data type to its corresponding extension
        type_to_extension = {
            'ck': 'bc',
            'dsk': 'bds',
            'spk': 'bsp',
            'fk': 'tpc',
            'mk': 'tm',
            'ik': 'ti',
            'lsk': 'tls',
            'pck': 'bpc',
            'sclk': 'tsc',
            'odf': 'odf',
            'ifms': 'tab',
            'dp2': 'tab',
            'dps': 'tab',
            'dpx': 'tab',
            'ion': 'ion',
            'tro': 'tro'
        }
        
        # Return the extension for the given data type if it's supported, else None
        return type_to_extension.get(data_type.lower())

#########################################################################################################
    
    def dynamic_download_url_files_single_time(self, input_mission, local_path, start_date, end_date, url):
        """
        This function downloads files for a single time interval based on specified mission parameters and date range.
        Parameters:
        - input_mission: Name of the mission (string).
        - local_path: Directory path where files are stored (string).
        - filename_format: The format of the filename to match (string).
        - start_date: The start date of the time interval (datetime).
        - end_date: The end date of the time interval (datetime).
        - url: The base URL for downloading files (string).
        - time_format: Format for the specified date (string).
        Returns:
        - self.relevant_files: List of relevant files that were found or downloaded.
        """

        self.relevant_files = []
        print('Checking URL:', url)
        data_type = url.split('/')[-2]      
        local_subfolder = local_path + data_type + '/'

        try: 
            supported_pattern = self.supported_patterns[input_mission][data_type]
            #print(f'supported_pattern: {supported_pattern}')
            
        except:
            raise ValueError('Pattern not found among supported patterns.')

        existing_files = self.check_existing_files(input_mission, data_type, local_subfolder, supported_pattern, start_date, end_date)

        if existing_files and self.check == False:
            print(f'--------------------------------------- EXISTING FILES CHECK ----------------------------------------------')
            print(f'The following files already exist in the folder:\n\n {existing_files}\n\n and will not be downloaded.')
            self.relevant_files.extend(existing_files)

        self.check = True
        all_dates = [start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)]
        # Initialize a dictionary to hold files from the HTML response
        files_url_dict = {}    

        try:
            reqs = requests.get(url)
            reqs.raise_for_status()
        except:
            raise ValueError(f"Error fetching data from {url}: {e}")

        # Parse links from the HTML response
        for link in BeautifulSoup(reqs.text, 'html.parser').find_all('a'):
            if not '[To Parent Directory]' in link:
                full_link = link.get('href')  # Extract the href attribute
                if isinstance(full_link, str):
                    # Check if the link is a full URL or a relative path
                    if ('/') in full_link:
                        # It's a full URL
                        filename = os.path.basename(full_link)  # Get the filename from the URL
                    else:
                        # It's a relative path
                        filename = full_link
                    
                    # Check if the filename matches the specified format
                    if self.match_type_extension(data_type, filename):  
                        try:
                            RS_dict, RS_underscores = self.parse_filename(input_mission, data_type, filename)
                            filename_to_download = self.reconstruct_filename(RS_dict, RS_underscores)
    
                            # Determine the date string format
                            date_key = RS_dict["date_file"][:5] if input_mission.lower() == 'mex' else RS_dict["date_file"][:8]
                            files_url_dict[date_key] = filename_to_download
                            
                        except Exception as e:
                            print(f"Could not parse file: {filename} - Error: {e}")
                            continue
                    else:
                        continue
            else:
                continue
                
        # Download missing files
        for date in all_dates:
            if input_mission.lower() == 'mex' :
                date_string = f"{date.year % 100:02d}{date.timetuple().tm_yday:03d}"
            elif input_mission.lower() == 'mro':
                format_key = "YYYY_jjj"
                date_string = LoadPDS.generic_strftime(LoadPDS,date, format_key)

            else:
                print('Error: I dont know what mission you are talking about!')
                
            if date_string in files_url_dict:
                download_file = files_url_dict[date_string]
                full_download_url = os.path.join(url, download_file)

                full_local_path = os.path.join(local_subfolder, download_file)
                if existing_files: 
                    if full_local_path not in existing_files:
                        try:
                            print('Downloading:', full_download_url)
                            urlretrieve(full_download_url, os.path.join(local_subfolder, download_file))
                            self.relevant_files.append(full_local_path)
                        except Exception as e:
                            print(f"Failed to download {full_download_url}: {e}")
                else:
                    try:
                        print('Downloading:', full_download_url)
                        urlretrieve(full_download_url, os.path.join(local_subfolder, download_file))
                        self.relevant_files.append(full_local_path)
                    except Exception as e:
                        print(f"Failed to download {full_download_url}: {e}")                    

        if len(self.relevant_files) == 0:
            print('Nothing to download.')
    
        print('...Done.')
            
        return self.relevant_files

#########################################################################################################
    
    def check_existing_files(self, input_mission, data_type, local_subfolder, supported_pattern, start_date, end_date):
        # Prepare the date range for searching existing files
        all_dates = [start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)]

        # Get all existing files that match the filename format
        ext = self.get_extension_for_data_type(data_type)
        self.existing_files = [f for f in glob.glob(f'{local_subfolder}*') if re.search(rf'\.{ext}$', f, re.IGNORECASE)]
        
        #print(f'filtered existing: {self.existing_files}')
        if self.existing_files:
            return self.existing_files
#########################################################################################################
    
    def dynamic_download_url_files_time_interval(self, input_mission, local_path, start_date, end_date, url): 
        """
        Downloads files within a specified time interval if they are not already present in the local path.
    
        Parameters:
        - input_mission: Name of the mission (string).
        - local_path: Directory path of mission where subfolders are stored (string).
        - start_date: The start date of the time interval (datetime).
        - end_date: The end date of the time interval (datetime).
        - url: The base URL(s) for downloading files (string).
    
        Returns:
        - List of files that were found or downloaded.
        """

        self.relevant_files = []     
        print('Checking URL:', url)
        data_type = url.split('/')[-2]      
        local_subfolder = local_path + data_type + '/'
    
        # Prepare the date range for searching existing files
        all_dates = [start_date + timedelta(days=x) for x in range((end_date - start_date).days + 1)]

        try: 
            supported_pattern = self.supported_patterns[input_mission][data_type]
        except:
            raise ValueError(f'Pattern not found among supported patterns.')

        existing_files = self.check_existing_files(input_mission, data_type, local_subfolder, supported_pattern, start_date, end_date)
        if existing_files and self.check == False:
            print(f'--------------------------------------- EXISTING FILES CHECK ---------------------------------------------\n')
            print(f'The following files already exist in the folder:\n\n {existing_files}\n\n and will not be downloaded.')
            self.relevant_files.extend(existing_files)

        self.check = True
   
        # Initialize a dictionary to hold files from the HTML response
        files_url_dict = {}  
    
        # Get the content of the URL
        reqs = requests.get(url)
    
        # Parse links from the HTML response
        for link in BeautifulSoup(reqs.text, 'html.parser').find_all('a'):
            full_link = link.get('href')  # Extract the href attribute
            if isinstance(full_link, str):
                # Check if the link is a full URL or a relative path
                if ('/') in full_link:
                    # It's a full URL
                    filename = os.path.basename(full_link)  # Get the filename from the URL
                else:
                    # It's a relative path
                    filename = full_link                
                # Check if the filename matches the specified format
                if self.match_type_extension(data_type, filename):
                    try:
                        # Parse the information from the filename
                        dictionary, underscores = self.parse_filename(input_mission, data_type, filename)
                        # Reconstruct the filename and get time intervals
                        filename_to_download = self.reconstruct_filename(dictionary, underscores)
                        #print(f'filename_to_download: {filename_to_download}')
                        start_time = dictionary["start_date_utc"]
                        end_time = dictionary["end_date_utc"]
    
                        if start_time is None or end_time is None:
                            print(f'Unwanted filename found: {filename_to_download}. Skipping... [Do not worry! ;)]')
                            continue
    
                        # Store the filename in the dictionary with its time interval
                        files_url_dict[(start_time, end_time)] = filename_to_download
                        
                    except:
                        continue  # Skip to the next link

        # Download files for all intervals from the HTML response
        for new_interval, filename_to_download in files_url_dict.items():
            full_local_path = local_subfolder + filename_to_download
            #print(f'full_local_path: {full_local_path}')

            if existing_files:
                if full_local_path not in existing_files:
                    if any(self.is_date_in_intervals(date, [new_interval]) for date in all_dates):
                        print('Downloading', filename_to_download)  # Print which file is being downloaded
                        urlretrieve(url + filename_to_download, full_local_path)  # Download the file
                        self.relevant_files.append(full_local_path)
                        #spice.load_kernel(full_local_path)  # Load the downloaded file
            else:
                if any(self.is_date_in_intervals(date, [new_interval]) for date in all_dates):
                    print('Downloading', filename_to_download)  # Print which file is being downloaded
                    urlretrieve(url + filename_to_download, full_local_path)  # Download the file
                    self.relevant_files.append(full_local_path)
                    #spice.load_kernel(full_local_path)  # Load the downloaded file              

        if len(self.relevant_files) == 0:
            print('Nothing to download.')
    
        print('...Done.')  # Indicate completion

        self.check = False
        return self.relevant_files

#########################################################################################################
    
    def parse_filename(self, input_mission, data_type, filename):
        
        data_type_lower = data_type.lower()
        supported_radio_science_types = ['odf', 'ifms', 'dp2', 'dps', 'dpx']
        supported_spice_kernels_types = ['ck', 'dsk', 'spk', 'fk', 'mk', 'ik', 'lsk', 'pck', 'sclk']
        supported_ancillary_types = ['ion','tro']
        
        all_supported_types = supported_radio_science_types + supported_spice_kernels_types + supported_ancillary_types
        
        
        if input_mission in self.supported_patterns:
            level_one_object = self.supported_patterns[input_mission]
        else:
            raise ValueError('Selected Mission Not Supported (yet!) Aborting ...')   
            
        if data_type_lower == "radioscience":
            data_type_types = supported_radio_science_types
                    # Attempt to retrieve the pattern for the specified mission directly

        elif data_type_lower  == "kernels":
            data_type_types = supported_spice_kernels_types

        elif data_type_lower == "all":
            data_type_types = all_supported_types
            
        elif data_type_lower  in all_supported_types:
            data_type_types = [data_type_lower]

        else:
            print(f'Specified data type: {data_type} not supported.')
        
        """ 
        This function parses the filename into its components and tracks the position of underscores. 
        The underscore position is indexed based on the groups' placement.
        """

        # Initialize the dictionary and underscore index list
        dictionary = {}
        underscore_indices = []
    
        if type(level_one_object) is not str: #deals with multiple keys
            for data_type_type in level_one_object.keys():
                for data_type_type in data_type_types:
                    if self.match_type_extension(data_type_type, filename):
                        pattern = level_one_object[data_type_type]
                        match = re.match(pattern, filename)
    
                        if match:
                            # Populate the dictionary with matched groups
                            dictionary = match.groupdict()
                    
                            # Add underscore index tracking based on group positions
                            last_pos = 0
                            group_index = 0
                            for key in dictionary.keys():
                                current_pos = match.start(key)  # Get the start position of the matched group
                                if dictionary[key] is not None:
                                    if last_pos != current_pos and last_pos != 0:  # Check if there's a gap since last valid group
                                        underscore_indices.append(group_index)  # Track the index of the underscore
                                    last_pos = match.end(key)  # Move to the end of the current matched group
                                group_index += 1  # Increment the group index
                        
                            # If present, convert start and end dates in utc (both will only be present in spice file names).
                            if "start_date_file" in dictionary and "end_date_file" in dictionary:
                                if dictionary['start_date_file'] is not None:
                                    
                                    if dictionary['start_date_file'][0] == 'P': #this deals with MEX ORMF files
                    
                                        self.start_date_utc = LoadPDS.generic_strptime(LoadPDS, dictionary["start_date_file"][1:7]) 
                                        self.end_date_utc = self.start_date_utc + relativedelta(months = 1)
                                        dictionary["start_date_utc"] = self.start_date_utc 
                                        dictionary["end_date_utc"] = self.end_date_utc
                                        
                                    elif (len(dictionary['start_date_file']) == 4 and dictionary['purpose'] == 'SA'):
                                        mex_sa_date_trick = dictionary["start_date_file"] + '0101_'
                                        self.start_date_utc = LoadPDS.generic_strptime(LoadPDS, mex_sa_date_trick) 
                                        self.end_date_utc = self.start_date_utc + relativedelta(years = 1) 
                                        dictionary["start_date_utc"] = self.start_date_utc
                                        dictionary["end_date_utc"] = self.end_date_utc
                                        
                                    else:
                                        self.start_date_utc = LoadPDS.generic_strptime(LoadPDS, dictionary["start_date_file"]) 
                                        self.end_date_utc = LoadPDS.generic_strptime(LoadPDS, dictionary["end_date_file"]) if dictionary["end_date_file"] != '000000' else self.start_date_utc + relativedelta(months = 1) #this deals with MEX ORMM files
                                        dictionary["start_date_utc"] = self.start_date_utc
                                        dictionary["end_date_utc"] = self.end_date_utc

                            # If present, convert date_file in utc (only one date is present in the Radio Science file names)
                            elif "date_file" in dictionary:
                                self.date_utc = LoadPDS.generic_strptime(LoadPDS, dictionary["date_file"])  
                                dictionary["date_utc"] = self.date_utc

                            return dictionary, underscore_indices
              
                        else:
                            #print(f'Filename: {filename} does not match any supported patterns.')
                            continue
                    else:
                        continue
#########################################################################################################

    def reconstruct_filename(self, dictionary, underscore_indices):
        """ 
        This function reconstructs the original filename based on the parsed groups and underscore indices.
        Args:
            dictionary (dict): Dictionary with the parsed components of the filename.
            underscore_indices (list): List of group indices where underscores should be placed.
        Returns:
            str: The reconstructed filename.
        """
        # Extract the keys (group names) from the dictionary in order
        group_values = list(dictionary.values())
        
        # Initialize an empty list to build the filename
        reconstructed_filename = []
        
        # Iterate through the group values and append to the list
        for i, group in enumerate(group_values):
            if isinstance(group, str):
                reconstructed_filename.append(group)  # Add the current group
            
            # Check if the current index is in the underscore_indices list
            if i + 1 in underscore_indices:  # Underscore comes *after* the current group (between i and i+1)
                reconstructed_filename.append("_")
        
        # Join the list into a single string (the reconstructed filename)
        return "".join(reconstructed_filename)

#########################################################################################################

    def get_mission_files(self, input_mission, start_date, end_date, wanted_files = None):
     
        spice.clear_kernels()
        all_dates = [start_date+timedelta(days=x) for x in range((end_date-start_date).days+1)]
        print(f'======================================= Downloading {input_mission.upper()} Data ==============================================\n')   
        print(f'=========================================== Folder(s) Creation ==================================================')
        local_folder = f'{input_mission}_archive/'
        if (os.path.exists(local_folder) == False):
            print(f'Creating Local Folder: {local_folder} and its subfolders (kernels and radio)...')
            os.mkdir(local_folder)
            # Define the subfolders to create
        else:
            print(f'Folder: {local_folder} already exists and will not be overwritten.') 
            
        subfolders = ['ck', 'spk', 'fk', 'sclk', 'lsk', 'ifms', 'odf', 'dp2', 'dps', 'dpx', 'tro', 'ion']     
        # Create each subfolder inside the main folder
        for subfolder in subfolders:
            if (os.path.exists(local_folder + subfolder) == False):
                os.makedirs(os.path.join(local_folder, subfolder))
                print(f'Creating folder: {subfolder}.') 
            else:
                print(f'Folder: {subfolder} already exists and will not be overwritten.') 
        
        print('...Done.')
        print(f'===============================================================================================================\n')
        if input_mission == 'mex':
            self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load = self.get_mex_files(local_folder, start_date, end_date) 
        elif input_mission == 'juice':
            self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load =              self.get_juice_files(local_folder, start_date, end_date) 
        elif input_mission == 'mro':
            self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load = self.get_mro_files(local_folder, start_date, end_date) 

        if self.kernel_files_to_load:
            for kernel_type, kernel_files in self.kernel_files_to_load.items():
                for kernel_file in kernel_files:  # Iterate over each file in the list
                    spice.load_kernel(kernel_file)  # Load each file individually

            n_kernels = spice.get_total_count_of_kernels_loaded()
            std_kernels = spice.load_standard_kernels()
            n_standard_kernels = spice.get_total_count_of_kernels_loaded() - n_kernels
            print(f'===============================================================================================================')
            print(f'Number of Loaded Existing + Downloaded Kernels: {n_kernels}')
            print(f'Number of Loaded Standard Kernels: {n_standard_kernels}')
            print(f'===============================================================================================================')
            
        else:
            print('No Kernel Files to Load.')
            
        return self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load, 
            
########################################################################################################################################
############################################# START OF MEX SECTION #####################################################################
########################################################################################################################################
    
    def get_mex_files(self, local_folder, start_date, end_date):

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {} #empty for now
        
        input_mission = 'mex'
        all_dates = [start_date+timedelta(days=x) for x in range((end_date-start_date).days+1)]
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Radio Science Kernels:')
        url_radio_science_files = self.get_url_mex_radio_science_files(start_date, end_date) 
        for url_radio_science_file_new in url_radio_science_files:
            for closed_loop_type in ['ifms/dp2/', 'dsn/dps/', 'dsn/dpx/']:
                try:
                    url_radio_science_file = url_radio_science_file_new + 'data/level02/closed_loop/' + closed_loop_type
                    files = self.dynamic_download_url_files_single_time(input_mission,
                        local_path=local_folder, start_date=start_date,end_date=end_date,
                        url=url_radio_science_file)
                    key = f"{closed_loop_type.split('/')[0]}_{closed_loop_type.split('/')[1]}"
                    self.radio_science_files_to_load[key] = files
                except:
                    continue
                
        # Clock files
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Clock Kernels:')
        url_clock_files="https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/sclk/"
        wanted_clock_files = ["MEX_241031_STEP.TSC"] #latest mex tsc file
        clock_files_to_load = self.get_kernels(input_mission, url_clock_files, wanted_clock_files, local_folder)
        
        if clock_files_to_load:
            self.kernel_files_to_load['sclk'] = clock_files_to_load
        else:
            print('No sclk files to download this time.')   

        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Frame Kernels:')
        url_frame_files="https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/fk/"
        wanted_frame_files = ["MEX_V16.TF","MEX_SCI_V01.TF", 
                              "MEX_RELAY_LOCATIONS_V03.TF",
                              "MEX_DSK_SURFACES_V04.TF",
                              "MEX_PFS_ROIS_V02.TF"] 
        frame_files_to_load = self.get_kernels(input_mission, url_frame_files, wanted_frame_files, local_folder)

        if frame_files_to_load:
            self.kernel_files_to_load['fk'] = frame_files_to_load
        else:
            print('No fk files to download this time.')   

        # Spk files
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} SPK Kernels:')
        url_spk_files = ["https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/spk/"]
        
        spk_files_to_load = []
        
        if len(url_spk_files) == 1: 
            spk_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
            local_path=local_folder, start_date=start_date,end_date=end_date,
            url=url_spk_files[0])  
            
        else:            
            for url_spk_file in url_spk_files:
                spk_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=url_spk_file)


        if spk_files_to_load:
            self.kernel_files_to_load['spk'] = spk_files_to_load
        else:
            print('No spk files to download this time.')   

    
        # Orientation files
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} CK Kernels:')
        url_ck_files = ["https://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/ck/"]
        
        ck_files_to_load = []
        if len(url_ck_files) == 1:
            ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
            local_path=local_folder, start_date=start_date,end_date=end_date,
            url=url_ck_files[0])    

        else:
            for url_ck_file in url_ck_files:
                ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=url_ck_file)

        if ck_files_to_load:
            self.kernel_files_to_load['ck'] = ck_files_to_load
        else:
            print('No spk files to download this time.')  

        # Tropospheric corrections
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Tropospheric and Ionospheric Corrections Files')
        url_tropo_files = self.get_url_mex_radio_science_files(start_date, end_date) 
        for url_tropo_file_new in url_tropo_files:                
            for folder_type in ['calib/closed_loop/ifms/met/','calib/closed_loop/dsn/ion/', 'calib/closed_loop/dsn/tro/']:
                url_flag = False
                url_tropo_file = url_tropo_file_new + folder_type
                response = requests.get(url_tropo_file)
                if response.status_code == 200:
                    url_flag = True
                    key = folder_type.split('/')[-2]  
                    html = response.text
                    # Parse the HTML with BeautifulSoup
                    soup = BeautifulSoup(html, 'html.parser')
                    # Extract file links and their names
                    wanted_tropo_files = []
                    for link in soup.find_all('a'):
                        href = link.get('href')
                        if href.endswith('.tab') or href.endswith('.aux') :
                            wanted_tropo_files.append(href.split('/')[-1])
                    tropo_files_to_load = self.get_kernels(input_mission, url_tropo_file, wanted_tropo_files, local_folder)
                else:
                    print(f'URL: {url_tropo_file} does not exist.') 
                    continue
                    
                if tropo_files_to_load:
                    key = folder_type.split('/')[-2]
                    self.ancillary_files_to_load[key] = tropo_files_to_load
        
                else:
                    print('No tropospheric or ionospheric files to download this time.') 

        print(f'-----------------------------------------------------------------------------------------------------------')
        print('All requested, relevant and previously non-existing MEX files have been now downloaded. Enjoy!')

        return self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load

#########################################################################################################

    def get_mex_volume_ID(self, start_date, end_date, interval_dict):
        self.volume_id_list = []
        for (key_start, key_end), item in interval_dict.items():
            # Check if either start or end of the input interval overlaps with the dictionary interval
            if (start_date <= key_end and end_date >= key_start):
                self.volume_id_list.append(item['volume_id'])

        if self.volume_id_list:
            return self.volume_id_list
        else:
            raise ValueError(f'No MEX Volume_ID found associated to input interval: {start_date} - {end_date}.')
            
#########################################################################################################
    
    def get_url_mex_radio_science_files(self, start_date_mex, end_date_mex):

        url = "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-v1/mexmrs_0735/aareadme.txt"
        radio_science_base_url = "https://pds-geosciences.wustl.edu/mex/mex-m-mrs-1_2_3-v1/"
        mapping_dict = self.get_mex_volume_ID_mapping(url)

        self.radio_science_urls = []
        volume_ID_list = self.get_mex_volume_ID(start_date_mex, end_date_mex, mapping_dict)
        if self.get_mex_volume_ID(start_date_mex, end_date_mex, mapping_dict):
            for volume_ID in volume_ID_list:
                volume_ID_url = radio_science_base_url + volume_ID + '/'
                self.radio_science_urls.append(volume_ID_url)
                
        if len(self.radio_science_urls) > 0:
            return self.radio_science_urls
        else:
            raise ValueError(f'No url available for MEX radio science files. Please check the mapping.')
            
#########################################################################################################  
    
    def get_mex_volume_ID_mapping(self, url):
        # Step 1: Fetch content from the URL
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors
        aareadme_text = response.text
    
        # Step 2: Parse content using regex to extract the table entries
        self.mapping_dict = {}
        pattern = re.compile(r'^\s*(MEXMRS_\d{4})\s+(\d{4}-\d{2}-\d{2})\s+(\d{4}-\d{2}-\d{2})\s+(.+?)\s*$', re.MULTILINE)
    
        # Step 3: Find all matches and populate the dictionary
        for match in pattern.finditer(aareadme_text):
            volume_id = match.group(1)
            start_date_file = match.group(2)
            end_date_file = match.group(3)
            start_date_utc = LoadPDS.generic_strptime(LoadPDS, match.group(2)) 
            end_date_utc = LoadPDS.generic_strptime(LoadPDS, match.group(3)) 
            interval_key_for_retrieval = (start_date_utc, end_date_utc)
            observation_type = match.group(4).strip()
    
            # Add entry to dictionary
            self.mapping_dict[interval_key_for_retrieval] = {
                "volume_id": volume_id,
                "start_date_file": start_date_file,
                "end_date_file": end_date_file,
                "observation_type": observation_type
            }
    
        return self.mapping_dict

########################################################################################################################################
################################################# END OF MEX SECTION ###############################################################
########################################################################################################################################
    
#--------------------------------------------------------------------------------------------------------------------------------------#
    
########################################################################################################################################
################################################# START OF JUICE SECTION ###############################################################
########################################################################################################################################

    def get_juice_files(self, local_folder, start_date, end_date):

        self.radio_science_files_to_load = {} #empty for now, since we need fdets
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {} #empty for now
        
        input_mission = 'juice'
        all_dates = [start_date+timedelta(days=x) for x in range((end_date-start_date).days+1)]
    
        # Clock files

        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Clock Files:')
        url_clock_files="https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/sclk/"
        wanted_clock_files =["juice_step_20160326_v03.tsc", "juice_fict_160326_v02.tsc"] #latest fict and step tsc files
        clock_files_to_load = self.get_kernels(input_mission, url_clock_files, wanted_clock_files, local_folder)
        
        if clock_files_to_load:
            self.kernel_files_to_load['sclk'] = clock_files_to_load
        else:
            print('No sclk files to download this time.')   

        # Frame Kernels
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Frame Files:')
        url_frame_files="https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/fk/"
        wanted_frame_files=["juice_v41.tf",
                                           "juice_events_crema_5_1_150lb_23_1_v02.tf", 
                                           "juice_ops_v11.tf",
                                           "juice_stations_topo_v01.tf",
                                           "juice_dsk_surfaces_v11.tf",
                                           "juice_sci_v17.tf",
                                           "juice_roi_v02.tf"] 
        
        frame_files_to_load = self.get_kernels(input_mission, url_frame_files, wanted_frame_files, local_folder)
        if frame_files_to_load:
            self.kernel_files_to_load['fk'] = frame_files_to_load
        else:
            print('No fk files to download this time.')   

        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Orientation Kernels:')
        ck_files_to_load = []
        wanted_ck_files=["juice_mga_crema_5_1_150lb_23_1_baseline_v04.bc",
                                           "juice_sa_crema_5_1_150lb_23_1_baseline_v04.bc", 
                                           "juice_sc_crema_5_1_150lb_23_1_baseline_v03.bc",] 
        
        url_planned_ck_files="https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/"
        planned_ck_files_to_load = self.get_kernels(input_mission, url_planned_ck_files, wanted_ck_files, local_folder)
        
        if planned_ck_files_to_load:
            for file in planned_ck_files_to_load:
                ck_files_to_load.append(file)
        else:
            print('No Crema Planned CK files to download this time.')   
    

        measured_url_ck_files = ["https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/ck/"]
        if len(measured_url_ck_files) == 1:
            measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_ck_files[0])
        else:
            for measured_url_ck_file in measured_url_ck_files:
                measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_ck_file)    

        if measured_ck_files_to_load:
            for file in measured_ck_files_to_load:
                ck_files_to_load.append(file)
        else:
            print('No Measured ck files to download this time.')   

        if len(ck_files_to_load) > 0:
            self.kernel_files_to_load['ck'] = ck_files_to_load
        else:
            print('No Overall CK files to download this time.')  

        # SPK files
        spk_files_to_load = []
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} SPK Kernels:')
        measured_url_spk_files = ["https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/"]

        if len(measured_url_spk_files) == 1:
            measured_spk_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_spk_files[0])
        else:
            for measured_url_spk_file in measured_url_spk_files:
                measured_spk_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_spk_file)    

        if measured_spk_files_to_load:
            for file in measured_spk_files_to_load:
                spk_files_to_load.append(file)
        else:
            print('No spk files to download this time.')   

        print('---------------------------------------------')

        if len(spk_files_to_load) > 0:
            self.kernel_files_to_load['spk'] = spk_files_to_load
        else:
            print('No Overall SPK files to download this time.')      


        # Tropospheric corrections
        #print(f'===========================================================================================================\n') 
        #print(f'Download {input_mission.upper()} tropospheric corrections files\n')
        #url_tropo_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/tro/"
        #tropo_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
        #    local_path=local_folder, start_date=start_date,
        #    end_date=end_date, url=url_tropo_files)
    
        #if tropo_files_to_load:
        #    self.ancillary_files_to_load['tro'] = tropo_files_to_load
        #else:
        #    print('No tropospheric files to download this time.') 
            
        # Ionospheric corrections
        #print(f'===========================================================================================================\n') 
        #print(f'Download {input_mission.upper()} ionospheric corrections files\n')
        #url_ion_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/ion/"
        #ion_files_to_load = self.dynamic_download_url_files_time_interval(input_mission, 
        #    local_path=local_folder, start_date=start_date, end_date=end_date, 
        #    url=url_ion_files)

        #if ion_files_to_load:
        #    self.ancillary_files_to_load['ion'] = ion_files_to_load
        #else:
        #    print('No ionospheric files to download this time.') 
        
        return self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load

########################################################################################################################################
################################################### END OF JUICE SECTION ###############################################################
########################################################################################################################################
    
#--------------------------------------------------------------------------------------------------------------------------------------#
    
########################################################################################################################################
################################################### START OF MRO SECTION ###############################################################
########################################################################################################################################

    def get_mro_files(self,local_folder, start_date, end_date):

        self.radio_science_files_to_load = {}
        self.kernel_files_to_load = {}
        self.ancillary_files_to_load = {}
        
        input_mission = 'mro'
        all_dates = [start_date+timedelta(days=x) for x in range((end_date-start_date).days+1)]

        # ODF files
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} ODF files:')
        url_radio_science_files = ["https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/odf/"]
        for url_radio_science_file in url_radio_science_files:
                try:
                    files = self.dynamic_download_url_files_single_time(input_mission,
                        local_path=local_folder, start_date=start_date,end_date=end_date,
                        url=url_radio_science_file)
                    key = "odf"
                    self.radio_science_files_to_load[key] = files
                except:
                    continue
        
        if not self.radio_science_files_to_load:
            print('No Radio Science files to download this time.')
        
        #Clock Kernels
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Clock Kernels:')
        url_clock_files = "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/sclk/"
        wanted_clock_files= ["mro_sclkscet_00112_65536.tsc"]
        clock_files_to_load = self.get_kernels(input_mission, url_clock_files, wanted_clock_files, local_folder)

        if clock_files_to_load:
            self.kernel_files_to_load['sclk'] = clock_files_to_load
        else:
            print('No sclk files to download this time.')   
            
        #Frame Kernels
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Frame Kernels:')
        url_frame_files="https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/fk/"
        wanted_frame_files=["mro_v16.tf"] 
        frame_files_to_load = self.get_kernels(input_mission, url_frame_files, wanted_frame_files, local_folder)

        if frame_files_to_load:
            self.kernel_files_to_load['fk'] = frame_files_to_load
        else:
            print('No fk files to download this time.')   
            
        # Planetary and Ephemeris Kernels
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} SPK Kernels:')
        url_spk_files ="https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/spk/"
        wanted_spk_files = self.get_url_mro_spk_files(start_date, end_date)
        spk_files_to_load = self.get_kernels(input_mission, url_spk_files, wanted_spk_files, local_folder)
        wanted_struct_files = ["mro_struct_v10.bsp"]
        struct_files_to_load = self.get_kernels(input_mission, url_spk_files, wanted_struct_files, local_folder)  
        spk_files_to_load.extend(struct_files_to_load)
        
        if spk_files_to_load:
            self.kernel_files_to_load['spk'] = spk_files_to_load
        else:
            print('No spk files to download this time.')    
                   
        # Orientation Kernels
        ck_files_to_load = []
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Orientation Kernels:')
        measured_url_ck_files =["https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/ck/"]

        if len(measured_url_ck_files) == 1:
            measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_ck_files[0])
        else:
            for measured_url_ck_file in measured_url_ck_files:
                measured_ck_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
                local_path=local_folder, start_date=start_date,end_date=end_date,
                url=measured_url_ck_file)   

        if measured_ck_files_to_load:
            self.kernel_files_to_load['ck'] = measured_ck_files_to_load
        else:
            print('No ck files to download this time.')   
              
        # Tropospheric corrections
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Tropospheric Corrections Files')
        url_tropo_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/tro/"
        tropo_files_to_load = self.dynamic_download_url_files_time_interval(input_mission,
            local_path=local_folder, start_date=start_date,
            end_date=end_date, url=url_tropo_files)
    
        if tropo_files_to_load:
            self.ancillary_files_to_load['tro'] = tropo_files_to_load
        else:
            print('No tropospheric files to download this time.') 
            
        # Ionospheric corrections
        print(f'===========================================================================================================') 
        print(f'Download {input_mission.upper()} Ionospheric Corrections Files')
        url_ion_files = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/mrors_0xxx/ancillary/ion/"
        ion_files_to_load = self.dynamic_download_url_files_time_interval(input_mission, 
            local_path=local_folder, start_date=start_date, end_date=end_date, 
            url=url_ion_files)

        if ion_files_to_load:
            self.ancillary_files_to_load['ion'] = ion_files_to_load
        else:
            print('No ionospheric files to download this time.') 
        
        return self.kernel_files_to_load, self.radio_science_files_to_load, self.ancillary_files_to_load

########################################################################################################################################

    def get_mapping_dict_spk_mro(self,url):
        # Step 1: Fetch content from the URL
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors
        file_text = response.text
    
        # Step 2: Define a dictionary to hold date intervals and phase names
        self.phase_dict = {}
    
        # Step 3: Regex pattern to capture phase name and dates
        pattern = re.compile(
            r'^\s*(\w+)\s+(\d{4} \w{3} \d{2})\s+(\d{4} \w{3} \d{2})\s*$',
            re.MULTILINE
        )
    
        # Step 4: Extract all matches and populate the dictionary
        for match in pattern.finditer(file_text):
            phase = match.group(1)
            start_date_str = match.group(2)
            end_date_str = match.group(3)
    
            # Convert string dates to datetime objects for easier manipulation if needed
            start_date = datetime.strptime(start_date_str, '%Y %b %d')
            end_date = datetime.strptime(end_date_str, '%Y %b %d')
    
            # Use the date interval as the key and the phase as the value
            self.phase_dict[(start_date, end_date)] = {'spk_ID': phase}
    
        return self.phase_dict

########################################################################################################################################

    def get_mro_spk_ID(self, start_date, end_date, mapping_dict):
        self.spk_id_list = []
        for (key_start, key_end), item in mapping_dict.items():
            # Check if either start or end of the input interval overlaps with the dictionary interval
            if (start_date <= key_end and end_date >= key_start):
                self.spk_id_list.append(item['spk_ID'])

        if self.spk_id_list:
            return self.spk_id_list
        else:
            raise ValueError(f'No MRO SPK_ID found associated to input interval: {start_date} - {end_date}.')
            
########################################################################################################################################

    def get_url_mro_spk_files(self, start_date, end_date):

        url = "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/mrosp_1000/data/spk/spkinfo.txt"
        mapping_dict = self.get_mapping_dict_spk_mro(url)
        
        self.spk_ID_urls = []
        spk_ID_list = self.get_mro_spk_ID(start_date, end_date, mapping_dict)
        if spk_ID_list:
            for spk_ID in spk_ID_list:
                spk_ID_url = f"mro_{spk_ID}.bsp"
                self.spk_ID_urls.append(spk_ID_url)
                
        if len(self.spk_ID_urls) > 0:
            return self.spk_ID_urls
        else:
            raise ValueError(f'No url available for MEX radio science files. Please check the mapping.')
########################################################################################################################################
##################################################### END OF MRO SECTION ###############################################################
########################################################################################################################################
    