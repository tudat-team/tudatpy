import os
import sys
sys.path.insert(0, '/home/fdahmani/tudatcompile/tudat-bundle/cmake-build-release/tudatpy')
import inspect
import csv
import numpy as np
import datetime
from matplotlib import pyplot as plt
from math import *
from astropy.time import Time
from astropy.coordinates import SkyCoord, ICRS,FK4, FK5

from tudatpy.io import save2txt
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation  
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import estimation, estimation_setup
from tudatpy.kernel.numerical_simulation.estimation_setup import observation
#from HorizonsLoad import jpl_horizons_ephemeris_table 
from tudatpy.kernel.astro import time_conversion, element_conversion


# Utility functions
def moon_ID (moon, central_body):     #Numbers, e.g. 5, 5 for Amalthea
    '''
    Converts a combination of moon and central body into a spice-recognizable ID

    Parameters
    -----------
    moon: str, float or int
        number of moon around the central body
    central_body: str, float or int
        number of planet around sun 


    Examples
    -----------
    .. code-block:: python
        # Define required data
        central_body = 5    #Jupiter
        moon = 5            #Amalthea

        body = moon_ID(moon, central_body)


    Returns
    -----------
    Body_ID: str
        spice recognisable number of the body in a string

    '''
    if central_body.isalpha():
        if central_body.lower() == 'jupiter':
            central_body = 5
        elif central_body.lower() == 'saturn':
            central_body = 6
        elif central_body.lower() == 'uranus':
            central_body = 7
        elif central_body.lower() == 'neptune':
            central_body = 8        
        else:
            print('Central body not found, please specify by number')
    return str(100*int(central_body)+int(moon))

def observatory_info(Observatory):
    '''
    Retrieves observatory info from the observatory code

    Parameters
    -----------
    Observatory: str
        observatory identifier


    Examples
    -----------
    .. code-block:: python
        # Define required data
        Observatory = 'B18'
        
        longitude, latitude, altitude = observatory_info(Observatory)


    Returns
    -----------
    longitude: float
        longitude of the observatory [rad]
    latitude: float
        latitude of the observatory [rad]
    altitude: float
        altitude of the observatory [m]
        
    '''
    if len(Observatory) == 2:                   #Making sure 098 and 98 are the same
        Observatory = '0' + Observatory
    elif len(Observatory) == 1:                   #Making sure 098 and 98 are the same
        Observatory = '00' + Observatory
    with open('Observatories.txt', 'r') as file:    #https://www.projectpluto.com/obsc.htm, https://www.projectpluto.com/mpc_stat.txt
        lines = file.readlines()
        for line in lines[1:]:  # Ignore the first line
            columns = line.split()
            if columns[1] == Observatory:
                longitude = float(columns[2])
                latitude = float(columns[3])
                altitude = float(columns[4])
                return np.deg2rad(longitude),  np.deg2rad(latitude), altitude
        print('No matching Observatory found')

def read_file_settings(filename):
    '''
    Read the meta data out of the first lines of the file

    Parameters
    -----------
    filename: str
        name of the nsdc datafile to read


    Examples
    -----------
    .. code-block:: python
        # Define required data
        filename = 'ji0001.txt'    

        meta_data = read_file_settings(filename)


    Returns
    -----------
    data_type: str
        describes type of data: absolute, separation arc, relative tangential
    moon_entry: str or int
        the entry of the target body for every observation. Can also be a string with a moon's name if only 1 moon is observed in the set
    time_entry: int
        the first entry of the time of observation
    amount_time_entries: int
        the amount of time entries per observation, necessary in processing. IMPORTANT: assumes year-month-day hour-minute-second framework. Different order would require manual fixes
    time_type: str
        time format: JD, MJD, J2000, DATE, MINsinceJDaddingtheJD
    time_scale: str
        time standards: utc, ut1, ...
    time_delta: float
        offset between the time standard and the reference epoch for the observation set
    RA_entry: int
        the first entry describing the right ascension of the observed body
    DEC_entry: int
        the first entry describing declination of the observed body
    number_observatories: int
        amount of observatories employed in the dataset. Mostly relevant for the difference between 1 or multiple
    observatory: str
        only relevant if only 1 observatory: the ID of the observatory
    observatory_entrace: int
        only relevant if multiple observatories: the entry of the observatory ID
    relative_body_entry: str
        the entry of the relative body for every observation. Can also be a string with a body's name if only 1 relative body is used in the set
    process_relative: bool
        switch whether to process the relative observations to absolute
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...
    central_body: str
        identifier or name of the central body

    '''
    first_words = []
    with open(filename, 'r') as datafile:
        #Extracting meta data
        for reading in range(16):
            first_word = datafile.readline().strip().split()[0]
            first_words.append(first_word)
    
    #Coupling meta data to correct variables
    data_type = first_words[0]              #ABS (absolute), SEP (separation), REL (tangential), DIF (differential coorinates) or CPV (camera pointing vector)
    moon_entry = first_words[1]                     #Target body
    time_entry = int(first_words[2])                       #First time entry
    amount_time_entries = int(first_words[3])              #Amount of time entries
    time_type = first_words[4]                  #JD, MJD, J2000, DATE, MINsinceJDaddingtheJD
    time_scale = first_words[5]                  #time scale
    time_delta = float(first_words[6])                     #hours away from standard time (i.e. utc-3 is -3)
    RA_entry = int(first_words[7])                         #First RA or X entry
    DEC_entry = int(first_words[8])                        #First DEC or Y entry
    number_observatories = int(first_words[9])             #Amount of observatories, needs extra definition which one if 1
    observatory = first_words[10]                 #Ignored if more than 1 observatory
    observatory_entrace = int(first_words[11])              #Ignored if 1 observatory
    relative_body_entry = first_words[12]    #Ignored if ABS, can either be a name or entry position
    process_relative = bool(first_words[13])              #Switch for processing relative observations or not
    orientation = first_words[14]
    central_body =  first_words[15]

    return (data_type,moon_entry,time_entry,amount_time_entries,time_type,time_scale,time_delta,RA_entry,DEC_entry,number_observatories,observatory,observatory_entrace,relative_body_entry,process_relative,orientation,central_body)

def read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body):    #retrieving data from line
    '''
    Retrieve data from a single observation/line in the observation set

    Parameters
    -----------
    observation: str
        one line of the datafile, describing 1 observation
    number_observatories: int
        amount of observatories employed in the dataset. Mostly relevant for the difference between 1 or multiple
    observatory: str
        only relevant if only 1 observatory: the ID of the observatory. If more than 1 observatory, it is overwritten in the function
    observatory_entrace: int
        only relevant if multiple observatories: the entry of the observatory ID
    moon_entry: str or int
        the entry of the target body for every observation. Can also be a string with a moon's name if only 1 moon is observed in the set
    central_body: str
        identifier or name of the central body


    Examples
    -----------
    .. code-block:: python
        with open(filename, 'r') as datafile:
            # Create lists to store the data
            #Skip first 14 lines of metadata
            for skips in range(15):
                next(datafile)

            # Read each remaining line in the file
            for observation in datafile:    
                body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)


    Returns
    -----------
    body: str
        the target body either by name or by code
    observatory: str
        the ID of the observatory
    values: list[float]
        a list with all the entries in the datafile line, describing the observation
    

    '''
    entries = observation.split()
    values = [entry if not (entry.replace(".", "", 1).isdigit() or ((entry[0] == '-' or entry[0] == '+') and entry[1:].replace(".", "", 1).isdigit())) else float(entry) for entry in entries]
    if number_observatories != 1:
        observatory = values[observatory_entrace]
        if type(observatory) == float:
            observatory = str(int(observatory))
    if moon_entry.isalpha():
        body = moon_entry
    else:
        body = moon_ID(values[int(moon_entry)],central_body)  

    return body, observatory, values

def save_to_csv(times, RADEC, moons, observatories, studyname, diflist, output_folder):
    '''
    Sorts observations by observed targed and observatory and writes them to a csv per unique combination. Based on seperate lists where the ith entry of every list corresponds to the same observation
    Saves the time and values of every observation

    Parameters
    -----------
    times: list[float]
        list of observation epochs in seconds since J2000
    RADEC: list[list[float]]
        nested list of 2 entries per observation containing the RA and DEC of the observed body
    moons: list[str]
        observed target bodies
    observatories: list[str]
        observatory used for the observation
    studyname: str
        name of the study or equivalent identifier to retrace the data's origin
    diflist: list[list[float]]
        List of differences between observation as spice, which can be used as noise level later
    output_folder: str
        folder to write the processed csv's to

    Examples
    -----------
    .. code-block:: python
        with open(filename, 'r') as datafile:
            studyname = filename.split("\\")[-1][:-4]
            # Create lists to store the data
            data, times, RADEC, moons,observatories = [],[],[],[],[]

            #Skip first 16 lines of metadata
            for skips in range(16):
                next(datafile)

            # Read each remaining line in the file
            for observation in datafile:    
                body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)
                epoch = return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type)
                RA_deg, DEC_deg = return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body)
                RA_EJ2000, DEC_EJ2000 = return_standardised_rotation(RA_deg,DEC_deg,epoch,orientation)

                ### Saving extracted values
                data.append(values)
                times.append(epoch.epoch())
                RADEC.append([np.rad2deg(RA_EJ2000),np.rad2deg(DEC_EJ2000)])
                moons.append(body)
                observatories.append(observatory)

        save_to_csv(times, RADEC, moons, observatories, studyname)


    Returns
    -----------
    /

    '''
    #Sort by moon and observatory
    unique_moons = list(set(moons))

    #loop over moons and select list of only that moon
    for moon in unique_moons:
        indices = [i for i, x in enumerate(moons) if x == moon]
        times_permoon = [times[i] for i in indices]
        RADEC_permoon = [RADEC[i] for i in indices]
        diflist_permoon = [diflist[i] for i in indices]
        
        observatories_permoon = [observatories[i] for i in indices]

        #Sort by observatory
        unique_observatories = list(set(observatories_permoon))

        #Loop over observatories and select list of only that observatory
        for observatory in unique_observatories:
            indices2 = [counter for counter, y in enumerate(observatories_permoon) if y == observatory]
            epochs = [times_permoon[counter] for counter in indices2]
            observations = [RADEC_permoon[counter] for counter in indices2]         
            difs = [diflist_permoon[counter] for counter in indices2]

            #Repackage into csv file with correct name and time/observations
            with open( output_folder + str(int(moon)) + '_' + str(observatory) + '_' + studyname + '.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header row
                writer.writerow(['Seconds since J2000', 'RA ECLIPJ2000 [rad]','DEC ECLIPJ2000 [rad]','O-C RA ECLIPJ2000 [rad]','O-C DEC ECLIPJ2000 [rad]'])
                
                # Write each matching name and result to a new row
                for j in range(len(epochs)):
                    writer.writerow([epochs[j], observations[j][0],observations[j][1],difs[j][0],difs[j][1]])

def create_observation_set(file,bodies):
    '''
    Reads a csv file and converts it into a tudat-compatible observation set

    Parameters
    -----------
    file: str
        csv file name containing a sorted observation by moon and observatory
    bodies: obj
        system of bodies/environment relevant to the observation set 


    Examples
    -----------
    .. code-block:: python
        bodies => see tudatpy documentation of bodies setup
        observation_set_list = []
        observation_settings_list = []
        raw_observation_files = LIST OF FILES

        #Start loop over every csv in folder
        for file in raw_observation_files:
            observation_set, observation_settings = create_observation_set(file,bodies)
            observation_set_list.append(observation_set)
            observation_settings_list.append(observation_settings)
        observations = estimation.ObservationCollection( observation_set_list )

    Returns
    -----------
    observation_set: obj
        a single tudat-compatible observation set
    observation_settings: obj
        describes the (type of the) corresponding observation 

    '''
    #Reading information from file name
    string_to_split = file.split("\\")[-1]
    split_string = string_to_split.split("_")

    #Extracting moon and observatory numbers from csv title
    Moon = int(split_string[0])
    Observatory = split_string[1]

    #Reading observational data from file
    angles = []
    times = []
    with open(file, 'r') as f:
        csv_reader = csv.reader(f)
        next(csv_reader)
        for row in csv_reader:
            times.append(float(row[0]))
            angles.append(np.asarray([float(row[1]),float(row[2])]))

    # Define the position of the observatory on Earth
    observatory_altitude, observatory_latitude, observatory_longitude = observatory_info (Observatory)

    # Add the ground station to the environment
    environment_setup.add_ground_station(
        bodies.get_body("Earth"),
        "Observatory",
        [observatory_altitude, observatory_latitude, observatory_longitude],
        element_conversion.geodetic_position_type)
    
    # Define link ends
    link_ends = dict()                  
    link_ends[observation.transmitter] = observation.body_origin_link_end_id(Moon)
    link_ends[observation.receiver] = observation.body_reference_point_link_end_id("Earth", "Observatory")
    link_definition = observation.LinkDefinition(link_ends)

    # Create observation set 
    observation_set = estimation.single_observation_set(
        observation.angular_position_type, 
        link_definition,
        angles,
        times, 
        observation.receiver 
    )

    # Create observation settings for each link/observable                                                                                          
    observation_settings = observation.angular_position(link_definition)

    return observation_set, observation_settings

def analysis(times,observations, moons, observatories, bodies, plotting,standard_orientation):
    '''
    Analysing the processed observations

    Parameters
    -----------
    times: list[float]
        list of observation epochs in seconds since J2000
    observations: list[list[float]]
        nested list of 2 entries per observation containing the RA and DEC of the observed body
    moons: list[str]
        observed target bodies
    observatories: list[str]
        observatory used for the observation
    bodies: obj
        system of bodies/environment relevant to the observation set 
    plotting: bool
        whether the user also wants a plot of the difference between spice and observed data as part of the output
    standard_orientation: str
        orientation which will be used for operations, and to which all data will be converted

    Examples
    -----------
    .. code-block:: python
        bodies = load_standard_jupiter_EJ2000()
        means, stddevs, diflist = analysis(times, RADEC, moons, observatories, bodies, True, 'ECLIPJ2000')
        


    Returns
    -----------
    means: list[float]
        mean difference between spice and observed for RA and DEC
    stddevs: list[float]
        standard deviation of difference between spice and observed for RA and DEC
    diflist: list[list[float]]
        nested lists each containing the difference between spice and observed for RA and DEC for every observation
    '''
    diflist = []
    for i in range(len(moons)):
        observatory_ephemerides = environment_setup.create_ground_station_ephemeris(
        bodies.get_body("Earth"),
        'obs' + str(observatories[i]),
        )
        RA_spice, DEC_spice = get_angle_rel_body(time_conversion.date_time_from_epoch(times[i]),'ECLIPJ2000',observatory_ephemerides, moons[i],standard_orientation)
        RA, DEC = observations[i]
        diflist.append([RA-RA_spice,DEC-DEC_spice])
    means = np.mean(diflist, axis=0)
    stddevs = np.std(diflist, axis=0)
    
    
    if plotting == True:
        plt.scatter(times,np.asarray(diflist)[:,0])
        plt.xlabel('Observation epoch [s since J2000]')
        plt.ylabel('spice-observed RA [rad]')
        plt.grid()
        plt.show()

        plt.scatter(times,np.asarray(diflist)[:,1])
        plt.xlabel('Observation epoch [s since J2000]')
        plt.ylabel('spice-observed DEC [rad]')
        plt.grid()
        plt.show()       
    return means, stddevs, diflist 

def remove_outliers(times, observations, moons, observatories, bodies,standard_orientation, outlier_limit = 3, plotting= True):
    '''
    Removing outliers from the dataset

    Parameters
    -----------
    times: list[float]
        list of observation epochs in seconds since J2000
    observations: list[list[float]]
        nested list of 2 entries per observation containing the RA and DEC of the observed body
    moons: list[str]
        observed target bodies
    observatories: list[str]
        observatory used for the observation
    bodies: obj
        system of bodies/environment relevant to the observation set 
    standard_orientation: str
        orientation which will be used for operations, and to which all data will be converted
    outlier_limit: float, optional
        amount of standard deviations points can be away from mean to be considered outliers
    plotting: bool
        whether the user also wants a plot of the difference between spice and observed data as part of the output


    Examples
    -----------
    .. code-block:: python
        bodies = load_standard_jupiter_EJ2000()
        filtered_times, filtered_observations, filtered_moons, filtered_observatories, filtered_diflist, rms, means, stddevs = remove_outliers(times, observations, moons, observatories, bodies, 'ECLPIJ2000', 3, True)


    Returns
    -----------
    filtered_times: list[float]
        list of observation epochs in seconds since J2000 with outliers removed
    filtered_observations: list[list[float]]
        nested list of 2 entries per observation containing the RA and DEC of the observed body  with outliers removed
    filtered_moons: list[str]
        observed target bodies with outliers removed
    filtered_observatories: list[str]
        observatory used for the observation with outliers removed
    filtered_diflist: list[list[float]]
        nested lists each containing the difference between spice and observed for RA and DEC for every observation with outliers removed
    rms: float
        root mean squared for difference between spice and observed for RA and DEC, highest of RA and DEC  
    means: list[float]
        mean difference between spice and observed for RA and DEC
    stddevs: list[float]
        standard deviation of difference between spice and observed for RA and DEC
    '''
    
    means, stddevs, diflist = analysis(times, observations, moons, observatories, bodies, False,standard_orientation)

    outlier_indices = [i for i, nested_list in enumerate(diflist) if any(abs((x - mean) / std_dev) > outlier_limit for x, mean, std_dev in zip(nested_list, means, stddevs))]
    filtered_diflist = [x for i, x in enumerate(diflist) if i not in outlier_indices]
    filtered_observations = [x for i, x in enumerate(observations) if i not in outlier_indices]
    filtered_times = [x for i, x in enumerate(times) if i not in outlier_indices]
    filtered_moons = [x for i, x in enumerate(moons) if i not in outlier_indices]
    filtered_observatories = [x for i, x in enumerate(observatories) if i not in outlier_indices]
    rms = [np.sqrt(np.mean(np.square(np.asarray(filtered_diflist)[:,0]))),np.sqrt(np.mean(np.square(np.asarray(filtered_diflist)[:,1])))]

    print('removed', len(outlier_indices), 'outliers, with indices:', outlier_indices)

    means, stddevs, diflist = analysis(filtered_times, filtered_observations, filtered_moons, filtered_observatories, bodies, False,standard_orientation)
    print('Means:', means, 'Stddevs:', stddevs)

    if plotting == True:
        plt.scatter(filtered_times,np.asarray(filtered_diflist)[:,0])
        plt.xlabel('Observation epoch [s since J2000]')
        plt.ylabel('spice-observed RA [rad]')
        plt.grid()
        plt.show()

        plt.scatter(filtered_times,np.asarray(filtered_diflist)[:,1])
        plt.xlabel('Observation epoch [s since J2000]')
        plt.ylabel('spice-observed DEC [rad]')
        plt.grid()
        plt.show()    
    return filtered_times, filtered_observations, filtered_moons, filtered_observatories, filtered_diflist, max(rms), means, stddevs

#Rotation functions
def J2000_2_ECLIPJ2000(RA_J2000,DEC_J2000): 
    '''
    Rotates a set of RA and DEC from J2000 to ECLIPJ2000

    Parameters
    -----------
    RA_J2000: float
        right ascension in the J2000 reference frame in rad
    DEC_J2000: float
        declination in the J2000 reference frame in rad


    Examples
    -----------
    .. code-block:: python
        RA_EJ2000, DEC_EJ2000 = J2000_2_ECLIPJ2000(RA_J2000,DEC_J2000)


    Returns
    -----------
    RA_ECLIPJ2000: float
        right ascension in the ECLIPJ2000 reference frame in rad
    DEC_ECLIPJ2000: float
        declination in the ECLIPJ2000 reference frame in rad

    '''
    time = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime.datetime(2000,1,1,12,0,0)))   #Dummy values, no effect but required input
    A = spice.compute_rotation_matrix_between_frames("J2000","ECLIPJ2000", time)
    #obliquity = np.deg2rad(23.439292)
    
    #Convert J2000 to Unit Cartesian
    X = np.cos(DEC_J2000) * np.cos(RA_J2000)
    Y = np.cos(DEC_J2000) * np.sin(RA_J2000)
    Z = np.sin(DEC_J2000)
    
    # #Rotate Unit Cartesian to ecliptic cartesian
    # Xecl = X
    # Yecl = Y * cos(obliquity) + Z * sin(obliquity)
    # Zecl = -Y * sin(obliquity) + Z * cos(obliquity)

    Xecl,Yecl,Zecl = A@np.array([X,Y,Z])

    # #Convert ecliptic cartesian to ECLIPJ2000
    RA_ECLIPJ2000 = np.arctan2(Yecl, Xecl)
    DEC_ECLIPJ2000 = np.arcsin(Zecl)

    return RA_ECLIPJ2000, DEC_ECLIPJ2000

def ECLIPJ2000_2_J2000(RA_ECLIPJ2000,DEC_ECLIPJ2000):
    '''
    Rotates a set of RA and DEC from J2000 to ECLIPJ2000

    Parameters
    -----------
    RA_ECLIPJ2000: float
        right ascension in the EJ2000 reference frame in rad
    DEC_ECLIPJ2000: float
        declination in the EJ2000 reference frame in rad


    Examples
    -----------
    .. code-block:: python
        RA_J2000,DEC_J2000 = ECLIPJ2000_2_J2000(RA_EJ2000, DEC_EJ2000)


    Returns
    -----------
    RA_J2000: float
        right ascension in the J2000 reference frame in rad
    DEC_J2000: float
        declination in the J2000 reference frame in rad

    '''
    time = time_conversion.julian_day_to_seconds_since_epoch(time_conversion.calendar_date_to_julian_day(datetime.datetime(2000,1,1,12,0,0)))   #Dummy values, no effect but required input
    A = spice.compute_rotation_matrix_between_frames("ECLIPJ2000","J2000", time)
    #Convert ECLIPJ2000 to Unit Cartesian in the ecliptic
    Xecl = np.cos(DEC_ECLIPJ2000) * np.cos(RA_ECLIPJ2000)
    Yecl = np.cos(DEC_ECLIPJ2000) * np.sin(RA_ECLIPJ2000)
    Zecl = np.sin(DEC_ECLIPJ2000)

    #Rotate ecliptic cartesian to cartesian in the equitoreal
    # X = Xecl
    # Y = Yecl * cos(-obliquity) + Zecl * sin(-obliquity)
    # Z = Yecl * -sin(-obliquity) + Zecl * cos(-obliquity)

    X,Y,Z = A@np.array([Xecl,Yecl,Zecl])
    #Convert ecliptic cartesian to ECLIPJ2000
    RA_J2000 = np.arctan2(Y, X)
    DEC_J2000 = np.arcsin(Z)

    return RA_J2000, DEC_J2000

def Byear_2_J2000(RA,DEC,epoch,orientation):            #Converts a different orientation (e.g. B1950) to J2000 (actually ICRS)
    '''
    Rotates a set of RA and DEC from a specified orientation to J2000

    Parameters
    -----------
    RA_: float
        right ascension in the specified reference frame in rad
    DEC: float
        declination in the specified reference frame in rad
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...


    Examples
    -----------
    .. code-block:: python
        RA_J2000, DEC_J2000 = Byear_2_J2000(RA_deg,DEC_deg,epoch,orientation)


    Returns
    -----------
    RA_J2000: float
        right ascension in the J2000 reference frame in rad
    DEC_J2000: float
        declination in the J2000 reference frame in rad

    '''
    if epoch.year<1970:
        frame = 'fk4'
    else:
        frame = 'fk5'

    if orientation == 'date':
        equinox = str(epoch.year) + '-' + str(epoch.month) + '-' + str(epoch.day)
    elif orientation == 'year':
        equinox = 'B'+ str(epoch.year)
    else:
        equinox = orientation
    coord = SkyCoord(ra=RA, dec=DEC,unit='rad', obstime=Time(epoch.julian_day(),format ='jd'), frame=frame, equinox=equinox)
    coord_j2000 = coord.transform_to(FK5(equinox='J2000'))
    # Extract the J2000.0 coordinates
    RA_J2000 = coord_j2000.ra.rad  # Right Ascension in J2000.0
    DEC_J2000 = coord_j2000.dec.rad  # Declination in J2000.0
    return RA_J2000, DEC_J2000


#Time functions
def JD_to_tudat(values,time_entry,time_scale):
    '''
    Convert a julian date into a tudat-compatible time object

    Parameters
    -----------
    values: list[float]
        entries in one line (one observation) from the dataset
    time_entry: int
        the first entry of the time of observation
    time_scale: str
        time standards: utc, ut1, ...


    Examples
    -----------
    .. code-block:: python
        if time_type == 'JD':
            epoch = JD_to_tudat(values,time_entry,time_scale) 


    Returns
    -----------
    epoch: obj
        tudat-compatible time object

    '''
    epoch = time_conversion.datetime_to_tudat(time_conversion.julian_day_to_calendar_date(Time(values[time_entry],scale=time_scale,format='jd').tdb.value))
    return(epoch)

def date_to_tudat(values,time_entry,time_scale,amount_time_entries,time_delta):
    '''
    Convert a calendar date into a tudat-compatible time object

    Parameters
    -----------
    values: list[float]
        entries in one line (one observation) from the dataset
    time_entry: int
        the first entry of the time of observation
    time_scale: str
        time standards: utc, ut1, ...
    amount_time_entries: int
        the amount of time entries per observation, necessary in processing. IMPORTANT: assumes year-month-day hour-minute-second framework. Different order would require manual fixes
    time_delta: float
        offset between the time standard and the reference epoch for the observation set


    Examples
    -----------
    .. code-block:: python
        if time_type == 'DATE': 
            epoch = date_to_tudat(values,time_entry,time_scale,amount_time_entries,time_delta) 


    Returns
    -----------
    epoch: obj
        tudat-compatible time object

    '''
    if time_delta != 0 and int(amount_time_entries)>3:  
        if values[time_entry+3] - time_delta >= 24:
            values[time_entry+3] = values[time_entry+3] - 24
            values[time_entry+2] = values[time_entry+2] + 1

        elif values[time_entry+3] - time_delta < 0:
            values[time_entry+3] = values[time_entry+3] + 24
            values[time_entry+2] = values[time_entry+2] - 1

    if int(amount_time_entries)==1:
        print('1 date entry not supported, check the correct conversion to days')                 #1 and 2 not fully done yet like below
    elif int(amount_time_entries)==2:                                      
        print('2 date entries not supported, check the correct conversion to days')
    elif int(amount_time_entries)==3:
        date = time_conversion.DateTime(int(values[time_entry]),int(values[time_entry+1]),floor(values[time_entry+2]),floor((values[time_entry+2]%1)*24-time_delta),floor((values[time_entry+2]*24)%1*60),floor((values[time_entry+2]*24*60)%1*60*10**6)/10**6)
    elif int(amount_time_entries)==4:
        date = time_conversion.DateTime(int(values[time_entry]),int(values[time_entry+1]),int(values[time_entry+2]),floor(values[time_entry+3]-time_delta),floor((values[time_entry+3]%1)*60),floor((values[time_entry+3]*60)%1*60*10**6)/10**6)
    elif int(amount_time_entries)==5:
        date = time_conversion.DateTime(int(values[time_entry]),int(values[time_entry+1]),int(values[time_entry+2]),int(values[time_entry+3]-time_delta),floor(values[time_entry+4]),floor((values[time_entry+4]%1)*60*10**6)/10**6)
    elif int(amount_time_entries)==6:
        date = time_conversion.DateTime(int(values[time_entry]),int(values[time_entry+1]),int(values[time_entry+2]),int(values[time_entry+3]-time_delta),int(values[time_entry+4]),values[time_entry+5])

    #Conversion from utc to tdb to julian day to seconds
    #epoch = time_conversion.time_from_julian_day(Time(date.julian_day(),scale=time_scale,format='jd').tdb.value)  #Seconds since 2000

    if date.year <1960: #utc only starts at 1960, before sacrificing precision for functionality
        string = str(date.year) + ' ' + str(date.month) + ' ' + str(date.day) + ' ' + str(date.hour) + ' ' + str(date.minute) + ' ' + str(date.seconds) + ' ' + time_scale
        epoch = time_conversion.date_time_from_epoch(spice.convert_date_string_to_ephemeris_time(string))
    else:    
        converter = time_conversion.default_time_scale_converter()
        TimeScales = time_conversion.TimeScales
        epoch = time_conversion.date_time_from_epoch(converter.convert_time(getattr(TimeScales, f"{time_scale}_scale", None),TimeScales.tdb_scale,date.epoch()))
    return epoch

def MJD_to_tudat(values,time_entry,time_scale):
    '''
    Convert a modified julian date into a tudat-compatible time object

    Parameters
    -----------
    values: list[float]
        entries in one line (one observation) from the dataset
    time_entry: int
        the first entry of the time of observation
    time_scale: str
        time standards: utc, ut1, ...


    Examples
    -----------
    .. code-block:: python
        if time_type == 'MJD':
            epoch = MJD_to_tudat(values,time_entry,time_scale) 


    Returns
    -----------
    epoch: obj
        tudat-compatible time object

    '''        
    epoch = time_conversion.datetime_to_tudat(time_conversion.julian_day_to_calendar_date(time_conversion.modified_julian_day_to_julian_day(Time(values[time_entry],scale=time_scale,format='mjd').tdb.value)))  #Seconds since 2000
    return epoch

def minsinceJD_to_tudat(values,time_entry,time_scale,JD):
    '''
    Convert 'minute since julian date' time type into a tudat-compatible time object

    Parameters
    -----------
    values: list[str]
        entries in one line (one observation) from the dataset
    time_entry: int
        the first entry of the time of observation
    time_scale: str
        time standards: utc, ut1, ...
    JD: float
        the reference Julian Date


    Examples
    -----------
    .. code-block:: python
        elif time_type[:10] == 'MINsinceJD':
            JD = 2451545
            epoch = minsinceJD_to_tudat(values,time_entry,time_scale,JD)


    Returns
    -----------
    epoch: obj
        tudat-compatible time object

    '''
    epoch = time_conversion.datetime_to_tudat(time_conversion.julian_day_to_calendar_date(Time(JD + values[time_entry]/(60*24),scale=time_scale,format='jd').tdb.value))
    return epoch


#Environment functions
def add_observatory(bodies,observatory):       #Generate ephemerides for the observatory
    '''
    Add an observatory to the system of bodies    

    Parameters
    -----------
    bodies: obj
        system of bodies relevant for the problem
    observatory: str
        the ID of the observatory


    Examples
    -----------
    .. code-block:: python
        bodies = generate_standard_environment('Jupiter','ECLIPJ2000)
        add_observatory(bodies, observatory)


    Returns
    -----------
    /

    '''
    #Create observation station
    # Define the position of the observatory on Earth
    observatory_longitude, observatory_latitude, observatory_altitude = observatory_info (observatory)

    # Add the ground station to the environment
    environment_setup.add_ground_station(
        bodies.get_body("Earth"),
        'obs' + str(observatory),
        [observatory_altitude, observatory_latitude, observatory_longitude],
        element_conversion.geodetic_position_type)
    return

def get_angle_rel_body(epoch,orientation,observatory_ephemerides, relative_body, standard_orientation):                          #Returns the absolute angular position of the relative body from an observatory
    '''
    Retrieve the absolute angular position of a body, for example the reference body in relative observations or a reference ephemeris from spice    

    Parameters
    -----------
    epoch: obj
        tudat-compatible time object
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...
    observatory_ephemerides: obj
        ephemerides describing the observatory position
    relative_body: str, float or int
        name or ID of the body who's position is requested
    standard_orientation: str
        orientation which will be used for operations, and to which all data will be converted


    Examples
    -----------
    .. code-block:: python
        bodies = generate_standard_environment('Jupiter','ECLIPJ2000)
        add_observatory(bodies, observatory)
        RA_relative, DEC_relative = get_angle_rel_body(epoch,orientation,observatory_ephemerides,relative_body)


    Returns
    -----------
    RA_ECLIPJ2000: float
        absolute right ascension of the requested body
    DEC_ECLIPJ2000: float
        absolute declination of the requested body

    '''    
    #Relative body from spice, from observer, in spherical coordinates
    position = spice.get_body_cartesian_state_at_epoch(relative_body,"Earth",orientation,'LT',epoch.epoch())
    position_correction = observatory_ephemerides.cartesian_state(epoch.epoch())
    rot_mat = spice.compute_rotation_matrix_between_frames(standard_orientation,orientation,epoch.epoch())
    position_corrected = position-np.concatenate([rot_mat@position_correction[:3],rot_mat@position_correction[3:]])
    spherical = element_conversion.cartesian_to_spherical(position_corrected) 

    #Break down into RA and DEC
    RA_relative = spherical[2]
    DEC_relative = spherical[1]

    # RA_relative, DEC_relative = ECLIPJ2000_2_J2000(RA_relative,DEC_relative)

    return RA_relative, DEC_relative          #RADIANS      

def generate_standard_environment(central_body, standard_orientation, extra_bodies=[]):
    '''
    Generate a basic bodies object for a given central body with the major moons, which can be expanded by adding observatories or less significant moons    

    Parameters
    -----------
    central_body: obj
        tudat-compatible time object
    standard_orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...
    extra_bodies: obj
        ephemerides describing the observatory position
    relative_body: str, float or int


    Examples
    -----------
    .. code-block:: python
        bodies = generate_standard_environment('Jupiter','ECLIPJ2000)
        add_observatory(bodies, observatory)
        RA_relative, DEC_relative = get_angle_rel_body(epoch,orientation,observatory_ephemerides,relative_body)


    Returns
    -----------
    RA_ECLIPJ2000: float
        absolute right ascension of the requested body
    DEC_ECLIPJ2000: float
        absolute declination of the requested body

    '''    
    common_bodies = ["Sun", "Earth",central_body]
    print('Loading the standard', central_body, 'environment.')

    if central_body == 'Jupiter':
        standard_moons = ['501', '502', '503', '504', '505','514', '515', '516']

    elif central_body == 'Saturn':
        standard_moons = ['601', '602','603', '604', '605','606','607','608']
    
    elif central_body == 'Uranus':
        standard_moons = ['701', '702','703', '704', '705']

    elif central_body == 'Neptune':
        standard_moons = ['801']

    else:
        standard_moons = []
        print('central body was not correctly interpreted, no moons are loaded. Please check spelling and availability of your central body')

    # Create default body settings for bodies_to_create, with "Earth" as the global frame origin and 
    bodies_to_create = common_bodies+standard_moons+extra_bodies
    global_frame_origin = "Earth"
    global_frame_orientation = standard_orientation
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation)
    
    # Create system of bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)
    return bodies


#Angle process functions
def process_relative(values,RA_entry,DEC_entry,body,relative_body,orientation,observatory,epoch,data_type,bodies,standard_orientation):                   #Optional function for processing the relative data
    '''
    Converts a relative observation into an absolute one    

    Parameters
    -----------
    values: list[str]
        entries in one line (one observation) from the dataset
    RA_entry: int
        the first entry describing the right ascension of the observed body
    DEC_entry: int
        the first entry describing declination of the observed body
    body: str
        the target body either by name or by code
    relative_body: str, float or int
        name or ID of the body who's position is requested
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...
    observatory: str
        the ID of the observatory
    epoch: obj
        tudat-compatible time object
    data_type: str
        describes type of data: absolute, separation arc, relative tangential

        
    Examples
    -----------
    .. code-block:: python
        RA_rad, DEC_rad = process_relative(values,RA_entry,DEC_entry,body,relative_body,orientation,observatory,epoch,data_type)


    Returns
    -----------
    RA_rad: float
        absolute right ascension of the observed body in rad
    DEC_rad: float
        absolute declination of the observed body in rad

    ''' 
    observatory_ephemerides = environment_setup.create_ground_station_ephemeris(
        bodies.get_body("Earth"),
        'obs' + str(observatory),
        )
    RA_relative, DEC_relative = get_angle_rel_body(epoch,orientation,observatory_ephemerides,relative_body, standard_orientation)
    
    if data_type=='REL':
        Delta_DEC = 1/3600* values[DEC_entry]
        Delta_RA = 1/3600* values[RA_entry]/(np.cos(DEC_relative))
    elif data_type=='REL_wo_cos':
        Delta_DEC = 1/3600* values[DEC_entry]
        Delta_RA = 1/3600* values[RA_entry]
    elif data_type=='DIF':
        Delta_RA = 10**(-5) *values[RA_entry]               
        Delta_DEC = 10**(-5) * values[DEC_entry]

    #Correct RA and DEC based on spice of relative body
    RA_rad = (RA_relative + np.deg2rad(Delta_RA))%(2*np.pi)
    DEC_rad = DEC_relative + np.deg2rad(Delta_DEC)

    return RA_rad, DEC_rad                  #note: same orientation as input, has to be rotated as desired

def abs_uniform_deg_output(values,RA_entry,DEC_entry):      #Converts the RA and DEC input into a single value in degrees for ABS inputs
    '''
    Converts a absolue observation into an uniformised output to ease processing    

    Parameters
    -----------
    values: list[str]
        entries in one line (one observation) from the dataset
    RA_entry: int
        the first entry describing the right ascension of the observed body
    DEC_entry: int
        the first entry describing declination of the observed body


    Examples
    -----------
    .. code-block:: python
        elif data_type == 'ABS':
            RA_deg, DEC_deg = abs_uniform_deg_output(values,RA_entry,DEC_entry)


    Returns
    -----------
    RA_rad: float
        absolute right ascension of the observed body in rad
    DEC_rad: float
        absolute declination of the observed body in rad

    '''
    if abs(RA_entry - DEC_entry) != 1:
        
        DEC_sign, RA_sign = 1, 1
        if str(values[DEC_entry]).startswith("-") or values[DEC_entry+1] < 0 or values[DEC_entry+2] < 0:
            DEC_sign = -1
        if str(values[RA_entry]).startswith("-") or values[RA_entry+1] < 0 or values[RA_entry+2] < 0:
            RA_sign = -1
        RA_deg = (15*abs(values[RA_entry])+ 1/4*abs(values[RA_entry+1]) + 1/240*abs(values[RA_entry+2]))*RA_sign                  
        DEC_deg = (abs(values[DEC_entry])+ 1/60*abs(values[DEC_entry+1]) + 1/3600*abs(values[DEC_entry+2]))*DEC_sign
    else:
        RA_deg = values[RA_entry]
        DEC_deg = values[DEC_entry]
    RA_rad = np.deg2rad(RA_deg)
    DEC_rad = np.deg2rad(DEC_deg)
    return RA_rad, DEC_rad


#Meta Functions
def return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body,bodies,standard_orientation): #Meta-function for standardised output
    '''
    Meta function to convert every type of observation into uniformised absolute observations   

    Parameters
    -----------
    values: list[str]
        entries in one line (one observation) from the dataset
    RA_entry: int
        the first entry describing the right ascension of the observed body
    DEC_entry: int
        the first entry describing declination of the observed body
    body: str
        the target body either by name or by code
    relative_body: str, float or int
        name or ID of the body who's position is requested
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...
    observatory: str
        the ID of the observatory
    epoch: obj
        tudat-compatible time object
    data_type: str
        describes type of data: absolute, separation arc, relative tangential
    central_body: str, float or int
        number of planet around sun     

        
    Examples
    -----------
    .. code-block:: python
        for observation in datafile:    
            body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)
            epoch = return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type)
            RA_deg, DEC_deg = return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body)
            RA_EJ2000, DEC_EJ2000 = return_standardised_rotation(RA_deg,DEC_deg,epoch,orientation)


    Returns
    -----------
    RA_rad: float
        absolute right ascension of the observed body in rad
    DEC_rad: float
        absolute declination of the observed body in rad

    ''' 
    if not bodies.does_body_exist(body):
        body_setting = environment_setup.get_default_body_settings([body],"SSB",standard_orientation)
        body_setting.get(body).gravity_field_settings = environment_setup.gravity_field.central(10**(9))    #Dummy value, irrelevant for entire process
        new_body = environment_setup.create_system_of_bodies(body_setting)
        bodies.add_body(body_to_add= new_body.get(body), body_name = body)
        print("Adding", body)
    if not bodies.does_body_exist('obs' + str(observatory)):
        add_observatory(bodies, observatory)
        

    if data_type=='REL' or data_type=='DIF' or data_type == 'REL_wo_cos':
        #Retrieval of relative body
        if relative_body_entry.isalpha():
            relative_body = relative_body_entry
        else:
            relative_body = moon_ID(values[int(relative_body_entry)],central_body)
        
        
        if not bodies.does_body_exist(relative_body):
            body_setting = environment_setup.get_default_body_settings([relative_body],"SSB",standard_orientation)
            new_body = environment_setup.create_system_of_bodies(body_setting)
            body_setting.get(relative_body).gravity_field_settings = environment_setup.gravity_field.central(10**(9))    #Dummy value, irrelevant for entire process
            bodies.add_body(body_to_add= new_body.get(relative_body), body_name = relative_body)
            print("Adding", relative_body)
        RA_rad, DEC_rad = process_relative(values,RA_entry,DEC_entry,body,relative_body,orientation,observatory,epoch,data_type,bodies,standard_orientation)
        
    elif data_type == 'ABS':
        RA_rad, DEC_rad = abs_uniform_deg_output(values,RA_entry,DEC_entry)

    else:
        print('Unknown data type')
    return RA_rad, DEC_rad

def return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type):                    #Meta-function for standardised time
    '''
    Meta function to convert every type of time reporting into uniformised tudat-compatible time object  

    Parameters
    -----------
    values: list[str]
        entries in one line (one observation) from the dataset
    time_entry: int
        the first entry of the time of observation
    time_scale: str
        time standards: utc, ut1, ...
    amount_time_entries: int
        the amount of time entries per observation, necessary in processing. IMPORTANT: assumes year-month-day hour-minute-second framework. Different order would require manual fixes
    time_delta: float
        offset between the time standard and the reference epoch for the observation set  
    time_type: str
        time format: JD, MJD, J2000, DATE, MINsinceJDaddingtheJD
        

    Examples
    -----------
    .. code-block:: python
        for observation in datafile:    
            body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)
            epoch = return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type)
            RA_deg, DEC_deg = return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body)
            RA_EJ2000, DEC_EJ2000 = return_standardised_rotation(RA_deg,DEC_deg,epoch,orientation)


    Returns
    -----------
    epoch: obj
        tudat-compatible time object

    '''  
    if time_type == 'JD':
        epoch = JD_to_tudat(values,time_entry,time_scale) 
    elif time_type == 'DATE': 
        epoch = date_to_tudat(values,time_entry,time_scale,amount_time_entries,time_delta)
    elif time_type == 'MJD':
        epoch = MJD_to_tudat(values,time_entry,time_scale)
    elif time_type[:10] == 'MINsinceJD':
        JD = float(time_type[10:])
        epoch = minsinceJD_to_tudat(values,time_entry,time_scale,JD)
    else:
        print('Unknown Time type')
    return epoch

def return_standardised_rotation_ECLIPJ2000(RA_rad,DEC_rad,epoch,orientation):     #Convert from any orientation to ECLIPJ2000
    '''
    Meta function to convert every type of orientation into a uniform angular position oriented in the ECLIPJ2000 frame  

    Parameters
    -----------
    RA_rad: float
        absolute right ascension of the observed body in rad
    DEC_rad: float
        absolute declination of the observed body in rad
    epoch: obj
        tudat-compatible time object
    orientation: str
        the orientation used in reporting the observation: J2000, B1950, ...

        
    Examples
    -----------
    .. code-block:: python
        for observation in datafile:    
            body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)
            epoch = return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type)
            RA_rad, DEC_rad = return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body)
            RA_EJ2000, DEC_EJ2000 = return_standardised_rotation(RA_deg,DEC_deg,epoch,orientation)


    Returns
    -----------
    RA_ECLIPJ2000: float
        right ascension in the ECLIPJ2000 reference frame in rad
    DEC_ECLIPJ2000: float
        declination in the ECLIPJ2000 reference frame in rad

    '''  
    #any -> J2000/ICRS
    if orientation != 'J2000':
        RA_J2000, DEC_J2000 = Byear_2_J2000(RA_rad,DEC_rad,epoch,orientation)

    else:
        RA_J2000 = RA_rad
        DEC_J2000 = DEC_rad

    #Rotation from J2000 to ECLIPJ2000
    RA_EJ2000, DEC_EJ2000 = J2000_2_ECLIPJ2000(RA_J2000,DEC_J2000)
    return RA_EJ2000, DEC_EJ2000

def process_file(filename,analyse, standard_orientation):
    print('processing', filename)
    '''
    Meta function to fully process a file from reading to saving as seperate csv files per observed body-observatory combination

    Parameters
    -----------
    filename: str
        name of the nsdc datafile to read

        
    Examples
    -----------
    .. code-block:: python
        for filename in nsdc_folder:    
            process_file(filename)


    Returns
    -----------
    -

    '''  
    #Meta_data
    data_type,moon_entry,time_entry,amount_time_entries,time_type,time_scale,time_delta,RA_entry,DEC_entry,number_observatories,observatory,observatory_entrace,relative_body_entry,process_relative,orientation,central_body = read_file_settings(filename)
    with open(filename, 'r') as datafile:
        studyname = filename.split("\\")[-1][:-4]
        # Create lists to store the data
        data, times, RADEC, moons,observatories = [],[],[],[],[]

        #Skip first 16 lines of metadata
        for skips in range(16):
            next(datafile)

        bodies = generate_standard_environment(central_body, standard_orientation)
        # Read each remaining line in the file
        for observation in datafile: 
            body, observatory, values = read_line(observation,number_observatories,observatory,observatory_entrace,moon_entry,central_body)
            epoch = return_standardised_time(values,time_entry,time_scale,amount_time_entries,time_delta,time_type)
            RA_rad, DEC_rad = return_standardised_angles(values,RA_entry,DEC_entry,body,relative_body_entry,orientation,observatory,epoch,data_type,central_body,bodies,standard_orientation)
            RA_EJ2000, DEC_EJ2000 = return_standardised_rotation_ECLIPJ2000(RA_rad,DEC_rad,epoch,orientation)

            ### Saving extracted values
            data.append(values)
            times.append(epoch.epoch())
            RADEC.append([RA_EJ2000,DEC_EJ2000])
            moons.append(body)
            observatories.append(observatory)
    if analyse == True:
        times, RADEC, moons, observatories, diflist, rms, means, stddevs = remove_outliers(times, RADEC, moons, observatories, bodies,standard_orientation, 3, True)
    if (abs(means[0]) < stddevs[0] and abs(means[1]) <stddevs[1]):
        save_to_csv(times, RADEC, moons, observatories, studyname, diflist, 'ObservationsProcessed/ProcessedForThesisWeights/')
    return


