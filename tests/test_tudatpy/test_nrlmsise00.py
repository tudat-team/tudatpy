# Load standard modules
import os
import pickle
import numpy as np

# Load Tudatpy modules
from tudatpy.dynamics import environment_setup
from tudatpy.data import read_solar_activity_data, get_space_weather_path, get_atmosphere_tables_path
from tudatpy.astro.time_representation import iso_string_to_epoch
from tudatpy.astro.element_conversion import convert_geographic_to_geodetic_latitude


def test_nrlmsise00():

    # Load the validation data obtained with the pymsis library. 
    # The data contains the generated input data (date [iso string], longitude [deg], latitude [deg], altitude [meters]), 
    # pymsis resulting densities [kg/m^3], and the corresponding space weather parameters (F10.7, F10.7a, and daily Ap) used 
    # internally by pymsis for future references.
    with open( get_atmosphere_tables_path() + '/nrlmsise00_validation_data.pkl', "rb") as file:
        validation_data = pickle.load(file)

    # Initialise Tudat NRLMSISE00 model
    solar_weather_data = read_solar_activity_data(get_space_weather_path() + '/sw19571001.txt')
    NRLMSISE00 = environment_setup.atmosphere.NRLMSISE00Atmosphere(solar_weather_data,True,False,True)

    # Extract input data and expected densities from the list of dictionaries
    altitudes = np.array([data["altitude"] for data in validation_data])
    longitudes = np.array([data["longitude"] for data in validation_data])
    latitudes = np.array([data["latitude"] for data in validation_data])
    epochs = np.array([iso_string_to_epoch(data["date"]) for data in validation_data])
    expected_densities = np.array([data["density"] for data in validation_data])

    # Convert geographic latitudes to geodetic latitudes as the validation data was obtained by passing geodetic latitudes to pymsis
    geodetic_latitudes = np.zeros(len(latitudes))
    for i in range(len(latitudes)):
        geodetic_latitudes[i] = convert_geographic_to_geodetic_latitude(latitudes[i], 6378137.0, 1.0 / 298.257223563, altitudes[i], 1e-15, 3)

    # Compute densities using Tudat NRLMSISE00
    computed_densities = np.array([
        NRLMSISE00.get_density(alt, lon, lat, epoch)
        for alt, lon, lat, epoch in zip(altitudes, longitudes, geodetic_latitudes, epochs)
    ])
        
    # Compute the relative differences between computed and expected densities
    relative_differences = np.abs(computed_densities - expected_densities) / expected_densities

    assert np.all(relative_differences < 5e-6)
