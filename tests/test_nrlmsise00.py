# Load standard modules
import os
import pickle
import numpy as np

# Load Tudatpy modules
from tudatpy.kernel import numerical_simulation
from tudatpy.data import read_solar_activity_data, get_space_weather_path
from tudatpy.kernel.astro.time_conversion import epoch_from_date_time_iso_string


def test_nrlmsise00():

    # Load the validation data obtained with the pymsis library. 
    # The data contains the generated input data (date [iso string], longitude [deg], latitude [deg], altitude [meters]), 
    # pymsis resulting densities [kg/m^3], and the corresponding space weather parameters (F10.7, F10.7a, and daily Ap) used 
    # internally by pymsis for future references.
    with open(os.path.dirname(__file__) + '/nrlmsise00_validation_data.pkl', "rb") as file:
        validation_data = pickle.load(file)

    # Initialise Tudat NRLMSISE00 model
    solar_weather_data = read_solar_activity_data(get_space_weather_path() + '/sw19571001.txt')
    NRLMSISE00 = numerical_simulation.environment_setup.atmosphere.NRLMSISE00Atmosphere(solar_weather_data,True,False,True)

    # Extract input data and expected densities from the list of dictionaries
    altitudes = np.array([data["altitude"] for data in validation_data])
    longitudes = np.array([data["longitude"] for data in validation_data])
    latitudes = np.array([data["latitude"] for data in validation_data])
    epochs = np.array([epoch_from_date_time_iso_string(data["date"]) for data in validation_data])
    expected_densities = np.array([data["density"] for data in validation_data])

    # Compute densities using Tudat NRLMSISE00
    computed_densities = np.array([
        NRLMSISE00.get_density(alt, lon, lat, epoch)
        for alt, lon, lat, epoch in zip(altitudes, longitudes, latitudes, epochs)
    ])
        
    # Compute the relative differences between computed and expected densities
    relative_differences = np.abs(computed_densities - expected_densities) / expected_densities

    
    print(relative_differences.max())
    # Check if all relative differences are within 1e-2
    assert np.all(relative_differences < 1e-2)
