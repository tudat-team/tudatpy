# Load standard modules
import pickle
import numpy as np
import pytest

# Load Tudatpy modules
from tudatpy.dynamics import environment, environment_setup
from tudatpy.data import (
    get_space_weather_path,
    get_atmosphere_tables_path,
)
from tudatpy.astro.time_representation import iso_string_to_epoch
from tudatpy.astro.element_conversion import convert_geographic_to_geodetic_latitude


@pytest.fixture
def nrlmsise00_bodies():
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Earth")
    body_settings.get("Earth").atmosphere_settings = environment_setup.atmosphere.nrlmsise00(
        get_space_weather_path() + "/sw19571001.txt",
        use_storm_conditions=False,
        use_anomalous_oxygen=True,
    )
    return environment_setup.create_system_of_bodies(body_settings)


def test_nrlmsise00_model_exposure(nrlmsise00_bodies):
    model = nrlmsise00_bodies.get("Earth").atmosphere_model
    assert isinstance(model, environment.NRLMSISE00Atmosphere)
    assert isinstance(model, environment.AtmosphereModel)
    assert not hasattr(environment_setup.atmosphere, "NRLMSISE00Atmosphere")

    with pytest.raises(TypeError):
        environment.NRLMSISE00Atmosphere()
    with pytest.raises(TypeError):
        environment.NRLMSISE00Atmosphere({}, True, False, True)

    # The derived model can be passed back through an AtmosphereModel property.
    nrlmsise00_bodies.create_empty_body("Other")
    nrlmsise00_bodies.get("Other").atmosphere_model = model
    assert nrlmsise00_bodies.get("Other").atmosphere_model is model

    # Inherited queries and the model's existing time/latitude controls remain available.
    model.set_use_utc(True)
    assert model.get_use_utc()
    model.set_use_geodetic_latitude(True)
    assert model.get_use_geodetic_latitude()
    epoch = iso_string_to_epoch("2020-01-01T12:00:00")
    assert model.get_density(400e3, 0.0, 0.0, epoch) > 0.0
    assert model.get_temperature(400e3, 0.0, 0.0, epoch) > 0.0
    assert model.get_pressure(400e3, 0.0, 0.0, epoch) > 0.0
    assert model.get_speed_of_sound(400e3, 0.0, 0.0, epoch) > 0.0


def test_nrlmsise00(nrlmsise00_bodies):

    # Load the validation data obtained with the pymsis library.
    # The data contains the generated input data (date [iso string], longitude [deg], latitude [deg], altitude [meters]),
    # pymsis resulting densities [kg/m^3], and the corresponding space weather parameters (F10.7, F10.7a, and daily Ap) used
    # internally by pymsis for future references.
    with open(get_atmosphere_tables_path() + "/nrlmsise00_validation_data.pkl", "rb") as file:
        validation_data = pickle.load(file)

    # Retrieve the NRLMSISE-00 model created by the environment factory.
    NRLMSISE00 = nrlmsise00_bodies.get("Earth").atmosphere_model

    # Extract input data and expected densities from the list of dictionaries
    altitudes = np.array([data["altitude"] for data in validation_data])
    longitudes = np.array([data["longitude"] for data in validation_data])
    latitudes = np.array([data["latitude"] for data in validation_data])
    epochs = np.array([iso_string_to_epoch(data["date"]) for data in validation_data])
    expected_densities = np.array([data["density"] for data in validation_data])

    # Convert geographic latitudes to geodetic latitudes as the validation data was obtained by passing geodetic latitudes to pymsis
    geodetic_latitudes = np.zeros(len(latitudes))
    for i in range(len(latitudes)):
        geodetic_latitudes[i] = convert_geographic_to_geodetic_latitude(
            latitudes[i], 6378137.0, 1.0 / 298.257223563, altitudes[i], 1e-15, 3
        )

    # Compute densities using Tudat NRLMSISE00
    computed_densities = np.array(
        [
            NRLMSISE00.get_density(alt, lon, lat, epoch)
            for alt, lon, lat, epoch in zip(altitudes, longitudes, geodetic_latitudes, epochs)
        ]
    )

    # Compute the relative differences between computed and expected densities
    relative_differences = np.abs(computed_densities - expected_densities) / expected_densities

    assert np.all(relative_differences < 5e-6)
