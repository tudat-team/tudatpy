import numpy as np

from tudatpy.astro import element_conversion
from tudatpy.dynamics import environment_setup

JPL_TO_MPC_RADAR_STATIONS = {
    "-1": "251",
    "-14": "253",
}

RADAR_STATION_GEODETIC = {
    "JPL:-1": np.array([453.34, np.deg2rad(18.3442199), np.deg2rad(293.2473068)]),
    "JPL:-14": np.array([1001.39, np.deg2rad(35.4259009), np.deg2rad(243.1104618)]),
    "MPC:251": np.array([453.34, np.deg2rad(18.3442199), np.deg2rad(293.2473068)]),
    "MPC:253": np.array([1001.39, np.deg2rad(35.4259009), np.deg2rad(243.1104618)]),
    "251": np.array([453.34, np.deg2rad(18.3442199), np.deg2rad(293.2473068)]),
    "253": np.array([1001.39, np.deg2rad(35.4259009), np.deg2rad(243.1104618)]),
}


def normalize_radar_station_id(source: str, raw_station_id: str | int) -> str:
    source = str(source).strip().upper()
    station_id = str(raw_station_id).strip()
    if ":" in station_id:
        return station_id
    if source == "MPC":
        return station_id.zfill(3)
    return f"{source}:{station_id}"


def get_radar_station_geodetic_positions() -> dict[str, np.ndarray]:
    return {station: position.copy() for station, position in RADAR_STATION_GEODETIC.items()}


def add_radar_ground_stations(
    bodies,
    station_ids,
    station_body: str = "Earth",
    station_positions: dict[str, np.ndarray] | None = None,
) -> None:
    positions = get_radar_station_geodetic_positions()
    if station_positions is not None:
        positions.update(station_positions)

    for station_id in sorted({str(station) for station in station_ids}):
        if station_id in bodies.get(station_body).ground_station_list:
            continue
        if station_id not in positions:
            raise KeyError(
                f"No geodetic radar station position is available for '{station_id}'. "
                "Pass station_positions or add the station to the environment manually."
            )

        station_settings = environment_setup.ground_station.basic_station(
            station_id,
            positions[station_id],
            element_conversion.geodetic_position_type,
        )
        environment_setup.add_ground_station(
            bodies.get_body(station_body),
            station_settings,
        )
