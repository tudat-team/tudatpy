"""Reconstruct actual fixed-event optical delays and bound convention sensitivities."""

import hashlib
import json
import math
from pathlib import Path

import numpy as np

import llr_residuals as llr
import matera_case
from tudatpy import constants, kernel
from tudatpy.astro import element_conversion


def optical(latitude, height, temperature, pressure, humidity, wavelength, elevation):
    """Evaluate IERS2010 Eqs.9.3-9.10 with the recorded station units."""

    sigma2 = (1e-6 / wavelength) ** 2
    fh = (
        0.01
        * 0.99995995
        * (
            19990.975 * (238.0185 + sigma2) / (238.0185 - sigma2) ** 2
            + 579.55174 * (57.362 + sigma2) / (57.362 - sigma2) ** 2
        )
    )
    fw = 0.003101 * (
        295.235 + 3 * 2.6422 * sigma2 - 5 * 0.032380 * sigma2**2 + 7 * 0.004028 * sigma2**3
    )
    gravity = 1 - 0.00266 * math.cos(2 * latitude) - 0.00000028 * height
    vapor = humidity * 611.2 * math.exp(17.62 * (temperature - 273.15) / (temperature - 30.03))
    dry = 0.002416579 * fh * pressure / 100 / gravity
    wet = 0.0001 * (5.316 * fw - 3.759 * fh) * vapor / 100 / gravity
    tc = temperature - 273.15
    coefficients = [
        0.00121008 + 0.0000017295 * tc + 0.00003191 * math.cos(latitude) - 0.000000018478 * height,
        0.00304965 + 0.000002346 * tc - 0.0001035 * math.cos(latitude) - 0.00000001856 * height,
        0.068777 + 0.00001972 * tc - 0.003458 * math.cos(latitude) + 0.000000106 * height,
    ]
    a, b, c = coefficients
    se = math.sin(elevation)
    mapping = (1 + a / (1 + b / (1 + c))) / (se + a / (se + b / (se + c)))
    return np.array([dry * mapping, wet * mapping])


repo = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr/analysis/laser_ranging")
data = Path("/home/dominic/Tudat/laser-ranging-work/data/llr")
results = Path("/home/dominic/Tudat/laser-ranging-work/results")
apollo = json.loads((repo / "llr_dataset.json").read_text())
points = llr.read_points(data / apollo["normal_points"]["file"], apollo["normal_points"]["sha256"])
matera = matera_case.read_points(data / "apollo15_20220821.npt")
configurations = [
    ("apollo_survey", "apollo-20170309-cubic30", points),
    ("apollo_inpop17a", "apollo-20170309-inpop17a-gir-cubic30", points),
    ("matera", "matera-20220821-de440-cubic30-committed", matera),
]
output = {
    "kernel_path": kernel.__file__,
    "kernel_sha256": llr.sha256(Path(kernel.__file__)),
    "cases": {},
}
for case, directory, selected in configurations:
    summary = json.loads((results / directory / "summary.json").read_text())
    events = json.loads((results / directory / "event_diagnostics.json").read_text())
    xyz = np.asarray(summary["station_xyz_m"])
    height, latitude, longitude = element_conversion.convert_cartesian_to_geodetic_coordinates(
        xyz, 6378137.0, 1 / 298.257223563, 1e-8
    )
    radial_latitude = math.atan2(xyz[2], math.hypot(xyz[0], xyz[1]))
    times = [float(p["utc"]) for p in selected]
    rows = []
    for point, event in zip(selected, events):
        legs = []
        for leg, elevation_key in [
            ("uplink", "transmit_elevation_rad"),
            ("downlink", "receive_elevation_rad"),
        ]:
            complete = event["complete"]
            elapsed = (
                0
                if leg == "uplink"
                else float(
                    complete["link_end_epochs_tdb_s"][-1] - complete["link_end_epochs_tdb_s"][0]
                )
            )
            utc = float(point["utc"]) + elapsed
            temperature, pressure, humidity = [
                np.interp(utc, times, [p[key] for p in selected])
                for key in ["temperature_k", "pressure_pa", "humidity"]
            ]
            elevation = complete[elevation_key][0]
            arguments = [
                latitude,
                height,
                temperature,
                pressure,
                humidity,
                point["wavelength_m"],
                elevation,
            ]
            prediction = optical(*arguments)
            native = complete[leg + "_corrections_s"][1] * constants.SPEED_OF_LIGHT
            low_height = arguments.copy()
            low_height[1] -= 100.0
            radial = arguments.copy()
            radial[0] = radial_latitude
            legs.append(
                {
                    "leg": leg,
                    "elevation_deg": math.degrees(elevation),
                    "pressure_hpa": pressure / 100,
                    "temperature_k": temperature,
                    "humidity": humidity,
                    "independent_delay_m": float(sum(prediction)),
                    "native_delay_m": native,
                    "independent_minus_native_m": float(sum(prediction) - native),
                    "dry_delay_m": float(prediction[0]),
                    "wet_delay_m": float(prediction[1]),
                    "height_minus_100m_delay_change_m": float(
                        sum(optical(*low_height)) - sum(prediction)
                    ),
                    "geocentric_latitude_delay_change_m": float(
                        sum(optical(*radial)) - sum(prediction)
                    ),
                    "pressure_1hpa_delay_change_m": float(prediction[0] / (pressure / 100)),
                    "pressure_sensor_height_2p4m_sensitivity_m": float(
                        prediction[0] * 9.81 * 2.4 / (287.05 * temperature)
                    ),
                }
            )
        rows.append({"index": point["index"], "legs": legs})
    output["cases"][case] = {
        "source_summary_sha256": llr.sha256(results / directory / "summary.json"),
        "source_events_sha256": llr.sha256(results / directory / "event_diagnostics.json"),
        "height_m": height,
        "geodetic_latitude_deg": math.degrees(latitude),
        "geocentric_latitude_deg": math.degrees(radial_latitude),
        "rows": rows,
    }
output["probe_sha256"] = llr.sha256(Path(__file__))
print(json.dumps(output, indent=2))
