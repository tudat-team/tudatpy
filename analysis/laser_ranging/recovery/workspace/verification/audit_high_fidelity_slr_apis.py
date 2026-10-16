"""Read-only API/resource inventory for the high-fidelity SLR configuration."""

from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import sys

import tudatpy
import tudatpy.kernel
from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import atmospheric_loading, ocean_loading, station_loading
from tudatpy.dynamics import environment
from tudatpy.dynamics.environment_setup import rotation_model, shape_deformation
from tudatpy.estimation.observable_models_setup import model_settings

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
resources = Path(resource_paths.get_resource_path())
refs = Path("/home/dominic/Tudat/earth-station-model-work/references")
kernel = Path(tudatpy.kernel.__file__)
digest = hashlib.sha256(kernel.read_bytes()).hexdigest()
assert digest == "10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5"
eop = rotation_model.EarthOrientationEvaluationSettings()
atmosphere = shape_deformation.AtmosphericTideEvaluationSettings()
station = shape_deformation.StationLoadingEvaluationSettings()
apis = {
    "rotation_model": [
        rotation_model,
        "gcrs_to_itrs",
        "EarthOrientationEvaluationSettings",
        "EarthOrientationInterpolation",
        "OceanTideEopModel",
    ],
    "shape_deformation": [
        shape_deformation,
        "iers_2010_solid_body_tidal",
        "pole_tidal",
        "ocean_pole_tidal",
        "ocean_pole_tidal_from_data",
        "atmospheric_tide_from_grid_file",
        "atmospheric_tide_from_grid",
        "station_loading",
    ],
    "atmospheric_loading": [
        atmospheric_loading,
        "read_atmospheric_tide_grid",
        "AtmosphericTideMetadata",
        "AtmosphericTideSourceProfile",
    ],
    "ocean_loading": [ocean_loading, "read_ocean_pole_tide_grid"],
    "station_loading": [
        station_loading,
        "read_ephedisp_station_loading_data",
        "imls_merra2_atmospheric_loading_metadata",
    ],
    "model_settings": [model_settings, "two_way_range", "n_way_range"],
}
api_inventory = {}
for group, items in apis.items():
    module, *names = items
    api_inventory[group] = {
        name: {"available": hasattr(module, name), "doc": getattr(module, name).__doc__}
        for name in names
    }
files = [
    resources / "earth_orientation" / name
    for name in [
        "README_desaiSibois2016.txt",
        "desaiSibois2016PolarMotionOceanTidesAmplitudes.txt",
        "desaiSibois2016PolarMotionOceanTidesFundamentalArgumentMultipliers.txt",
        "desaiSibois2016Ut1OceanTidesAmplitudes.txt",
        "desaiSibois2016Ut1OceanTidesFundamentalArgumentMultipliers.txt",
        "eopc04_20_IAU2000.62-now.txt",
    ]
]
files += [
    resources / "earth_deformation" / name
    for name in [
        "diurnalDisplacementFrequencyDependence2.txt",
        "longPeriodDisplacementFrequencyDependence.txt",
    ]
]
files += [
    refs / "ocean-pole-tide" / "opoleloadcoefcmcor.txt",
    refs / "ocean-pole-tide" / "opoleloadcoefcmcor.txt.gz",
    refs / "atmospheric-tides" / "s1_s2_def_ce.dat",
    refs / "atmospheric-tides" / "S1_S2_CM_NOIB_VAPL.dat",
]
files += sorted((WORK / "research/slr-loading-202404").glob("*.eph"))
inventory = []
for path in files:
    item = {"path": str(path), "available": path.is_file()}
    if path.is_file():
        raw = path.read_bytes()
        item.update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        if "desaiSibois2016" in path.name:
            item["nonempty_rows"] = len(
                [line for line in raw.decode().splitlines() if line.strip()]
            )
    inventory.append(item)
rotation_methods = {
    name: getattr(environment.RotationalEphemeris, name).__doc__
    for name in [
        "body_fixed_to_inertial_rotation",
        "inertial_to_body_fixed_rotation",
        "time_derivative_body_fixed_to_inertial_rotation",
        "time_derivative_inertial_to_body_fixed_rotation",
    ]
}
rotation_methods["get_gcrs_to_itrs_rotation_angles"] = (
    environment.EarthOrientationAnglesCalculator.get_gcrs_to_itrs_rotation_angles.__doc__
)
result = {
    "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "interpreter": sys.executable,
    "package_origin": tudatpy.__file__,
    "kernel_origin": str(kernel),
    "kernel_sha256": digest,
    "scope": "API/docstring/default-settings and resource-byte inspection only; no numerical model tests, builds, or range solves.",
    "api_inventory": api_inventory,
    "resources": inventory,
    "native_defaults": {
        "eop_interpolation": str(eop.interpolation),
        "eop_maximum_gap_seconds": eop.maximum_gap_seconds,
        "ocean_eop_model": "iers_2010 (gcrs_to_itrs default argument)",
        "atmospheric_lookup": str(atmosphere.lookup_latitude_convention),
        "atmospheric_displacement_basis": str(atmosphere.displacement_basis),
        "atmospheric_convert_to_centre_of_mass": atmosphere.convert_to_centre_of_mass,
        "station_loading_maximum_gap_seconds": station.maximum_gap_seconds,
        "station_loading_station_bindings": station.station_bindings,
        "station_loading_spatial_selection": str(station.spatial_selection),
        "automatic_deformation_addons": "None: each desired shape deformation must be attached explicitly.",
    },
    "rotation_binding_docs": rotation_methods,
    "compatibility_findings": {
        "rotation_and_deformations": "Selected rotation owns one DS/four-point evaluator+converter, shared by solid Earth/pole/ocean-pole/atmospheric/station-loading factory paths.",
        "observable_tai_interval": "NOT SHARED. one/n-way factory calls no-arg setTimeScaleConverter(), which creates/reuses independent default converter with IERS2010 subdaily UT1 and linear daily EOP. No Python injection API found.",
        "precise_time_rotation": "Native GCRS rotation and derivative have Time overrides; Python full matrices/derivatives bind only double. Time-valued EarthOrientationAnglesCalculator input is exposed, but not full Time matrix/derivative overloads.",
        "eop_selection": "DS159 terms selectable, replaces conventional ocean-tide terms and retains libration; four-point daily interpolation selectable independently. Custom short-period files rejected when auto-selecting DS.",
        "deformation_data": "Official ocean-pole and two atmospheric grids cached externally; not auto-selected. Ray/Ponte grid identified by header/profile supports geocentric default and optional CE→CM translation; VAPL needs explicit conventions/origin, not filename inference.",
        "station_loading": "Read-only APIs present; native default explicit station bindings and86400s maxgap, source Cartesian frame must exactly match Earth target frame. Visible20240401 cached files alone do not establish full-week coverage.",
    },
    "source_evidence": {
        "src/tudat/simulation/environment_setup/createRotationModelEarthOrientation.cpp": "DS replacement and selected evaluator/UT1 converter wiring",
        "src/tudat/simulation/environment_setup/createBodyDeformationModel.cpp": "Earth deformation models use rotation-owned converter",
        "include/tudat/astro/observation_models/observationModel.h:351": "setTimeScaleConverter() calls createDefaultTimeConverter(nullptr,true)",
        "include/tudat/simulation/estimation_setup/createObservationModelFactory.h:1151": "n-way model factory invokes independent default converter",
        "src/tudat/astro/earth_orientation/terrestrialTimeScaleConverter.cpp:348": "default converter builds default linear evaluator and getDefaultUT1CorrectionCalculator()",
        "src/tudatpy/dynamics/environment/expose_environment.cpp:2815": "Python full rotation matrix binds double method",
        "include/tudat/astro/ephemerides/itrsToGcrsRotationModel.h:91": "native Time rotation override exists",
    },
}
(WORK / "logs/high-fidelity-slr-compatibility-audit.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(
    json.dumps(
        {
            "kernel_sha256": digest,
            "api_groups": {
                key: all(item["available"] for item in values.values())
                for key, values in api_inventory.items()
            },
            "native_defaults": result["native_defaults"],
            "compatibility_findings": result["compatibility_findings"],
            "resource_count": len(inventory),
        },
        indent=2,
    )
)
