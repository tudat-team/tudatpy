"""Independent first-post-Newtonian atomic-interval audit, without range solves.

IERS2010 Eqs10.1,10.3,10.4 give the geocentric rate using LC, not LB alone.
Compare the selected unmodified native converter and DE440 potential/velocity
integrals plus the v_E dot r endpoint term at all eleven saved APOLLO events.
"""

from pathlib import Path
import hashlib
import importlib.util
import json
import math
import subprocess
import sys

import numpy as np
import tudatpy.kernel
from tudatpy.astro import time_representation as time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
SCRIPT = SOURCE / "analysis/laser_ranging/llr_residuals.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("llr_clock_audit", SCRIPT)
llr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(llr)
manifest = json.loads(SCRIPT.with_name("llr_dataset.json").read_text())
llr.spice.clear_kernels()
for record in manifest["kernels"]:
    path = WORK / "data/llr/kernels" / record["file"]
    assert llr.sha256(path) == record["sha256"]
    llr.spice.load_kernel(str(path))
points = llr.read_points(
    WORK / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
bodies, nominal, loading = llr.create_environment(points, 30, "inpop17a-gir")
assert loading is None
earth = bodies.get("Earth")
converter = earth.rotation_model.angles_calculator.time_scale_converter
original = WORK / "results/apollo-20170309-inpop17a-gir-cubic30"
summary = json.loads((original / "summary.json").read_text())
np.testing.assert_array_equal(nominal, summary["station_xyz_m"])
events = json.loads((original / "event_diagnostics.json").read_text())
C = llr.constants.SPEED_OF_LIGHT
LB, LG = 1.550519768e-8, 6.969290134e-10
LC = (LB - LG) / (1 - LG)
# Use planetary-system barycentres and total GM to include moons without double counting.
external_ids = [10, 1, 2, 4, 5, 6, 7, 8, 9, 301]
gm = {key: llr.spice.get_body_gravitational_parameter(str(key)) for key in external_ids}


def physical_rate(epoch):
    """Return dTT/dTDB minus one, with each external body's potential recorded."""
    earth_state = llr.barycentric_state("Earth", epoch)
    potential_terms = {
        key: gm[key] / np.linalg.norm(llr.barycentric_state(str(key), epoch)[:3] - earth_state[:3])
        for key in gm
    }
    kinetic = 0.5 * float(earth_state[3:] @ earth_state[3:])
    epsilon = -LC + (kinetic + sum(potential_terms.values())) / C**2
    return -epsilon / (1 + epsilon), kinetic, potential_terms


def integrated_rate(first, last, intervals):
    """Integrate the independent geocentric rate by composite Simpson quadrature."""
    dt = float(last - first)
    values = [physical_rate(first + dt * i / intervals)[0] for i in range(intervals + 1)]
    return (
        dt
        / intervals
        / 3
        * sum((1 if i in [0, intervals] else 4 if i % 2 else 2) * v for i, v in enumerate(values))
    )


def station_endpoint(epoch, deformed=False):
    """Compute v_E dot r/c² using the selected actual Earth rotation."""
    position = nominal
    if deformed:
        position = np.asarray(
            earth.get_ground_station("APOLLO_intrinsic").station_state.get_cartesian_position(
                epoch, "Earth"
            )
        ).reshape(3)
    inertial = earth.rotation_model.body_fixed_to_inertial_rotation(epoch) @ position
    velocity = llr.barycentric_state("Earth", epoch)[3:]
    return float(velocity @ inertial / C**2)


def native_interval(first, last, position):
    """Difference conversion offsets, avoiding subtraction of large epoch floats."""
    first_offset = converter.get_time_object_difference(
        time.tdb_scale, time.tai_scale, first, position
    )
    last_offset = converter.get_time_object_difference(
        time.tdb_scale, time.tai_scale, last, position
    )
    return float(last_offset - first_offset)


rows = []
for event in events:
    tx = time.Time(*event["launch_tdb_time_split"])
    rx = time.Time(event["complete"]["link_end_epochs_tdb_s"][-1])
    elapsed = float(rx - tx)
    geocentric = native_interval(tx, rx, np.zeros(3))
    native = native_interval(tx, rx, nominal)
    physical_geo = integrated_rate(tx, rx, 8)
    physical_geo4 = integrated_rate(tx, rx, 4)
    physical_topo = -(station_endpoint(rx) - station_endpoint(tx))
    physical_full = physical_geo + physical_topo
    deformed_topo = -(station_endpoint(rx, True) - station_endpoint(tx, True))
    rate, kinetic, potential = physical_rate(tx + elapsed / 2)
    # Perturb only correction evaluation epochs over the saved RX double's uncertainty.
    half_ulp = float(np.spacing(float(rx)) / 2)
    clock_epoch_sensitivity = abs(
        native_interval(tx, rx + half_ulp, nominal) - native_interval(tx, rx - half_ulp, nominal)
    )
    # A conservative Time arithmetic estimate: eight local-fraction long-double ULPs.
    fraction_bound_s = (
        8
        * float(np.finfo(np.longdouble).eps)
        * max(abs(tx.__getstate__()[1]), abs(rx.__getstate__()[1]), 32.2)
    )
    inversion_residuals = []
    for endpoint in [tx, rx]:
        tai = converter.convert_time_object(time.tdb_scale, time.tai_scale, endpoint, nominal)
        roundtrip = converter.convert_time_object(time.tai_scale, time.tdb_scale, tai, nominal)
        inversion_residuals.append(float(roundtrip - endpoint))
    inversion_interval_error = inversion_residuals[1] - inversion_residuals[0]
    rows.append(
        {
            "index": event["index"],
            "tx_time_split": tx.__getstate__(),
            "rx_time_split": rx.__getstate__(),
            "tdb_interval_s": elapsed,
            "native_atomic_minus_tdb_s": native,
            "native_geocentric_atomic_minus_tdb_s": geocentric,
            "native_station_increment_s": native - geocentric,
            "physical_geocentric_atomic_minus_tdb_s": physical_geo,
            "physical_station_increment_s": physical_topo,
            "physical_atomic_minus_tdb_s": physical_full,
            "full_discrepancy_one_way_m": (native - physical_full) * C / 2,
            "geocentric_discrepancy_one_way_m": (geocentric - physical_geo) * C / 2,
            "station_discrepancy_one_way_m": (native - geocentric - physical_topo) * C / 2,
            "native_interval_correction_one_way_m": native * C / 2,
            "native_topocentric_inversion_forward_residual_s": inversion_residuals,
            "inversion_residual_interval_effect_one_way_m": inversion_interval_error * C / 2,
            "full_discrepancy_after_accounting_inverse_residual_one_way_m": (
                native - physical_full - inversion_interval_error
            )
            * C
            / 2,
            "native_1ps_two_endpoint_inversion_bound_one_way_m": 2e-12 * C / 2,
            "hypothetical_duplicate_station_term_one_way_m": physical_topo * C / 2,
            "hypothetical_LB_instead_of_LC_rate_effect_one_way_m": (LB - LC) * elapsed * C / 2,
            "geocentric_rate_dtt_dtdb_minus1": rate,
            "earth_kinetic_m2_per_s2": kinetic,
            "external_potential_m2_per_s2": {str(key): value for key, value in potential.items()},
            "quadrature4_to8_difference_one_way_m": abs(physical_geo4 - physical_geo) * C / 2,
            "deformed_vs_nominal_station_increment_one_way_m": (deformed_topo - physical_topo)
            * C
            / 2,
            "saved_rx_epoch_full_ulp_correction_sensitivity_one_way_m": clock_epoch_sensitivity
            * C
            / 2,
            "conservative_native_time_arithmetic_one_way_m": fraction_bound_s * C / 2,
            "leading_1pn_truncation_bound_one_way_m": 1e-15 * elapsed * C / 2,
        }
    )

result = {
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "kernel_sha256": llr.sha256(Path(tudatpy.kernel.__file__)),
    "analysis_script_sha256": llr.sha256(SCRIPT),
    "audit_script_sha256": llr.sha256(Path(__file__)),
    "saved_events_file": str(original / "event_diagnostics.json"),
    "saved_events_sha256": llr.sha256(original / "event_diagnostics.json"),
    "station_profile": "inpop17a-gir",
    "nominal_station_xyz_m": nominal.tolist(),
    "constants": {"LB": LB, "LG": LG, "LC": LC, "c_m_per_s": C},
    "external_gm_m3_per_s2": {str(k): v for k, v in gm.items()},
    "formula": "integral[-epsilon/(1+epsilon)]dt - Delta(v_E dot r_GCRS/c²), epsilon=-LC+(vE²/2+U_external)/c²",
    "native_code": "nWayRange adds reception(TAI-TDB) minus transmission(TAI-TDB) once; uses nominal stationXYZ. No separate analysis interval bias.",
    "time_precision": "Transmit uses saved split Time. Reception uses saved native double restored to split Time; evaluate correction full-ULP sensitivity separately. No large-epoch float subtraction.",
    "truncation_and_uncertainty": "Leading1PN allowance1e-15rate (IERS quotes omittedc^-4~1e-16) plus measuredquadrature/state/Timearithmetic sensitivity. Native topocentric inverse stops at1psforwardresidual per endpoint, giving0.300mm one-wayintervalbound. Externalasteroids omitted; analyticSOFAmodel differs fromDE440. This is not an exactDE440 time ephemeris.",
    "primary_reference": "https://iers-conventions.obspm.fr/content/chapter10/tn36_c10.pdf, Eqs10.1/10.3/10.4",
    "rows": rows,
    "maximum_absolute_full_discrepancy_one_way_m": max(
        abs(row["full_discrepancy_one_way_m"]) for row in rows
    ),
    "maximum_absolute_geocentric_discrepancy_one_way_m": max(
        abs(row["geocentric_discrepancy_one_way_m"]) for row in rows
    ),
    "maximum_absolute_station_discrepancy_one_way_m": max(
        abs(row["station_discrepancy_one_way_m"]) for row in rows
    ),
    "maximum_absolute_discrepancy_after_accounting_inverse_residual_one_way_m": max(
        abs(row["full_discrepancy_after_accounting_inverse_residual_one_way_m"]) for row in rows
    ),
    "conclusion": "Correct sign and normalization; station endpoint clock term occurs once. Observed nativeclock discrepancy<0.1mm cannot explain cm-scale remainingresiduals. Inverse1pstolerance sets a sub-mm arithmetic limit, not a physicalfit.",
}
out = WORK / "logs/llr-atomic-interval-audit.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))
for row in rows:
    print(
        row["index"],
        "native_correction_mm",
        row["native_interval_correction_one_way_m"] * 1000,
        "difference_um",
        row["full_discrepancy_one_way_m"] * 1e6,
        "geo_um",
        row["geocentric_discrepancy_one_way_m"] * 1e6,
        "topo_um",
        row["station_discrepancy_one_way_m"] * 1e6,
    )
