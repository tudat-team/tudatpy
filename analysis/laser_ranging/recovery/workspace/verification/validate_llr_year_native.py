from pathlib import Path
import sys, json, collections
from datetime import date

sys.path.insert(0, "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
import llr_residuals as model
from tudatpy.data_input.tracking_data import slr

r = Path("/home/dominic/Tudat/laser-ranging-work/research/llr-year-2017")
a = json.loads((r / "apollo_subset_sources.json").read_text())
p = model.read_points(r / a["file"], a["sha256"], (date(2017, 1, 1), date(2018, 1, 1)), 311)
required = dict(
    cdp_pad_id=7845,
    system_number=78,
    occupancy=1,
    station_epoch_time_scale=4,
    range_type=2,
    data_quality_alert=0,
    tropospheric_correction_applied=0,
    center_of_mass_correction_applied=0,
    receive_amplitude_correction_applied=0,
    station_system_delay_applied=1,
    spacecraft_system_delay_applied=0,
)
counts = collections.Counter()
flags = collections.Counter()
weather = collections.Counter()
errors = []
for path in sorted((r / "crd").glob("*.npt")):
    for session in slr.read_slr_data([str(path)]):
        c = session.configuration
        counts[str(c.cdp_pad_id)] += len(session.data.normal_point_data)
        if c.cdp_pad_id != 7845:
            continue
        assert all(getattr(c, k) == v for k, v in required.items()), (path.name, c.cdp_pad_id)
        assert len(slr.extract_normal_point_measurements(session)) == len(
            session.data.normal_point_data
        )
        for m in session.data.meteorological_data:
            assert (
                m.origin_of_values == 0
                and 0 < m.pressure < 1100
                and 200 < m.temperature < 340
                and 0 <= m.humidity <= 100
            )
            weather["measured_records"] += 1
        for n in session.data.normal_point_data:
            assert n.epoch_event == 2
            w = c.system_configurations[n.system_configuration_id].transmit_wavelength_nm
            assert w in [532.1, 1064.2]
            flags[str(w)] += 1
out = dict(
    native_kernel_sha256=model.sha256(
        Path("/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so")
    ),
    station_np_counts=dict(counts),
    grasse_wavelength_counts=dict(flags),
    weather=dict(weather),
    apollo=dict(
        count=len(p),
        first_utc=min(x["utc_label"] for x in p),
        last_utc=max(x["utc_label"] for x in p),
        nights=len({x["utc_label"][:10] for x in p}),
        reflector_counts=dict(collections.Counter(str(x["reflector"]) for x in p)),
        source_metadata=dict(
            station_ids={"70610": len(p)},
            wavelength_angstroms=dict(
                collections.Counter(str(round(x["wavelength_m"] * 1e10)) for x in p)
            ),
            quality_grades=dict(collections.Counter(x["quality"] for x in p)),
            pressure_hPa_range=[
                min(x["pressure_pa"] for x in p) / 100,
                max(x["pressure_pa"] for x in p) / 100,
            ],
            temperature_C_range=[
                min(x["temperature_k"] for x in p) - 273.15,
                max(x["temperature_k"] for x in p) - 273.15,
            ],
            relative_humidity_percent_range=[
                min(x["humidity"] for x in p) * 100,
                max(x["humidity"] for x in p) * 100,
            ],
            source_rtt_sigma_seconds_range=[
                min(x["sigma_rtt_s"] for x in p),
                max(x["sigma_rtt_s"] for x in p),
            ],
        ),
    ),
)
(r / "native-reader-validation.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
