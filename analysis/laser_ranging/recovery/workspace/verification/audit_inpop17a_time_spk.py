"""Read three TT−TDB Chebyshev records; no native observation or clock changes."""

from pathlib import Path
import hashlib
import json
import numpy as np
import erfa
import spiceypy as spice

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
path = WORK / "research/llr-ephemeris-controls/inpop17a-spice/inpop17a_TDB_m100_p100_spice_time.bsp"
launch_file = WORK / "logs/llr-absolute-launch-clock-audit.json"
launch = json.loads(launch_file.read_text())
handle = spice.dafopr(str(path))
spice.dafbfs(handle)
assert spice.daffna()
dc, ic = spice.dafus(spice.dafgs(), 2, 6)
name = spice.dafgn()
assert list(ic[:4]) == [1000000001, 1000000000, 1, 2]
assert name == "TT-TDB"
first, last = map(int, ic[-2:])
init, interval, size, number = spice.dafgda(handle, last - 3, last)
size, number = int(size), int(number)
assert size == 38 and (size - 2) % 3 == 0
count = (size - 2) // 3
rows = []
spice.furnsh(str(path))
for point in launch["rows"]:
    hours, seconds = point["launch_tdb_time_split"]
    epoch = np.longdouble(hours) * 3600 + np.longdouble(seconds)
    index = int(np.floor((epoch - init) / interval))
    record = spice.dafgda(handle, first + index * size, first + (index + 1) * size - 1)
    midpoint, radius = record[:2]
    coefficients = record[2:].reshape(3, count)
    assert np.all(coefficients[1:] == 0)
    argument = (epoch - np.longdouble(midpoint)) / np.longdouble(radius)
    value = float(np.polynomial.chebyshev.chebval(argument, coefficients[0].astype(np.longdouble)))
    derivative = float(
        np.polynomial.chebyshev.chebval(
            argument, np.polynomial.chebyshev.chebder(coefficients[0].astype(np.longdouble))
        )
        / radius
    )
    reference = -float(erfa.dtdb(2451545.0, float(epoch / 86400), 0.0, 0.0, 0.0, 0.0))
    state, _ = spice.spkgeo(1000000001, float(epoch), "J2000", 1000000000)
    assert abs(float(state[0]) - value) < 1e-15
    rows.append(
        dict(
            station=point["station"],
            utc=point["utc"],
            launch_tdb_time_split=[hours, seconds],
            record_index_zero_based=index,
            chebyshev_argument=float(argument),
            inpop_tt_minus_tdb_s=value,
            inpop_dtt_minus_tdb_dtdb=derivative,
            sofa_tt_minus_tdb_s=reference,
            inpop_minus_sofa_s=value - reference,
            static_manual_minus_spice_s=value - float(state[0]),
            absolute_epoch_sensitivity_at_100m_per_s_m=100 * abs(value - reference),
        )
    )
spice.dafcls(handle)
spice.kclear()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
result = dict(
    scope="Static provider type2 coefficient evaluation at three saved launch epochs; no native module import, observation solve, model edit, or clock replacement.",
    provider_file=str(path),
    provider_sha256=sha(path),
    provider_bytes=path.stat().st_size,
    provider_release="INPOP17a TDB m100_p100 SPICE archive, official August2017 release",
    provider_readme="research/llr-ephemeris-controls/INPOP17a-README.txt: TT−TDB in seconds, TDB product argument",
    segment_name=name,
    segment_integer_descriptor=ic.tolist(),
    segment_tdb_coverage_s=dc.tolist(),
    record_layout=dict(
        init_tdb_s=float(init),
        interval_s=float(interval),
        record_words=size,
        records=number,
        coefficients_per_component=count,
    ),
    interpretation="Pseudo-position firstcomponent is TT−TDB in seconds; second/third components identicallyzero in inspectedrecords. No kilometre-to-second factor. Direct manualChebyshev agrees with unmodifiedSPICEpolynomial evaluation.",
    launch_source_file=str(launch_file),
    launch_source_sha256=sha(launch_file),
    erfa_version=erfa.__version__,
    spiceypy_version=spice.__version__,
    audit_source_sha256=sha(__file__),
    sign_and_clock_scope="Compare geocentric TT−TDB only. Native topocentric stationterm must remain separate. Loadingtime.bsp doesnotautomaticallyreplaceSPICE str2et ornativeSOFAconverter.",
    rows=rows,
    maximum_abs_inpop_minus_sofa_s=max(abs(row["inpop_minus_sofa_s"]) for row in rows),
    conclusion="Usable independent integratedgeocentric timecomparison, with knownsign/units/coverage. Threevalue comparisons cannotestablishintervalrateagreement overeachRTT orproveprismconvention. Do notsubstitute intoobservableclocks without explicitmodelreview.",
)
output = WORK / "logs/inpop17a-time-spk-audit.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
