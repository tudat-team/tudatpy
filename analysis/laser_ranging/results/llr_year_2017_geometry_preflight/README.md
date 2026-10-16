All **1582 original source points** (311 APOLLO, 1271 Grasse) were screened
before the annual rerun, including the three subsequently excluded identities.
No point was selected by residual, and no observation was solved. This is a
geometry/domain check, not an atmospheric-accuracy test.

The screen uses the pinned coherent DE440 kernels, direct SPICE states, the
existing native GCRS-to-ITRS rotation with ICRF base axes, WGS84 geodetic Up and
the published INPOP17a station nominal XYZ. Station drift, deformation and
spatial mapping are omitted at this screening stage. UTC is converted with the
native Earth converter; approximate bounce is transmit TDB plus half the source
RTT, and approximate receive is transmit plus the source RTT. At each terminal,
the rotated station-to-bounce-Moon-centre direction is projected onto geodetic
Up. There are no annual cubic tables or atmosphere calls.

The predeclared candidate threshold is **4°**, leaving a 1° margin above the
native FCULa 3° minimum for lunar radius and approximate event geometry.
Exactly two points fall below that threshold:

| Grasse A15 source identity | Verbatim UTC | Source RTT (s) | Minimum Moon-centre elevation (°) |
| --- | --- | ---: | ---: |
| `apollo15_201702.npt:168` | 2017-02-20T13:23:08.387557778100 | 0.000000024976 | −18.741654 |
| `apollo15_201702.npt:211` | 2017-02-20T13:23:20.476420038600 | 0.000000017271 | −18.776436 |

These are the independently identified impossible lunar-path RTT records,
excluded explicitly by source identity; low elevation is additional evidence,
not the exclusion criterion. Their approximate nanosecond RTT endpoints should
not be interpreted as a physical lunar light path.

Every other source point has both approximate terminal elevations ≥12.560756°.
The minimum is Grasse L2, `luna21_201705.npt:69`,
2017-05-02T23:25:02.961114373500 UTC: transmit 12.567662°, receive 12.560756°.
Thus every one of the 1579 retained points has a Moon-centre margin exceeding 9°
above the 3° domain boundary, comfortably beyond the predeclared 1° screen margin.
No reflector-level refinement, new observation control or guard change was
needed after review. This does not claim exact observable endpoint elevations.

[summary.json](summary.json) retains the two candidate rows, minimum other row,
epochs, kernel/input/profile hashes and the SHA256 of all 1582 screened rows.
The full CSV, original summary and one-off script remain in the external work
cache at `results/llr-year-2017-geometry-preflight/` and
`verification/llr_year_geometry_preflight.py`; their hashes pin the audit without
duplicating a large diagnostic table or adding another analysis framework.
