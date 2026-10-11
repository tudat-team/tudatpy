"""Offline tests for tudatpy.data.earth_orientation (no downloads)."""

import pytest

from tudatpy.data import earth_orientation as eo

C04_20 = """# EOP (IERS) 20 C04 TIME SERIES  consistent with ITRF 2020 - sampled at 0h UTC
# YR  MM  DD  HH       MJD        x(")        y(")  UT1-UTC(s)       dX(")       dY(")  xrt("/day)  yrt("/day)      LOD(s)        x Er        y Er  UT1-UTC Er       dX Er       dY Er      xrt Er      yrt Er      LOD Er
2026   9   9   0  61292.00    0.201420    0.334180  -0.0003593    0.000329   -0.000100   -0.001190   -0.000310   0.0010781    0.000039    0.000047   0.0000335    0.000567    0.000459    0.000094    0.000090   0.0000088
2026   9  10   0  61293.00    0.200225    0.333865  -0.0014496    0.000335   -0.000098   -0.001205   -0.000372   0.0010906    0.000039    0.000047   0.0000335    0.000567    0.000459    0.000094    0.000090   0.0000088
"""

# Real finals2000A rows (IERS, 2026-10-08 file): 61293 overlaps 20 C04 and must be skipped; 61294-61295 observed
FINALS = [
    '26 910 61293.00 I  0.200189 0.000019  0.333881 0.000014  I-0.0015239 0.0000173  1.0982 0.0121  I     0.457    0.492    -0.198    0.116',
    '26 911 61294.00 I  0.198562 0.000019  0.333421 0.000014  I-0.0026584 0.0000171  1.1727 0.0146  I     0.456    0.492    -0.190    0.116',
    '26 912 61295.00 I  0.196426 0.000017  0.332628 0.000013  I-0.0038581 0.0000235  1.2136 0.0124  I     0.455    0.128    -0.182    0.160',
]
PREDICTED = '2610 9 61322.00 P  0.167231 0.000646  0.322270 0.000469  P-0.0302856 0.0001080                 P     0.146    0.128     0.259    0.160'   # no LOD in prediction rows


def _check_finals_columns():
    for line in FINALS:
        assert line[7:15].strip() and line[57] in "IP"


def test_c04_20_rows_in_14_layout():
    rows = eo.parse_c04_20(C04_20)
    assert sorted(rows) == [61292, 61293]
    t = rows[61293].split()
    assert len(t) == 16                       # tudat EOPReader only reads 16-token rows
    assert t[3] == "61293"
    assert float(t[4]) == pytest.approx(0.200225)
    assert float(t[6]) == pytest.approx(-0.0014496)
    assert float(t[7]) == pytest.approx(0.0010906)   # LOD [s]
    assert float(t[8]) == pytest.approx(0.000335)    # dX [arcsec]


def test_finals_appended_after_c04():
    _check_finals_columns()
    text, summary = eo.build_eop_file_text(C04_20, "\n".join(FINALS))
    lines = text.splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.split() and l.split()[0] == "Date")
    data = [l.split() for l in lines[i0 + 1:] if len(l.split()) == 16]
    assert [int(t[3]) for t in data] == [61292, 61293, 61294, 61295]
    assert float(data[1][4]) == pytest.approx(0.200225)      # 20 C04 kept on the overlapping day
    assert float(data[2][4]) == pytest.approx(0.198562)      # finals x_p [arcsec]
    assert float(data[2][6]) == pytest.approx(-0.0026584)    # finals UT1-UTC [s]
    assert float(data[2][7]) == pytest.approx(0.0011727)     # finals LOD, ms -> s
    assert float(data[2][8]) == pytest.approx(0.000456)      # finals dX, mas -> arcsec
    assert summary["last_c04_20_day"] == "2026-09-10"
    assert summary["last_observed_ut1_day"] == "2026-09-12"
    assert summary["last_day"] == "2026-09-12"
    head = lines[i0].split()
    assert head[6] == "dX" and head[7] == "dY"              # EOPReader rejects dPsi/dEps files


def test_prediction_row_without_lod():
    rows, last_observed = eo.parse_finals2000a(PREDICTED, after_mjd=0)
    t = rows[61322].split()
    assert len(t) == 16 and float(t[7]) == 0.0 and last_observed is None


def test_gap_is_rejected():
    with pytest.raises(ValueError):
        eo.build_eop_file_text(C04_20, FINALS[2])           # 61294 missing
