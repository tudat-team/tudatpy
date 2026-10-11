"""Build an up-to-date Earth orientation parameter (EOP) file for tudat.

tudat reads Earth orientation from ``<earth_orientation_path>/eopc04_14_IAU2000.62-now.txt`` in the
16-column IERS 14 C04 layout (``EOPReader``: a header line starting with ``Date``, then rows of exactly
16 whitespace-separated tokens; MJD, x_p, y_p [arcsec], UT1-UTC, LOD [s], dX, dY [arcsec] in columns
4-10). The IERS stopped updating the 14 C04 series: the file at
https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt ends on 2026-01-05.

This module writes the same layout from the two current IERS products:

* EOP 20 C04 (final values, consistent with ITRF2020, about one month latency),
  https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now
* finals2000A (IERS Rapid Service / Bulletin A: rapid values and about one year of predictions),
  appended after the last 20 C04 day, https://datacenter.iers.org/data/9/finals2000A.all

Error columns are copied from 20 C04 and set to zero for finals2000A rows (tudat does not read them).
"""

from datetime import date, timedelta
from pathlib import Path
import shutil
from typing import Dict, List, Optional, Tuple

C04_20_URL = "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now"
FINALS2000A_URL = "https://datacenter.iers.org/data/9/finals2000A.all"
EOP_FILE_NAME = "eopc04_14_IAU2000.62-now.txt"

_HEADER = """               EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS OBSERVATORY)
                      INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE
                     EOP (IERS) 20 C04 TIME SERIES, IERS finals2000A (rapid + predicted) after its end
               Written in the 14 C04 layout by tudatpy.data.earth_orientation
               Sources: {c04}
                        {finals}
             FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))
##################################################################################

      Date      MJD      x          y        UT1-UTC       LOD         dX        dY        x Err     y Err   UT1-UTC Err  LOD Err     dX Err       dY Err
                         "          "           s           s          "         "           "          "          s         s            "           "
     (0h UTC)
"""


def _format_row(mjd: int, x, y, ut1, lod, dx, dy, err=(0.0,) * 6) -> str:
    d = date(1858, 11, 17) + timedelta(days=mjd)
    return (
        f"{d.year:4d}{d.month:4d}{d.day:4d}{mjd:7d}{x:11.6f}{y:11.6f}{ut1:12.7f}{lod:12.7f}{dx:11.6f}{dy:11.6f}"
        f"{err[0]:11.6f}{err[1]:11.6f}{err[2]:11.7f}{err[3]:11.7f}{err[4]:12.6f}{err[5]:12.6f}"
    )


def parse_c04_20(text: str) -> Dict[int, str]:
    """Rows of an EOP 20 C04 file (``#`` header, 21 columns), keyed by MJD, in the 14 C04 layout."""
    rows = {}
    for line in text.splitlines():
        t = line.split()
        if line.startswith("#") or len(t) != 21:
            continue
        mjd = int(round(float(t[4])))
        x, y, ut1, dx, dy, lod = (float(t[k]) for k in (5, 6, 7, 8, 9, 12))
        err = tuple(float(t[k]) for k in (13, 14, 15, 20, 16, 17))  # x, y, UT1-UTC, LOD, dX, dY
        rows[mjd] = _format_row(mjd, x, y, ut1, lod, dx, dy, err)
    return rows


def parse_finals2000a(text: str, after_mjd: int) -> Tuple[Dict[int, str], Optional[int]]:
    """Rows of a finals2000A file after ``after_mjd``, and the last MJD whose UT1-UTC is observed (flag I).

    Fixed columns (IERS readme): MJD 8-15, x 19-27, y 38-46 [arcsec], UT1-UTC 59-68 [s], LOD 80-86 [ms],
    dX 98-106, dY 117-125 [mas]; UT1-UTC flag at column 58 (I = IERS observed, P = predicted).
    Rows stop where the prediction runs out of polar motion / UT1 values.
    """
    rows, last_observed = {}, None
    for line in text.splitlines():
        if len(line) < 68 or not line[7:15].strip():
            continue
        mjd = int(float(line[7:15]))
        try:
            x, y, ut1 = float(line[18:27]), float(line[37:46]), float(line[58:68])
        except ValueError:
            break
        if line[57] == "I":
            last_observed = mjd
        if mjd <= after_mjd:
            continue
        num = lambda a, b, scale: float(line[a:b]) * scale if line[a:b].strip() else 0.0
        rows[mjd] = _format_row(
            mjd, x, y, ut1, num(79, 86, 1e-3), num(97, 106, 1e-3), num(116, 125, 1e-3)
        )
    return rows, last_observed


def build_eop_file_text(
    c04_20_text: str, finals2000a_text: Optional[str] = None
) -> Tuple[str, Dict[str, object]]:
    """Merged file content and a summary (first/last day, last 20 C04 day, last observed UT1-UTC day)."""
    rows = parse_c04_20(c04_20_text)
    if not rows:
        raise ValueError("No EOP 20 C04 rows found (expected 21 columns per data row).")
    last_c04 = max(rows)
    last_observed = None
    if finals2000a_text is not None:
        extra, last_observed = parse_finals2000a(finals2000a_text, last_c04)
        rows.update(extra)
    mjds = sorted(rows)
    if mjds != list(range(mjds[0], mjds[-1] + 1)):
        raise ValueError("Merged EOP series has missing days.")
    to_date = lambda m: (date(1858, 11, 17) + timedelta(days=m)).isoformat()
    summary = {
        "first_day": to_date(mjds[0]),
        "last_day": to_date(mjds[-1]),
        "last_c04_20_day": to_date(last_c04),
        "last_observed_ut1_day": to_date(last_observed) if last_observed else None,
        "rows": len(mjds),
    }
    header = _HEADER.format(
        c04=C04_20_URL, finals=FINALS2000A_URL if finals2000a_text is not None else "-"
    )
    return header + "\n".join(rows[m] for m in mjds) + "\n", summary


def update_earth_orientation_file(
    output_path: Optional[str] = None,
    include_rapid_and_predicted: bool = True,
    backup: bool = True,
    timeout: float = 120.0,
) -> Dict[str, object]:
    """Download EOP 20 C04 (+ finals2000A) and write tudat's EOP file.

    Parameters
    ----------
    output_path
        File to write. Default: ``eopc04_14_IAU2000.62-now.txt`` in ``tudatpy.data.get_earth_orientation_path()``,
        the file tudat reads.
    include_rapid_and_predicted
        Append finals2000A rows (rapid values and predictions) after the last 20 C04 day.
    backup
        Copy an existing file to ``<name>.bak-<YYYY-MM-DD>`` before overwriting it.
    timeout
        Download timeout per file [s].

    Returns
    -------
    dict
        Summary: output path, first/last day, last 20 C04 day, last day with observed UT1-UTC, rows.
    """
    import requests

    if output_path is None:
        from tudatpy.data import get_earth_orientation_path

        output_path = str(Path(get_earth_orientation_path()) / EOP_FILE_NAME)
    out = Path(output_path)

    def get(url: str) -> str:
        reply = requests.get(url, timeout=timeout)
        reply.raise_for_status()
        return reply.text

    text, summary = build_eop_file_text(
        get(C04_20_URL), get(FINALS2000A_URL) if include_rapid_and_predicted else None
    )
    if backup and out.is_file():
        shutil.copy2(out, out.with_name(f"{out.name}.bak-{date.today().isoformat()}"))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(out)
    summary["output_path"] = str(out)
    return summary
