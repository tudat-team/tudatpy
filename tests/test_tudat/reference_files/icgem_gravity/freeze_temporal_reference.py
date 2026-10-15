"""Freeze independent scalar G10 references from printed coefficients, without Tudat.

Run from any directory with Python's standard library. Dates here use the
explicit test conventions, not an inferred EIGEN time-scale/date interpretation.
"""

from datetime import datetime
from decimal import Decimal, localcontext
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
D = Decimal
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
EPS = D(2) ** -52
JULIAN_YEAR = D(31557600)


def seconds(date):
    """Uniform source-calendar seconds since J2000 noon, with no implicit time-scale conversion."""
    return D(int((datetime.fromisoformat(date) - datetime(2000, 1, 1, 12)).total_seconds()))


def decimal_year(date):
    """Exact rational Gregorian year coordinate from the explicitly declared source calendar."""
    time = datetime.fromisoformat(date)
    first = datetime(time.year, 1, 1)
    last = datetime(time.year + 1, 1, 1)
    return D(time.year) + D(int((time - first).total_seconds())) / D(
        int((last - first).total_seconds())
    )


def sincos(angle):
    """Evaluate both quadratures by 80-digit Decimal Taylor series."""
    angle %= 2 * PI
    if angle > PI:
        angle -= 2 * PI
    sine, cosine, sterm, cterm = angle, D(1), angle, D(1)
    for k in range(1, 180):
        sterm *= -(angle * angle) / D(2 * k * (2 * k + 1))
        cterm *= -(angle * angle) / D((2 * k - 1) * 2 * k)
        sine += sterm
        cosine += cterm
    return sine, cosine


def gamma(count):
    """Conservative binary64 arithmetic bound using epsilon rather than half epsilon."""
    return count * EPS / (1 - count * EPS)


def half_printed_unit(token):
    """Half the final printed decimal unit, separate from computational roundoff."""
    return D(5).scaleb(D(token).as_tuple().exponent - 1)


source = ROOT / "eigen6s4_v2_c20_1950_2006.gfc"
assert (
    hashlib.sha256(source.read_bytes()).hexdigest()
    == "308a3f950dab2901aaf691f5d39bb754852b00408442852f3c5e4168fd3569ed"
)
assert (
    hashlib.sha256((ROOT / "controlled_temporal_c21.gfc").read_bytes()).hexdigest()
    == "5ad9bc4f2d63cf857c817b44ef63198faad7dd80428cf8b2581562a392f5eb2a"
)
rows = [
    line.split()
    for line in source.read_text().splitlines()
    if line.startswith(("gfct", "trnd", "acos", "asin"))
]

# Select published records explicitly. The periodic origin remains 1950, independent of the active 2000/2001 baseline.
eigen_cases = [
    ("2000-07-02 00:00:00", 0, "20000101.0000", "19500101.0000"),
    ("2001-07-02 00:00:00", 0, "20010101.0000", "19500101.0000"),
    ("2003-01-01 00:00:00", -1, "20020101.0000", "19500101.0000"),
    ("2003-01-01 00:00:00", 0, "20030101.0000", "20030101.0000"),
    ("2003-07-02 00:00:00", 0, "20030101.0000", "20030101.0000"),
    ("2004-12-26 01:00:00", -1, "20040101.0000", "20040101.0000"),
    ("2004-12-26 01:00:00", 0, "20041226.0060", "20041226.0060"),
]


def record_date(token):
    """Resolve the one exact-token caller override; no rollover or fractional-day inference."""
    if token == "20041226.0060":
        return "2004-12-26 01:00:00"
    return f"{token[:4]}-{token[4:6]}-{token[6:8]} 00:00:00"


result = [
    "# kind source_seconds side fractional_year C S computation_budget_C computation_budget_S abs_rate_C abs_rate_S printed_rounding_C printed_rounding_S"
]
with localcontext() as context:
    context.prec = 80
    for date, side, baseline_origin, periodic_origin in eigen_cases:
        baseline = next(row for row in rows if row[0] == "gfct" and row[7] == baseline_origin)
        trend = next(row for row in rows if row[0] == "trnd" and row[7] == baseline_origin)
        elapsed = (seconds(date) - seconds(record_date(baseline_origin))) / JULIAN_YEAR
        period_elapsed = (seconds(date) - seconds(record_date(periodic_origin))) / JULIAN_YEAR
        value = D(baseline[3]) + D(trend[3]) * elapsed
        magnitude = abs(D(baseline[3])) + abs(D(trend[3]) * elapsed)
        printed = half_printed_unit(baseline[3]) + half_printed_unit(trend[3]) * abs(elapsed)
        rate = abs(D(trend[3])) / JULIAN_YEAR
        phase_budget = D(0)
        for row in rows:
            if row[0] not in ("asin", "acos") or row[7] != periodic_origin:
                continue
            phase = 2 * PI * period_elapsed / D(row[9])
            sine, cosine = sincos(phase)
            amplitude = D(row[3])
            value += amplitude * (sine if row[0] == "asin" else cosine)
            magnitude += abs(amplitude)
            printed += half_printed_unit(row[3])
            rate += abs(amplitude) * 2 * PI / (D(row[9]) * JULIAN_YEAR)
            phase_budget += abs(amplitude) * (gamma(12) * abs(phase) + 4 * EPS)
        nominal = D.from_float(-4.84165e-4)
        budget = gamma(32) * (magnitude + abs(nominal)) + phase_budget
        value -= nominal
        result.append(
            f"EIGEN {seconds(date)} {side} 0 {value:.35e} 0 {budget:.17e} 0 {rate:.17e} 0 {printed:.17e} 0"
        )

    # Analytic unnormalized C21/S21 with an independent HHMM trend origin and annual/semiannual quadratures.
    # N21 = sqrt(5/3), so normalization multiplies by sqrt(3/5), followed by the source/output GM and squared-radius ratios.
    scale = (D(3) / D(5)).sqrt() * D(4) / D(5) * (D(8) / D(7)) ** 2
    for fractional in (False, True):
        for date in ("2000-07-02 12:34:00", "2001-07-02 12:34:00"):
            if fractional:
                elapsed = decimal_year(date) - decimal_year("2000-01-01 00:00:00")
                trend_elapsed = decimal_year(date) - decimal_year("2000-01-01 12:34:00")
                year_seconds = seconds(
                    f"{datetime.fromisoformat(date).year+1}-01-01 00:00:00"
                ) - seconds(f"{datetime.fromisoformat(date).year}-01-01 00:00:00")
            else:
                elapsed = (seconds(date) - seconds("2000-01-01 00:00:00")) / JULIAN_YEAR
                trend_elapsed = (seconds(date) - seconds("2000-01-01 12:34:00")) / JULIAN_YEAR
                year_seconds = JULIAN_YEAR
            values, budgets, rates = [], [], []
            for baseline, trend, ac, asi, sc, ssi, nominal in (
                ("2e-6", "4e-8", "5e-8", "-6e-8", "2e-8", "1e-8", 1e-6),
                ("-3e-6", "7e-8", "-2e-8", "3e-8", "4e-8", "-5e-8", -2e-6),
            ):
                value = D(baseline) + D(trend) * trend_elapsed
                magnitude = abs(D(baseline)) + abs(D(trend) * trend_elapsed)
                rate = abs(D(trend)) / year_seconds
                phase_budget = D(0)
                for period, acos, asin in ((D(1), D(ac), D(asi)), (D("0.5"), D(sc), D(ssi))):
                    phase = 2 * PI * elapsed / period
                    sine, cosine = sincos(phase)
                    value += acos * cosine + asin * sine
                    magnitude += abs(acos) + abs(asin)
                    rate += (abs(acos) + abs(asin)) * 2 * PI / (period * year_seconds)
                    phase_budget += (abs(acos) + abs(asin)) * (gamma(12) * abs(phase) + 4 * EPS)
                value = scale * value - D.from_float(nominal)
                values.append(value)
                # Include normalization, GM/radius scaling and nominal subtraction, plus calendar-coordinate cancellation near year 2000.
                budgets.append(
                    gamma(32) * (scale * magnitude + abs(D.from_float(nominal)))
                    + scale * phase_budget
                    + EPS * D(2002) * scale * D("1e-6")
                )
                rates.append(scale * rate)
            result.append(
                f"CONTROLLED {seconds(date)} 0 {int(fractional)} {values[0]:.35e} {values[1]:.35e} {budgets[0]:.17e} {budgets[1]:.17e} {rates[0]:.17e} {rates[1]:.17e} 0 0"
            )
(ROOT / "temporal_reference.txt").write_text("\n".join(result) + "\n")
print("Wrote", len(result) - 1, "independent scalar reference rows.")
