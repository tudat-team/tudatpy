"""Plot unchanged one-way LLR O-C by reflector and UTC; summarize descriptive trends."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

REFLECTORS = ("A11", "L1", "A14", "A15", "L2")
COLORS = {"APOLLO": "#0072B2", "GRASSE": "#D55E00", "MATERA": "#009E73"}


def read_points(path, station=None):
    """Read the shared analysis CSV schema, retaining every finite row and its input identity."""

    points = []

    with path.open(newline="") as stream:
        for line, row in enumerate(csv.DictReader(stream), 2):
            epoch = datetime.fromisoformat(row["utc_launch"].replace("Z", "+00:00"))
            epoch = (
                epoch.replace(tzinfo=timezone.utc)
                if epoch.tzinfo is None
                else epoch.astimezone(timezone.utc)
            )
            name = row.get("station") or station
            residual = float(row["one_way_residual_m"])
            sigma = float(row["one_way_sigma_m"]) if row.get("one_way_sigma_m") else None

            if not name or row["reflector"] not in REFLECTORS or not np.isfinite(residual):
                raise ValueError(f"Invalid station, reflector or residual: {path}:{line}")

            if sigma is not None and (not np.isfinite(sigma) or sigma <= 0):
                raise ValueError(f"Invalid provider normal-point sigma: {path}:{line}")

            points.append(
                dict(
                    epoch=epoch,
                    station=name,
                    reflector=row["reflector"],
                    residual=residual,
                    sigma=sigma,
                )
            )

    if not points:
        raise ValueError(f"No observations in {path}")

    return points


def nightly_means(points):
    """Form unweighted raw nightly means at their mean observation times, without altering points."""

    groups = defaultdict(list)

    for point in points:
        groups[point["epoch"].date()].append(point)

    return [
        (
            np.mean([mdates.date2num(p["epoch"]) for p in group]),
            np.mean([p["residual"] for p in group]),
        )
        for _, group in sorted(groups.items())
    ]


def describe(points):
    """Keep any unweighted nightly-mean slope separate from raw residual statistics."""

    result = []

    for station, reflector in sorted({(p["station"], p["reflector"]) for p in points}):
        group = [p for p in points if (p["station"], p["reflector"]) == (station, reflector)]
        values = np.array([p["residual"] for p in group])
        nights = nightly_means(group)
        slope = None

        # A two-point sequence, or one night, does not establish a multi-night trend.
        if len(nights) >= 3:
            days, means = np.array(nights).T
            centred = days - days.mean()
            slope = float(centred @ (means - means.mean()) / (centred @ centred) * 1000)

        result.append(
            dict(
                station=station,
                reflector=reflector,
                points=len(group),
                utc_nights=len(nights),
                first_utc=min(p["epoch"] for p in group).isoformat(),
                last_utc=max(p["epoch"] for p in group).isoformat(),
                raw_mean_mm=float(values.mean() * 1000),
                raw_rms_mm=float(np.sqrt(np.mean(values**2)) * 1000),
                raw_peak_to_peak_mm=float(np.ptp(values) * 1000),
                descriptive_nightly_mean_slope_mm_per_day=slope,
                trend_scope=(
                    "Unweighted linear description of raw nightly means; no fitted correction or period inference"
                    if slope is not None
                    else "Fewer than three UTC nights; no multi-night trend assessment, including any two-point sequence"
                ),
            )
        )

    return result


def plot(points, output, title, connect_means):
    """Save full-range reflector panels with raw points, optional provider sigmas and nightly means."""

    fig, axes = plt.subplots(5, 1, sharex=True, figsize=(11, 12), constrained_layout=True)
    stations = sorted({p["station"] for p in points})
    first, last = min(p["epoch"] for p in points), max(p["epoch"] for p in points)

    for reflector, ax in zip(REFLECTORS, axes):
        ax.axhline(0, color="0.65", linewidth=0.7)

        for index, station in enumerate(stations):
            group = sorted(
                (p for p in points if p["station"] == station and p["reflector"] == reflector),
                key=lambda p: p["epoch"],
            )

            if not group:
                continue

            color = COLORS.get(station, f"C{index}")
            dates = [p["epoch"] for p in group]
            residuals = [p["residual"] * 1000 for p in group]
            ax.scatter(
                dates,
                residuals,
                s=25,
                color=color,
                label=f"{station} raw (n={len(group)})",
                zorder=3,
            )
            supplied = [p for p in group if p["sigma"] is not None]

            if supplied:
                ax.errorbar(
                    [p["epoch"] for p in supplied],
                    [p["residual"] * 1000 for p in supplied],
                    yerr=[p["sigma"] * 1000 for p in supplied],
                    fmt="none",
                    ecolor=color,
                    alpha=0.6,
                    capsize=2,
                    label=f"{station} provider NP ±1σ",
                )

            if connect_means:
                means = np.array(nightly_means(group))
                ax.plot(
                    mdates.num2date(means[:, 0]),
                    means[:, 1] * 1000,
                    ":D",
                    color=color,
                    markersize=4,
                    linewidth=1,
                    label=f"{station} nightly raw mean (unweighted)",
                )

        ax.set_ylabel(f"{reflector}\nO−C (mm)")
        ax.grid(alpha=0.2)
        handles, labels = ax.get_legend_handles_labels()

        if handles:
            ax.legend(handles, labels, fontsize=7, loc="best")
        else:
            ax.text(0.5, 0.5, "No observations", transform=ax.transAxes, ha="center", color="0.4")

    # Include all epochs honestly; optional station detail figures preserve this full-span figure.
    span = max(mdates.date2num(last) - mdates.date2num(first), 1 / 1440)
    axes[-1].set_xlim(mdates.date2num(first) - span * 0.04, mdates.date2num(last) + span * 0.04)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=9)
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator, tz=timezone.utc))
    axes[-1].set_xlabel("UTC launch time — gaps retain their full duration")
    fig.suptitle(title + "\nRaw one-way O−C; no offset or trend removed", fontsize=12)

    for extension in ("png", "pdf"):
        fig.savefig(output.with_suffix("." + extension), dpi=170)

    plt.close(fig)


def main():
    """Plot one or more compatible CSVs, recording their hashes and all-point counts."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument(
        "--station",
        action="append",
        help="One station label per input, needed only when its station column is absent",
    )
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--title", default="Fixed-ephemeris LLR residuals")
    parser.add_argument("--nightly-means", action="store_true")
    parser.add_argument(
        "--station-details",
        action="store_true",
        help="Also retain labeled per-station time detail figures",
    )
    args = parser.parse_args()

    if args.station and len(args.station) != len(args.csv):
        parser.error("Supply exactly one --station per CSV, or use the station column")

    points, inputs = [], []

    for path, station in zip(args.csv, args.station or [None] * len(args.csv)):
        group = read_points(path, station)
        points.extend(group)
        inputs.append(
            dict(
                path=str(path),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                rows=len(group),
                station_override=station,
            )
        )

    args.output.mkdir(parents=True, exist_ok=False)
    plot(points, args.output / "residuals_full_span", args.title, args.nightly_means)

    if args.station_details:
        for station in sorted({p["station"] for p in points}):
            plot(
                [p for p in points if p["station"] == station],
                args.output / f"residuals_{station.lower()}_time_detail",
                args.title + f" — {station} time detail",
                args.nightly_means,
            )

    summary = dict(
        inputs=inputs,
        points_plotted=len(points),
        statistics=describe(points),
        plot_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation="Every raw point retained; no offset, trend or outlier removed. Slopes describe unweighted nightly means separately by station and reflector. Two months do not establish a robust monthly period. Plot timestamps use datetime display precision; source CSV bytes are unchanged.",
    )
    (args.output / "plot_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Saved {len(points)} unchanged raw points to {args.output}")


if __name__ == "__main__":
    main()
