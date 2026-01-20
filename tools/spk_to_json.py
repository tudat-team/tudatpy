#!/usr/bin/env python3
"""
SPK to JSON Converter for WASM Ephemeris Support

Converts binary SPICE SPK kernel data to JSON format that can be loaded
in WebAssembly environments where binary kernel loading is not supported.

Usage:
    python spk_to_json.py <spk_file> <target> <observer> <frame> <start_date> <end_date> [options]

Example:
    python spk_to_json.py de430.bsp Earth Sun J2000 2020-01-01 2030-01-01 --step-hours=1 -o earth_sun_j2000.json

Dependencies:
    pip install spiceypy numpy

Author: Tudat Team
"""

import argparse
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

try:
    import spiceypy as spice
    import numpy as np
except ImportError as e:
    print(f"Error: Missing dependency - {e}")
    print("Install with: pip install spiceypy numpy")
    sys.exit(1)


def parse_date(date_str: str) -> datetime:
    """Parse date string in various formats."""
    formats = [
        "%Y-%m-%d",
        "%Y-%m-%dT%H:%M:%S",
        "%Y/%m/%d",
        "%d-%b-%Y",
    ]
    for fmt in formats:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    raise ValueError(f"Unable to parse date: {date_str}")


def datetime_to_et(dt: datetime) -> float:
    """Convert datetime to SPICE ephemeris time (seconds since J2000)."""
    date_str = dt.strftime("%Y-%m-%dT%H:%M:%S")
    return spice.str2et(date_str)


def extract_ephemeris(
    spk_file: str,
    target: str,
    observer: str,
    frame: str,
    start_date: str,
    end_date: str,
    step_hours: float = 1.0,
    lsk_file: str = None,
) -> dict:
    """
    Extract ephemeris data from SPK file.

    Args:
        spk_file: Path to SPK binary kernel file
        target: Target body name (e.g., "Earth", "Mars")
        observer: Observer body name (e.g., "Sun", "Earth")
        frame: Reference frame (e.g., "J2000", "ECLIPJ2000")
        start_date: Start date string
        end_date: End date string
        step_hours: Time step in hours
        lsk_file: Optional leap seconds kernel (uses built-in if not provided)

    Returns:
        Dictionary with metadata and states array
    """
    # Load kernels
    spice.furnsh(spk_file)

    # Load leap seconds kernel if provided, otherwise try common locations
    if lsk_file:
        spice.furnsh(lsk_file)
    else:
        # Try to find a leap seconds kernel
        common_lsk_paths = [
            Path(spk_file).parent / "naif0012.tls",
            Path(spk_file).parent / "latest_leapseconds.tls",
            Path.home() / ".spice" / "naif0012.tls",
        ]
        lsk_loaded = False
        for lsk_path in common_lsk_paths:
            if lsk_path.exists():
                spice.furnsh(str(lsk_path))
                lsk_loaded = True
                print(f"Loaded leap seconds kernel: {lsk_path}")
                break

        if not lsk_loaded:
            print("Warning: No leap seconds kernel found. Using approximate conversion.")
            print("For accurate results, provide --lsk-file or place naif0012.tls in the SPK directory.")

    # Parse dates and convert to ET
    start_dt = parse_date(start_date)
    end_dt = parse_date(end_date)

    start_et = datetime_to_et(start_dt)
    end_et = datetime_to_et(end_dt)

    step_seconds = step_hours * 3600.0

    # Generate epochs
    num_points = int((end_et - start_et) / step_seconds) + 1
    epochs = np.linspace(start_et, end_et, num_points)

    print(f"Extracting {num_points} states from {start_date} to {end_date}")
    print(f"Target: {target}, Observer: {observer}, Frame: {frame}")

    # Extract states
    states = []
    for i, epoch in enumerate(epochs):
        if i % 10000 == 0:
            progress = (i / num_points) * 100
            print(f"Progress: {progress:.1f}%", end="\r")

        try:
            # spkezr returns state in km and km/s
            state, _ = spice.spkezr(target, epoch, frame, "NONE", observer)

            # Convert to meters and m/s for Tudat compatibility
            state_m = [
                epoch,                  # seconds since J2000
                state[0] * 1000.0,      # x (m)
                state[1] * 1000.0,      # y (m)
                state[2] * 1000.0,      # z (m)
                state[3] * 1000.0,      # vx (m/s)
                state[4] * 1000.0,      # vy (m/s)
                state[5] * 1000.0,      # vz (m/s)
            ]
            states.append(state_m)
        except Exception as e:
            print(f"\nError at epoch {epoch}: {e}")
            continue

    print(f"\nExtracted {len(states)} states successfully")

    # Clean up
    spice.kclear()

    return {
        "metadata": {
            "target": target,
            "observer": observer,
            "frame": frame,
            "startEpoch": start_et,
            "endEpoch": end_et,
            "stepSize": step_seconds,
            "numStates": len(states),
            "units": {
                "time": "seconds since J2000",
                "position": "m",
                "velocity": "m/s"
            },
            "source": Path(spk_file).name,
            "generatedAt": datetime.utcnow().isoformat() + "Z"
        },
        "states": states
    }


def main():
    parser = argparse.ArgumentParser(
        description="Convert SPK binary kernel to JSON for WASM ephemeris support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Extract Earth ephemeris relative to Sun for 2020-2030
    python spk_to_json.py de430.bsp Earth Sun J2000 2020-01-01 2030-01-01 -o earth_sun.json

    # Extract multiple planets with 30-minute resolution
    python spk_to_json.py de430.bsp Mars "Sun" J2000 2020-01-01 2025-01-01 --step-hours=0.5

    # Batch extract all major planets
    for planet in Mercury Venus Earth Mars Jupiter Saturn Uranus Neptune; do
        python spk_to_json.py de430.bsp $planet Sun J2000 2020-01-01 2030-01-01 -o ${planet,,}_sun.json
    done
"""
    )

    parser.add_argument("spk_file", help="Path to SPK binary kernel file")
    parser.add_argument("target", help="Target body name (e.g., Earth, Mars)")
    parser.add_argument("observer", help="Observer body name (e.g., Sun)")
    parser.add_argument("frame", help="Reference frame (e.g., J2000, ECLIPJ2000)")
    parser.add_argument("start_date", help="Start date (YYYY-MM-DD)")
    parser.add_argument("end_date", help="End date (YYYY-MM-DD)")

    parser.add_argument("--step-hours", type=float, default=1.0,
                        help="Time step in hours (default: 1.0)")
    parser.add_argument("--lsk-file", help="Path to leap seconds kernel (.tls)")
    parser.add_argument("-o", "--output", help="Output JSON file path")
    parser.add_argument("--compact", action="store_true",
                        help="Output compact JSON (no indentation)")
    parser.add_argument("--precision", type=int, default=6,
                        help="Decimal precision for state values (default: 6)")

    args = parser.parse_args()

    # Validate SPK file exists
    if not Path(args.spk_file).exists():
        print(f"Error: SPK file not found: {args.spk_file}")
        sys.exit(1)

    # Extract ephemeris
    data = extract_ephemeris(
        spk_file=args.spk_file,
        target=args.target,
        observer=args.observer,
        frame=args.frame,
        start_date=args.start_date,
        end_date=args.end_date,
        step_hours=args.step_hours,
        lsk_file=args.lsk_file,
    )

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        target_lower = args.target.lower().replace(" ", "_")
        observer_lower = args.observer.lower().replace(" ", "_")
        frame_lower = args.frame.lower()
        output_path = Path(f"{target_lower}_{observer_lower}_{frame_lower}.json")

    # Custom JSON encoder for controlled precision
    class PrecisionEncoder(json.JSONEncoder):
        def encode(self, obj):
            if isinstance(obj, list) and len(obj) > 0 and isinstance(obj[0], list):
                # Format states array with controlled precision
                lines = []
                for state in obj:
                    formatted = [f"{v:.{args.precision}f}" if isinstance(v, float) else str(v) for v in state]
                    lines.append("[" + ",".join(formatted) + "]")
                return "[\n" + ",\n".join(lines) + "\n]"
            return super().encode(obj)

    # Write JSON
    with open(output_path, "w") as f:
        if args.compact:
            json.dump(data, f, separators=(",", ":"))
        else:
            # Write metadata with normal formatting
            f.write('{\n"metadata": ')
            json.dump(data["metadata"], f, indent=2)
            f.write(',\n"states": ')
            # Write states with controlled precision
            encoder = PrecisionEncoder()
            f.write(encoder.encode(data["states"]))
            f.write('\n}')

    file_size = output_path.stat().st_size / (1024 * 1024)
    print(f"\nOutput written to: {output_path}")
    print(f"File size: {file_size:.2f} MB")
    print(f"States: {len(data['states'])}")


if __name__ == "__main__":
    main()
