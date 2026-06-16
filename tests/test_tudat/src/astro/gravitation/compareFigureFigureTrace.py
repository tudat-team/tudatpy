#!/usr/bin/env python3

import math
import sys

MODELS = ("schutz", "dmr", "oracle")


def parse_line(line):
    if not line.startswith("FFDBG|"):
        return None

    fields = {}
    for part in line.strip().split("|")[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key] = value

    required = ("case", "model", "step", "name", "shape", "values")
    if any(key not in fields for key in required):
        raise RuntimeError(f"Malformed FFDBG line: {line.rstrip()}")

    values = tuple(float(value) for value in fields["values"].split(",") if value)
    return fields["case"], fields["model"], fields["step"], fields["name"], fields["shape"], values


def max_abs_diff(left, right):
    return max(abs(a - b) for a, b in zip(left, right))


def rel_diff(left, right):
    scale = max(1.0e-300, max(abs(value) for value in right))
    return max_abs_diff(left, right) / scale


def format_values(values):
    return ",".join(f"{value:.17e}" for value in values)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: compareFigureFigureTrace.py <trace-output-file>")

    traces = {}
    step_order = []
    with open(sys.argv[1], "r", encoding="utf-8") as handle:
        for line in handle:
            parsed = parse_line(line)
            if parsed is None:
                continue
            case, model, step, name, shape, values = parsed
            traces[(case, step, name, model)] = (shape, values)
            if (step, name) not in step_order:
                step_order.append((step, name))

    cases = sorted({case for case, _, _, _ in traces})
    if not cases:
        raise SystemExit("No FFDBG canonical trace lines found")

    tolerance = 1.0e-20
    max_abs = 0.0
    max_rel = 0.0
    max_case = ""
    max_step = ""
    max_name = ""

    for case in cases:
        for step, name in sorted(step_order):
            entries = {}
            missing = []
            for model in MODELS:
                key = (case, step, name, model)
                if key not in traces:
                    missing.append(model)
                else:
                    entries[model] = traces[key]
            if missing:
                print(
                    f"FIRST_DIVERGENCE|case={case}|step={step}|name={name}|schutz=<missing:{','.join(missing)}>"
                    f"|dmr=<missing:{','.join(missing)}>|oracle=<missing:{','.join(missing)}>"
                    "|absDiffSchutz=nan|absDiffDmr=nan|relDiffSchutz=nan|relDiffDmr=nan"
                )
                return 2

            oracle_shape, oracle_values = entries["oracle"]
            schutz_shape, schutz_values = entries["schutz"]
            dmr_shape, dmr_values = entries["dmr"]
            if schutz_shape != oracle_shape or dmr_shape != oracle_shape:
                print(
                    f"FIRST_DIVERGENCE|case={case}|step={step}|name={name}|schutz=shape:{schutz_shape}"
                    f"|dmr=shape:{dmr_shape}|oracle=shape:{oracle_shape}"
                    "|absDiffSchutz=nan|absDiffDmr=nan|relDiffSchutz=nan|relDiffDmr=nan"
                )
                return 2
            if len(schutz_values) != len(oracle_values) or len(dmr_values) != len(oracle_values):
                print(
                    f"FIRST_DIVERGENCE|case={case}|step={step}|name={name}|schutz=len:{len(schutz_values)}"
                    f"|dmr=len:{len(dmr_values)}|oracle=len:{len(oracle_values)}"
                    "|absDiffSchutz=nan|absDiffDmr=nan|relDiffSchutz=nan|relDiffDmr=nan"
                )
                return 2

            abs_schutz = max_abs_diff(schutz_values, oracle_values)
            abs_dmr = max_abs_diff(dmr_values, oracle_values)
            rel_schutz = rel_diff(schutz_values, oracle_values)
            rel_dmr = rel_diff(dmr_values, oracle_values)
            current_abs = max(abs_schutz, abs_dmr)
            current_rel = max(rel_schutz, rel_dmr)
            if current_abs > max_abs:
                max_abs = current_abs
                max_rel = current_rel
                max_case = case
                max_step = step
                max_name = name

            if not math.isfinite(current_abs) or current_abs > tolerance:
                print(
                    f"FIRST_DIVERGENCE|case={case}|step={step}|name={name}"
                    f"|schutz={format_values(schutz_values)}|dmr={format_values(dmr_values)}"
                    f"|oracle={format_values(oracle_values)}|absDiffSchutz={abs_schutz:.17e}"
                    f"|absDiffDmr={abs_dmr:.17e}|relDiffSchutz={rel_schutz:.17e}|relDiffDmr={rel_dmr:.17e}"
                )
                return 1

    print(
        f"TRACE_OK|maxAbsDiff={max_abs:.17e}|maxRelDiff={max_rel:.17e}"
        f"|case={max_case}|step={max_step}|name={max_name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
