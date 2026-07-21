"""Validate any JSON contract against a Tudatpy module."""

import argparse
from collections.abc import Sequence

from tudatpy.json_interface import validate_all_contracts, validate_contract_against_module


def main(arguments: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("contract", nargs="?")
    parser.add_argument("module", nargs="?")
    parser.add_argument("--types", nargs="*", default=())
    arguments = parser.parse_args(arguments)
    if arguments.contract is None:
        count = validate_all_contracts()
        print(f"Validated all {count} contract entries")
        return

    count = validate_contract_against_module(arguments.contract, arguments.module, arguments.types)
    target = arguments.module or "its corresponding Tudatpy module"
    print(f"Validated {count} contract entries against {target}")


if __name__ == "__main__":
    main()
