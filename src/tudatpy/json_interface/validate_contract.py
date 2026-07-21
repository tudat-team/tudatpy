"""Validate any JSON contract against a Tudatpy module."""

import argparse
from collections.abc import Sequence

from tudatpy.json_interface import validate_contract_against_module


def main(arguments: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("contract")
    parser.add_argument("module")
    parser.add_argument("--types", nargs="*", default=())
    arguments = parser.parse_args(arguments)
    count = validate_contract_against_module(arguments.contract, arguments.module, arguments.types)
    print(f"Validated {count} contract entries against {arguments.module}")


if __name__ == "__main__":
    main()
