"""Validate any JSON contract against a Tudatpy module."""

import argparse
from collections.abc import Sequence

from tudatpy.json_interface import validate_all_contracts, validate_contract_against_module


def main(arguments: Sequence[str] | None = None) -> None:
    """Run contract-to-module validation from the command line.

    Parameters
    ----------
    arguments : collections.abc.Sequence[str] or None, optional
        Command-line arguments excluding the executable name. Supplying ``None``
        makes ``argparse`` read ``sys.argv``. With no contract argument, every
        packaged contract is validated. Otherwise the first positional argument
        is the contract path, the optional second argument is its factory module,
        and ``--types`` lists modules in which enum types may be resolved.

    Returns
    -------
    None
        Validation results are written to standard output. Invalid contracts
        raise an exception and cause command-line execution to fail.

    Raises
    ------
    JSONSettingsValidationError
        If a contract is malformed or differs from the selected Tudat module.
    ModuleNotFoundError
        If an explicitly selected factory or type module cannot be imported.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("contract", nargs="?")
    parser.add_argument("module", nargs="?")
    parser.add_argument("--types", nargs="*", default=())
    arguments = parser.parse_args(arguments)
    if arguments.contract is None:
        # With no positional arguments this remains a short repository-wide
        # consistency check suitable for tests and command-line use.
        count = validate_all_contracts()
        print(f"Validated all {count} contract entries")
        return

    count = validate_contract_against_module(arguments.contract, arguments.module, arguments.types)
    target = arguments.module or "its corresponding Tudatpy module"
    print(f"Validated {count} contract entries against {target}")


if __name__ == "__main__":
    main()
