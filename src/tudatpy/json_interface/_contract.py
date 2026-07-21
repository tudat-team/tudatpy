"""Generic contract parsing, validation, and Tudat settings construction."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from importlib import import_module
from inspect import Parameter, signature
from pathlib import Path
from types import ModuleType
from typing import Any, Sequence


class JSONSettingsValidationError(ValueError):
    """Raised when a contract or settings document is invalid."""


@dataclass(frozen=True)
class ContractEntry:
    """Validated definition of one settings factory."""

    properties: dict[str, str]
    optional: dict[str, Any]


_TYPE_PATTERN = re.compile(r"^(?P<base>[A-Za-z_][A-Za-z0-9_]*)(?P<dimensions>(?:\[(?:\d+)?\])*)$")
_DIMENSION_PATTERN = re.compile(r"\[(\d*)\]")
_SCALAR_TYPES = {"bool", "float", "int", "string"}


def _reject_non_finite(value: str) -> None:
    raise JSONSettingsValidationError(f"Non-finite JSON number {value!r} is not permitted")


def read_json_object(path: str | Path, document_name: str) -> dict[str, Any]:
    """Read a JSON file and require an object at its root."""

    path = Path(path)
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream, parse_constant=_reject_non_finite)
    except json.JSONDecodeError as error:
        raise JSONSettingsValidationError(
            f"Invalid JSON in {path}: line {error.lineno}, column {error.colno}: " f"{error.msg}"
        ) from error

    return expect_object(value, document_name)


def expect_object(value: Any, path: str) -> dict[str, Any]:
    """Require a JSON object and return it with a useful error path."""

    if not isinstance(value, dict):
        raise JSONSettingsValidationError(f"{path} must be an object, got {type(value).__name__}")
    return value


def _parse_type(type_name: Any, path: str) -> tuple[str, list[int | None]]:
    if not isinstance(type_name, str):
        raise JSONSettingsValidationError(f"{path} must contain a type name as a string")

    match = _TYPE_PATTERN.fullmatch(type_name)
    if match is None:
        raise JSONSettingsValidationError(
            f"{path} contains unsupported type expression {type_name!r}"
        )

    dimensions = [
        int(size) if size else None
        for size in _DIMENSION_PATTERN.findall(match.group("dimensions"))
    ]
    return match.group("base"), dimensions


def _find_named_type(type_name: str, type_modules: tuple[ModuleType, ...]) -> Any:
    matches: list[Any] = []
    for type_module in type_modules:
        candidate = getattr(type_module, type_name, None)
        if candidate is not None:
            matches.append(candidate)

        for attribute_name in dir(type_module):
            attribute = getattr(type_module, attribute_name)
            if isinstance(attribute, ModuleType):
                candidate = getattr(attribute, type_name, None)
                if candidate is not None:
                    matches.append(candidate)

    unique_matches = {id(match): match for match in matches}
    if not unique_matches:
        raise JSONSettingsValidationError(
            f"Enum type {type_name!r} is not exposed in the configured Tudat modules"
        )
    if len(unique_matches) > 1:
        raise JSONSettingsValidationError(
            f"Enum type {type_name!r} is ambiguous in the configured Tudat modules"
        )
    return next(iter(unique_matches.values()))


def _convert_scalar(
    value: Any,
    type_name: str,
    path: str,
    type_modules: tuple[ModuleType, ...],
) -> Any:
    if type_name == "int":
        if type(value) is not int:
            raise JSONSettingsValidationError(f"{path} must be an int, got {type(value).__name__}")
        return value

    if type_name == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise JSONSettingsValidationError(f"{path} must be a float, got {type(value).__name__}")
        converted = float(value)
        if not math.isfinite(converted):
            raise JSONSettingsValidationError(f"{path} must be finite")
        return converted

    if type_name == "bool":
        if type(value) is not bool:
            raise JSONSettingsValidationError(f"{path} must be a bool, got {type(value).__name__}")
        return value

    if type_name == "string":
        if not isinstance(value, str):
            raise JSONSettingsValidationError(
                f"{path} must be a string, got {type(value).__name__}"
            )
        return value

    if not isinstance(value, str):
        raise JSONSettingsValidationError(f"{path} must name a {type_name} member as a string")
    enum_type = _find_named_type(type_name, type_modules)
    try:
        return getattr(enum_type, value)
    except AttributeError as error:
        raise JSONSettingsValidationError(
            f"{path} contains unknown {type_name} member {value!r}"
        ) from error


def _convert_value(
    value: Any,
    type_expression: Any,
    path: str,
    type_modules: tuple[ModuleType, ...],
) -> Any:
    type_name, dimensions = _parse_type(type_expression, path)

    def convert(
        current_value: Any,
        remaining_dimensions: list[int | None],
        current_path: str,
    ) -> Any:
        if not remaining_dimensions:
            return _convert_scalar(current_value, type_name, current_path, type_modules)

        if not isinstance(current_value, list):
            raise JSONSettingsValidationError(
                f"{current_path} must be an array, " f"got {type(current_value).__name__}"
            )

        expected_size = remaining_dimensions[-1]
        if expected_size is not None and len(current_value) != expected_size:
            raise JSONSettingsValidationError(
                f"{current_path} must contain {expected_size} elements, "
                f"got {len(current_value)}"
            )

        return [
            convert(item, remaining_dimensions[:-1], f"{current_path}[{index}]")
            for index, item in enumerate(current_value)
        ]

    return convert(value, dimensions, path)


def load_contract(
    contract_path: str | Path,
    factory_module: ModuleType,
    type_modules: tuple[ModuleType, ...],
) -> dict[str, ContractEntry]:
    """Load and validate a contract against a Tudat factory module."""

    document = read_json_object(contract_path, "contract")
    contract: dict[str, ContractEntry] = {}

    for factory_name, raw_entry in document.items():
        entry_path = f"contract.{factory_name}"
        entry = expect_object(raw_entry, entry_path)
        unknown_sections = entry.keys() - {"properties", "optional"}
        if unknown_sections:
            names = ", ".join(sorted(unknown_sections))
            raise JSONSettingsValidationError(f"{entry_path} contains undeclared sections: {names}")

        properties = expect_object(entry.get("properties"), f"{entry_path}.properties")
        optional = expect_object(entry.get("optional", {}), f"{entry_path}.optional")

        unknown_optional = optional.keys() - properties.keys()
        if unknown_optional:
            names = ", ".join(sorted(unknown_optional))
            raise JSONSettingsValidationError(
                f"{entry_path}.optional contains undeclared properties: {names}"
            )

        typed_properties: dict[str, str] = {}
        for property_name, type_expression in properties.items():
            property_path = f"{entry_path}.properties.{property_name}"
            base_type, _ = _parse_type(type_expression, property_path)
            if base_type not in _SCALAR_TYPES:
                _find_named_type(base_type, type_modules)
            typed_properties[property_name] = type_expression

        for property_name, default_value in optional.items():
            _convert_value(
                default_value,
                typed_properties[property_name],
                f"{entry_path}.optional.{property_name}",
                type_modules,
            )

        factory = getattr(factory_module, factory_name, None)
        if factory is None or not callable(factory):
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} is not exposed by {factory_module.__name__}"
            )

        contract[factory_name] = ContractEntry(typed_properties, optional)

    return contract


def _first_signature(factory: Any, factory_name: str) -> tuple[set[str], set[str]]:
    try:
        parameters = signature(factory).parameters.values()
    except (TypeError, ValueError):
        documentation = getattr(factory, "__doc__", "") or ""
        match = re.search(
            rf"(?:^|\n)\s*{re.escape(factory_name)}\((.*?)\)\s*->",
            documentation,
            re.DOTALL,
        )
        if match is None:
            raise JSONSettingsValidationError(
                f"Cannot read the exposed signature of factory {factory_name!r}"
            )

        declarations: list[str] = []
        start = 0
        depth = 0
        for index, character in enumerate(match.group(1)):
            if character in "([{<":
                depth += 1
            elif character in ")]}>" and depth:
                depth -= 1
            elif character == "," and depth == 0:
                declarations.append(match.group(1)[start:index])
                start = index + 1
        declarations.append(match.group(1)[start:])

        parsed: list[tuple[str, bool]] = []
        for declaration in declarations:
            declaration = declaration.strip()
            if declaration in {"", "*", "/"}:
                continue
            name_match = re.match(r"[A-Za-z_][A-Za-z0-9_]*", declaration)
            if name_match is None:
                raise JSONSettingsValidationError(
                    f"Cannot parse the exposed signature of factory {factory_name!r}"
                )

            depth = 0
            has_default = False
            for character in declaration:
                if character in "([{<":
                    depth += 1
                elif character in ")]}>" and depth:
                    depth -= 1
                elif character == "=" and depth == 0:
                    has_default = True
                    break
            parsed.append((name_match.group(), has_default))
        names = {name for name, _ in parsed}
        optional = {name for name, has_default in parsed if has_default}
        return names, optional

    names: set[str] = set()
    optional: set[str] = set()
    for parameter in parameters:
        if parameter.kind in (Parameter.VAR_POSITIONAL, Parameter.VAR_KEYWORD):
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} has no finite named-parameter interface"
            )
        if parameter.kind is Parameter.POSITIONAL_ONLY:
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} exposes positional-only parameters"
            )
        names.add(parameter.name)
        if parameter.default is not Parameter.empty:
            optional.add(parameter.name)
    return names, optional


def validate_contract_against_module(
    contract_path: str | Path,
    factory_module_name: str,
    type_module_names: Sequence[str] = (),
) -> int:
    """Validate an arbitrary contract against an importable Tudatpy module."""

    factory_module = import_module(factory_module_name)
    type_modules = tuple(
        import_module(name) for name in (type_module_names or (factory_module_name,))
    )
    contract = load_contract(contract_path, factory_module, type_modules)

    for factory_name, entry in contract.items():
        exposed, exposed_optional = _first_signature(
            getattr(factory_module, factory_name), factory_name
        )
        contracted = set(entry.properties)
        if contracted != exposed:
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} properties differ: contract-only "
                f"{sorted(contracted - exposed)}, module-only "
                f"{sorted(exposed - contracted)}"
            )
        if set(entry.optional) != exposed_optional:
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} optional properties differ: "
                f"contract-only {sorted(set(entry.optional) - exposed_optional)}, "
                f"module-only {sorted(exposed_optional - set(entry.optional))}"
            )

    return len(contract)


def create_settings_object(
    definition: Any,
    contract: dict[str, ContractEntry],
    factory_module: ModuleType,
    type_modules: tuple[ModuleType, ...],
    path: str,
) -> Any:
    """Validate a one-key settings definition and invoke its Tudat factory."""

    definition = expect_object(definition, path)
    if len(definition) != 1:
        raise JSONSettingsValidationError(f"{path} must contain exactly one settings factory")

    factory_name, arguments = next(iter(definition.items()))
    if factory_name not in contract:
        raise JSONSettingsValidationError(
            f"{path} contains unknown settings factory {factory_name!r}"
        )

    arguments = expect_object(arguments, f"{path}.{factory_name}")
    entry = contract[factory_name]

    missing = entry.properties.keys() - entry.optional.keys() - arguments.keys()
    if missing:
        names = ", ".join(sorted(missing))
        raise JSONSettingsValidationError(
            f"{path}.{factory_name} is missing required properties: {names}"
        )

    unknown = arguments.keys() - entry.properties.keys()
    if unknown:
        names = ", ".join(sorted(unknown))
        raise JSONSettingsValidationError(
            f"{path}.{factory_name} contains undeclared properties: {names}"
        )

    converted_arguments = {
        name: _convert_value(
            value,
            entry.properties[name],
            f"{path}.{factory_name}.{name}",
            type_modules,
        )
        for name, value in arguments.items()
    }

    factory = getattr(factory_module, factory_name)
    try:
        return factory(**converted_arguments)
    except (TypeError, ValueError) as error:
        raise JSONSettingsValidationError(
            f"Tudat rejected {path}.{factory_name}: {error}"
        ) from error
