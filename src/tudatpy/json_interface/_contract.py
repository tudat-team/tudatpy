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
from typing import Any, Mapping, Sequence


class JSONSettingsValidationError(ValueError):
    """Raised when a contract or settings document is invalid."""


@dataclass(frozen=True)
class ContractEntry:
    """Validated definition of one settings factory."""

    properties: dict[str, str]
    optional: dict[str, Any]
    unsupported: frozenset[str]


_TYPE_PATTERN = re.compile(
    r"^(?P<base>[A-Za-z_][A-Za-z0-9_.]*)(?P<containers>(?:(?:\[(?:\d+)?\])|\{\})*)$"
)
_CONTAINER_PATTERN = re.compile(r"\[(\d*)\]|(\{\})")
_BUILTIN_TYPES = {"any", "bool", "complex", "float", "int", "null", "object", "string"}
CONTRACT_ROOT = Path(__file__).parent / "contracts"


def _reject_non_finite(value: str) -> None:
    raise JSONSettingsValidationError(f"Non-finite JSON number {value!r} is not permitted")


def _read_json_value(path: Path) -> Any:
    try:
        with path.open("r", encoding="utf-8") as stream:
            return json.load(stream, parse_constant=_reject_non_finite)
    except OSError as error:
        raise JSONSettingsValidationError(f"Cannot read JSON file {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise JSONSettingsValidationError(
            f"Invalid JSON in {path}: line {error.lineno}, column {error.colno}: " f"{error.msg}"
        ) from error


def _resolve_references(value: Any, source: Path, active: frozenset[Path]) -> Any:
    if isinstance(value, list):
        return [_resolve_references(item, source, active) for item in value]
    if not isinstance(value, dict):
        return value
    if "$ref" in value:
        if set(value) != {"$ref"} or not isinstance(value["$ref"], str):
            raise JSONSettingsValidationError(
                f"A JSON reference in {source} must contain only a string $ref"
            )
        referenced = (source.parent / value["$ref"]).resolve()
        if referenced in active:
            raise JSONSettingsValidationError(f"Circular JSON reference to {referenced}")
        return _resolve_references(_read_json_value(referenced), referenced, active | {referenced})
    return {key: _resolve_references(item, source, active) for key, item in value.items()}


def read_json_object(path: str | Path, document_name: str) -> dict[str, Any]:
    """Read a JSON object and resolve relative file references."""

    path = Path(path).resolve()
    value = _resolve_references(_read_json_value(path), path, frozenset({path}))
    return expect_object(value, document_name)


def expect_object(value: Any, path: str) -> dict[str, Any]:
    """Require a JSON object and return it with a useful error path."""

    if not isinstance(value, dict):
        raise JSONSettingsValidationError(f"{path} must be an object, got {type(value).__name__}")
    return value


def _parse_type(type_name: Any, path: str) -> tuple[str, list[int | None | str]]:
    if not isinstance(type_name, str):
        raise JSONSettingsValidationError(f"{path} must contain a type name as a string")

    match = _TYPE_PATTERN.fullmatch(type_name)
    if match is None:
        raise JSONSettingsValidationError(
            f"{path} contains unsupported type expression {type_name!r}"
        )

    containers = [
        "dict" if dictionary else (int(size) if size else None)
        for size, dictionary in _CONTAINER_PATTERN.findall(match.group("containers"))
    ]
    return match.group("base"), containers


def _find_named_type(type_name: str, type_modules: tuple[ModuleType, ...]) -> Any:
    matches: list[Any] = []
    pending = list(type_modules)
    visited: set[int] = set()
    while pending:
        type_module = pending.pop()
        if id(type_module) in visited:
            continue
        visited.add(id(type_module))
        candidate = getattr(type_module, type_name, None)
        if candidate is not None:
            matches.append(candidate)

        for attribute_name in dir(type_module):
            attribute = getattr(type_module, attribute_name)
            if isinstance(attribute, ModuleType) and attribute.__name__.startswith("tudatpy"):
                pending.append(attribute)

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
    contract_root: Path,
    factory_context: Mapping[str, Mapping[str, Any]] | None,
) -> Any:
    if value is None:
        if type_name in {"any", "object", "null"}:
            return None
        raise JSONSettingsValidationError(f"{path} must not be null")
    if type_name in {"any", "object"}:
        return value
    if type_name == "null":
        raise JSONSettingsValidationError(f"{path} must be null")
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

    if type_name == "complex":
        if (
            not isinstance(value, list)
            or len(value) != 2
            or any(isinstance(item, bool) or not isinstance(item, (int, float)) for item in value)
        ):
            raise JSONSettingsValidationError(f"{path} must be a [real, imaginary] number pair")
        converted = complex(value[0], value[1])
        if not math.isfinite(converted.real) or not math.isfinite(converted.imag):
            raise JSONSettingsValidationError(f"{path} must be finite")
        return converted

    if "." in type_name:
        contract_path = contract_root.joinpath(*type_name.split(".")).with_suffix(".json")
        factory_module = import_module(module_name_from_contract(contract_path))
        contract = load_contract(contract_path, factory_module, type_modules)
        return create_settings_object(
            value,
            contract,
            factory_module,
            type_modules,
            path,
            contract_root,
            factory_context,
        )

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
    contract_root: Path = CONTRACT_ROOT,
    factory_context: Mapping[str, Mapping[str, Any]] | None = None,
) -> Any:
    if isinstance(type_expression, str) and "|" in type_expression:
        errors = []
        for alternative in type_expression.split("|"):
            try:
                return _convert_value(
                    value,
                    alternative,
                    path,
                    type_modules,
                    contract_root,
                    factory_context,
                )
            except JSONSettingsValidationError as error:
                errors.append(str(error))
        raise JSONSettingsValidationError(
            f"{path} does not match {type_expression!r}: {'; '.join(errors)}"
        )

    type_name, dimensions = _parse_type(type_expression, path)

    def convert(
        current_value: Any,
        remaining_dimensions: list[int | None | str],
        current_path: str,
    ) -> Any:
        if not remaining_dimensions:
            return _convert_scalar(
                current_value,
                type_name,
                current_path,
                type_modules,
                contract_root,
                factory_context,
            )

        expected_size = remaining_dimensions[-1]
        if expected_size == "dict":
            if not isinstance(current_value, dict):
                raise JSONSettingsValidationError(
                    f"{current_path} must be an object, got {type(current_value).__name__}"
                )
            return {
                key: convert(item, remaining_dimensions[:-1], f"{current_path}.{key}")
                for key, item in current_value.items()
            }

        if not isinstance(current_value, list):
            raise JSONSettingsValidationError(
                f"{current_path} must be an array, " f"got {type(current_value).__name__}"
            )

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
        unknown_sections = entry.keys() - {"properties", "optional", "unsupported"}
        if unknown_sections:
            names = ", ".join(sorted(unknown_sections))
            raise JSONSettingsValidationError(f"{entry_path} contains undeclared sections: {names}")

        properties = expect_object(entry.get("properties"), f"{entry_path}.properties")
        optional = expect_object(entry.get("optional", {}), f"{entry_path}.optional")
        unsupported = entry.get("unsupported", [])
        if (
            not isinstance(unsupported, list)
            or any(not isinstance(name, str) for name in unsupported)
            or len(unsupported) != len(set(unsupported))
        ):
            raise JSONSettingsValidationError(
                f"{entry_path}.unsupported must be an array of unique property names"
            )

        unknown_optional = optional.keys() - properties.keys()
        if unknown_optional:
            names = ", ".join(sorted(unknown_optional))
            raise JSONSettingsValidationError(
                f"{entry_path}.optional contains undeclared properties: {names}"
            )

        unknown_unsupported = set(unsupported) - properties.keys()
        if unknown_unsupported:
            names = ", ".join(sorted(unknown_unsupported))
            raise JSONSettingsValidationError(
                f"{entry_path}.unsupported contains undeclared properties: {names}"
            )

        typed_properties: dict[str, str] = {}
        for property_name, type_expression in properties.items():
            property_path = f"{entry_path}.properties.{property_name}"
            if not isinstance(type_expression, str):
                raise JSONSettingsValidationError(f"{property_path} must be a string")
            alternatives = type_expression.split("|")
            for alternative in alternatives:
                base_type, _ = _parse_type(alternative, property_path)
                if base_type not in _BUILTIN_TYPES and "." not in base_type:
                    _find_named_type(base_type, type_modules)
            typed_properties[property_name] = type_expression

        for property_name, default_value in optional.items():
            # JSON null represents defaults that are None or not finite in Python.
            if default_value is not None:
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

        contract[factory_name] = ContractEntry(typed_properties, optional, frozenset(unsupported))

    return contract


def _first_signature(factory: Any, factory_name: str) -> tuple[set[str], set[str]]:
    try:
        parameters = signature(factory).parameters.values()
    except (TypeError, ValueError):
        documentation = getattr(factory, "__doc__", "") or ""
        prefix = (
            rf"(?:^|\n)\s*1\.\s*{re.escape(factory_name)}"
            if "Overloaded function." in documentation
            else rf"(?:^|\n)\s*{re.escape(factory_name)}"
        )
        match = re.search(
            rf"{prefix}\((.*?)\)\s*->",
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
    factory_module_name: str | None = None,
    type_module_names: Sequence[str] = (),
) -> int:
    """Validate an arbitrary contract against an importable Tudatpy module."""

    factory_module_name = factory_module_name or module_name_from_contract(contract_path)
    factory_module = import_module(factory_module_name)
    default_type_modules = (
        ("tudatpy.dynamics",)
        if factory_module_name.startswith("tudatpy.dynamics.")
        else (factory_module_name,)
    )
    type_modules = tuple(
        import_module(name) for name in (type_module_names or default_type_modules)
    )
    contract = load_contract(contract_path, factory_module, type_modules)

    for factory_name, entry in contract.items():
        exposed, exposed_optional = _first_signature(
            getattr(factory_module, factory_name), factory_name
        )
        contracted = set(entry.properties)
        required = exposed - exposed_optional
        if not required <= contracted or not contracted <= exposed:
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} properties differ: contract-only "
                f"{sorted(contracted - exposed)}, module-only "
                f"{sorted(required - contracted)}"
            )
        contracted_optional = contracted & exposed_optional
        if set(entry.optional) != contracted_optional:
            raise JSONSettingsValidationError(
                f"Factory {factory_name!r} optional properties differ: "
                f"contract-only {sorted(set(entry.optional) - contracted_optional)}, "
                f"module-only {sorted(contracted_optional - set(entry.optional))}"
            )

    return len(contract)


def validate_all_contracts(contract_root: str | Path = CONTRACT_ROOT) -> int:
    """Validate every discovered contract against its corresponding module."""

    total = 0
    for contract_path in sorted(Path(contract_root).rglob("*.json")):
        total += validate_contract_against_module(
            contract_path,
            type_module_names=("tudatpy.dynamics",),
        )
    return total


def create_settings_object(
    definition: Any,
    contract: dict[str, ContractEntry],
    factory_module: ModuleType,
    type_modules: tuple[ModuleType, ...],
    path: str,
    contract_root: Path = CONTRACT_ROOT,
    factory_context: Mapping[str, Mapping[str, Any]] | None = None,
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

    supplied_unsupported = arguments.keys() & entry.unsupported
    required_unsupported = entry.unsupported - entry.optional.keys()
    unsupported = supplied_unsupported or required_unsupported
    if unsupported:
        names = ", ".join(sorted(unsupported))
        raise JSONSettingsValidationError(
            f"{path}.{factory_name} is not yet supported by the JSON interface; "
            f"unsupported properties: {names}"
        )

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

    values = {**entry.optional, **arguments}
    converted_arguments = {
        name: _convert_value(
            value,
            entry.properties[name],
            f"{path}.{factory_name}.{name}",
            type_modules,
            contract_root,
            factory_context,
        )
        for name, value in values.items()
        if value is not None or name not in entry.optional or entry.optional[name] is not None
    }

    factory = getattr(factory_module, factory_name)
    qualified_name = f"{factory_module.__name__}.{factory_name}"
    context_arguments = (factory_context or {}).get(
        qualified_name, (factory_context or {}).get(factory_name, {})
    )
    overlap = converted_arguments.keys() & context_arguments.keys()
    if overlap:
        raise JSONSettingsValidationError(
            f"Context duplicates JSON properties for {qualified_name}: {sorted(overlap)}"
        )
    try:
        return factory(**converted_arguments, **context_arguments)
    except (TypeError, ValueError) as error:
        raise JSONSettingsValidationError(
            f"Tudat rejected {path}.{factory_name}: {error}"
        ) from error


def module_name_from_contract(contract_path: str | Path) -> str:
    """Derive a Tudatpy module name from a contract's location."""

    relative = Path(contract_path).resolve().relative_to(CONTRACT_ROOT.resolve())
    parts = relative.with_suffix("").parts
    prefix = (
        "tudatpy.dynamics" if parts[0] in {"environment_setup", "propagation_setup"} else "tudatpy"
    )
    for length in range(len(parts), 0, -1):
        candidate = prefix + "." + ".".join(parts[:length])
        try:
            import_module(candidate)
        except ModuleNotFoundError as error:
            if error.name != candidate:
                raise
        else:
            return candidate
    raise JSONSettingsValidationError(f"No Tudatpy module corresponds to {contract_path}")


def _convert_settings_tree(
    value: Any,
    contract: dict[str, ContractEntry],
    factory_module: ModuleType,
    type_modules: tuple[ModuleType, ...],
    path: str,
    factory_expected: bool = False,
    factory_context: Mapping[str, Mapping[str, Any]] | None = None,
) -> Any:
    if isinstance(value, list):
        return [
            _convert_settings_tree(
                item,
                contract,
                factory_module,
                type_modules,
                f"{path}[{index}]",
                True,
                factory_context,
            )
            for index, item in enumerate(value)
        ]
    if not isinstance(value, dict):
        return value
    if factory_expected:
        return create_settings_object(
            value,
            contract,
            factory_module,
            type_modules,
            path,
            factory_context=factory_context,
        )
    if len(value) == 1 and next(iter(value)) in contract:
        return create_settings_object(
            value,
            contract,
            factory_module,
            type_modules,
            path,
            factory_context=factory_context,
        )
    return {
        key: _convert_settings_tree(
            item,
            contract,
            factory_module,
            type_modules,
            f"{path}.{key}",
            factory_context=factory_context,
        )
        for key, item in value.items()
    }


def load_settings(
    settings_path: str | Path,
    contract_path: str | Path,
    module_name: str | None = None,
    factory_context: Mapping[str, Mapping[str, Any]] | None = None,
) -> Any:
    """Load settings recursively using an arbitrary contract."""

    from tudatpy import dynamics

    factory_module = import_module(module_name or module_name_from_contract(contract_path))
    type_modules = (dynamics,)
    contract = load_contract(contract_path, factory_module, type_modules)
    document = read_json_object(settings_path, "settings")
    return _convert_settings_tree(
        document,
        contract,
        factory_module,
        type_modules,
        "settings",
        factory_context=factory_context,
    )
