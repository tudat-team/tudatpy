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
    """Report invalid contracts, JSON settings, references, or factory input.

    This single public exception type gives callers one error boundary for
    failures detected by the JSON layer. The original Tudat ``TypeError`` or
    ``ValueError`` is retained as the exception cause when a factory rejects
    arguments after contract validation.
    """


@dataclass(frozen=True)
class ContractEntry:
    """Store the validated contract for one settings factory.

    Attributes
    ----------
    properties : dict[str, str]
        Mapping from exposed factory argument names to contract type
        expressions.
    optional : dict[str, Any]
        Mapping of optional property names to their JSON-representable defaults.
        A ``None`` value delegates to the Tudat factory default.
    unsupported : frozenset[str]
        Declared factory inputs that remain part of the shared data model but
        cannot currently be supplied through the JSON interface.
    """

    properties: dict[str, str]
    optional: dict[str, Any]
    unsupported: frozenset[str]


_TYPE_PATTERN = re.compile(
    r"^(?P<base>[A-Za-z_][A-Za-z0-9_.]*)(?P<containers>(?:(?:\[(?:\d+)?\])|\{\})*)$"
)
_CONTAINER_PATTERN = re.compile(r"\[(\d*)\]|(\{\})")
_BUILTIN_TYPES = {"any", "bool", "complex", "float", "int", "null", "object", "string"}
CONTRACT_ROOT = Path(__file__).parent / "contracts"
DEFAULT_TYPE_MODULE_NAMES = ("tudatpy.dynamics", "tudatpy.math", "tudatpy.astro")


def _reject_non_finite(value: str) -> None:
    """Reject a non-standard, non-finite JSON numeric token.

    Parameters
    ----------
    value : str
        Token encountered by ``json.load``; normally ``NaN``, ``Infinity``, or
        ``-Infinity``.

    Returns
    -------
    None
        This callback never returns normally.

    Raises
    ------
    JSONSettingsValidationError
        Always, because contracts and settings require portable finite JSON.
    """

    raise JSONSettingsValidationError(f"Non-finite JSON number {value!r} is not permitted")


def _read_json_value(path: Path) -> Any:
    """Read one JSON document without resolving its references.

    Parameters
    ----------
    path : pathlib.Path
        File to decode as UTF-8 JSON.

    Returns
    -------
    Any
        Decoded JSON value. The root is not required to be an object here
        because referenced documents may contain arrays or scalar values.

    Raises
    ------
    JSONSettingsValidationError
        If the file cannot be read, contains malformed JSON, or contains a
        non-finite numeric extension.
    """

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
    """Recursively replace relative ``$ref`` objects.

    Parameters
    ----------
    value : Any
        Decoded JSON subtree currently being inspected.
    source : pathlib.Path
        File containing ``value``. Relative references are resolved against its
        parent directory.
    active : frozenset[pathlib.Path]
        Files in the current reference chain, used to detect direct and indirect
        circular references.

    Returns
    -------
    Any
        Equivalent JSON value with every nested reference replaced by the fully
        resolved contents of its target file.

    Raises
    ------
    JSONSettingsValidationError
        If a reference has siblings, is not a string, points to invalid JSON, or
        creates a reference cycle.
    """

    if isinstance(value, list):
        return [_resolve_references(item, source, active) for item in value]
    if not isinstance(value, dict):
        return value
    if "$ref" in value:
        # A reference replaces its complete object. Sibling keys would make the
        # merge semantics ambiguous, so they are deliberately forbidden.
        if set(value) != {"$ref"} or not isinstance(value["$ref"], str):
            raise JSONSettingsValidationError(
                f"A JSON reference in {source} must contain only a string $ref"
            )
        referenced = (source.parent / value["$ref"]).resolve()
        if referenced in active:
            raise JSONSettingsValidationError(f"Circular JSON reference to {referenced}")
        return _resolve_references(_read_json_value(referenced), referenced, active | {referenced})
    return {key: _resolve_references(item, source, active) for key, item in value.items()}


def read_json_document(path: str | Path) -> Any:
    """Read any JSON value and resolve every relative file reference.

    Parameters
    ----------
    path : str or pathlib.Path
        Root JSON document to read.

    Returns
    -------
    Any
        Fully resolved JSON value. The root may be an object, array, string,
        number, boolean, or null.

    Raises
    ------
    JSONSettingsValidationError
        If reading or reference resolution fails.
    """

    path = Path(path).resolve()
    return _resolve_references(_read_json_value(path), path, frozenset({path}))


def read_json_object(path: str | Path, document_name: str) -> dict[str, Any]:
    """Read a JSON object and resolve every relative file reference.

    Parameters
    ----------
    path : str or pathlib.Path
        Root JSON document to read.
    document_name : str
        Human-readable name inserted into validation errors, such as
        ``"contract"`` or ``"settings"``.

    Returns
    -------
    dict[str, Any]
        Fully resolved root object.

    Raises
    ------
    JSONSettingsValidationError
        If reading or reference resolution fails, or the resolved root value is
        not a JSON object.
    """

    return expect_object(read_json_document(path), document_name)


def expect_object(value: Any, path: str) -> dict[str, Any]:
    """Require a value to be a JSON object.

    Parameters
    ----------
    value : Any
        Decoded JSON value to check.
    path : str
        Logical location included in a validation error.

    Returns
    -------
    dict[str, Any]
        The original value, narrowed to an object.

    Raises
    ------
    JSONSettingsValidationError
        If ``value`` is not a dictionary.
    """

    if not isinstance(value, dict):
        raise JSONSettingsValidationError(f"{path} must be an object, got {type(value).__name__}")
    return value


def _parse_type(type_name: Any, path: str) -> tuple[str, list[int | None | str]]:
    """Split a contract type expression into a base type and containers.

    Parameters
    ----------
    type_name : Any
        Contract value expected to contain an expression such as ``"float[3]"``,
        ``"string[]"``, or ``"environment_setup.shape{}"``.
    path : str
        Logical contract location used in validation errors.

    Returns
    -------
    tuple[str, list[int or None or str]]
        Base type name and ordered container descriptors. A positive integer is
        a fixed-size array, ``None`` is a variable-size array, and ``"dict"`` is
        a JSON object whose values have the remaining inner type.

    Raises
    ------
    JSONSettingsValidationError
        If the value is not a string or uses unsupported type-expression syntax.
    """

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
    """Find one unambiguous enum or named type below configured modules.

    Parameters
    ----------
    type_name : str
        Unqualified public type name recorded in a contract.
    type_modules : tuple[types.ModuleType, ...]
        Root modules whose Tudatpy submodules and re-exports are searched.

    Returns
    -------
    Any
        The unique exposed Python type object with the requested name.

    Raises
    ------
    JSONSettingsValidationError
        If no configured module exposes the name or multiple distinct type
        objects expose it.
    """

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

    # Public Tudat modules often re-export the same enum object. Identity-based
    # de-duplication accepts those aliases while still detecting distinct types
    # that happen to share a name.
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
    """Validate and convert one scalar or nested settings-factory value.

    Parameters
    ----------
    value : Any
        Decoded JSON value at the innermost level of its container expression.
    type_name : str
        Built-in type, enum name, or dotted nested-contract identifier.
    path : str
        Logical settings path used in validation errors.
    type_modules : tuple[types.ModuleType, ...]
        Modules searched when resolving enum names.
    contract_root : pathlib.Path
        Root directory used to translate dotted nested types into contract files.
    factory_context : Mapping[str, Mapping[str, Any]] or None
        Runtime-only keyword arguments available to nested factory calls.

    Returns
    -------
    Any
        Validated Python scalar, enum member, complex number, or Tudat settings
        object produced from a nested contract.

    Raises
    ------
    JSONSettingsValidationError
        If the value has the wrong JSON type, is non-finite, names an unknown
        enum member, or fails nested settings construction.
    """

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
        # A dotted type identifies another contract, not a Python class name.
        # Its value is therefore another one-key settings factory definition.
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
    """Convert a JSON value using a complete contract type expression.

    Parameters
    ----------
    value : Any
        Decoded JSON value to validate and convert. A bare value is accepted as
        shorthand for a variable-length array containing one item.
    type_expression : Any
        Contract expression, optionally containing ``|`` alternatives and array
        or object containers.
    path : str
        Logical settings path used in validation errors.
    type_modules : tuple[types.ModuleType, ...]
        Modules searched for enum and named-type definitions.
    contract_root : pathlib.Path, optional
        Root of nested contract files.
    factory_context : Mapping[str, Mapping[str, Any]] or None, optional
        Runtime-only keyword arguments made available to nested factories.

    Returns
    -------
    Any
        Converted scalar, container, enum member, or settings object matching the
        first successful alternative in the type expression.

    Raises
    ------
    JSONSettingsValidationError
        If no alternative accepts the value, a container has the wrong kind or
        length, or scalar conversion fails.
    """

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
        """Apply parsed containers recursively from outermost to innermost.

        Parameters
        ----------
        current_value : Any
            JSON subtree at the current recursion level.
        remaining_dimensions : list[int or None or str]
            Containers not yet consumed.
        current_path : str
            Logical location of ``current_value`` for precise errors.

        Returns
        -------
        Any
            Converted subtree, or the converted scalar after the final container
            has been consumed.

        Raises
        ------
        JSONSettingsValidationError
            If a container kind or fixed length does not match the contract.
        """

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
            # A variable-length array containing one value may omit its square
            # brackets. Fixed-size arrays retain their explicit JSON shape so
            # vectors and matrices remain unambiguous.
            if expected_size is None:
                current_value = [current_value]
            else:
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
    """Load and structurally validate a contract for a Tudat factory module.

    Parameters
    ----------
    contract_path : str or pathlib.Path
        JSON contract containing factory names and their property definitions.
    factory_module : types.ModuleType
        Imported module that must expose every factory named by the contract.
    type_modules : tuple[types.ModuleType, ...]
        Root modules searched for enum types used by supported properties.

    Returns
    -------
    dict[str, ContractEntry]
        Validated entries keyed by factory name.

    Raises
    ------
    JSONSettingsValidationError
        If the contract structure, type expressions, defaults, unsupported lists,
        or factory names are invalid.
    """

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
                # Unsupported inputs retain their intended type in the shared
                # contract even when that type is outside the modules currently
                # understood by this JSON layer.
                if (
                    property_name not in unsupported
                    and base_type not in _BUILTIN_TYPES
                    and "." not in base_type
                ):
                    _find_named_type(base_type, type_modules)
            typed_properties[property_name] = type_expression

        for property_name, default_value in optional.items():
            # JSON null represents defaults that are None or not finite in
            # Python. Unsupported optional inputs use the Tudat default.
            if property_name not in unsupported and default_value is not None:
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
    """Read argument names from the first exposed factory overload.

    Parameters
    ----------
    factory : Any
        Python or pybind callable whose named arguments are inspected.
    factory_name : str
        Public name used to locate a generated signature in pybind documentation.

    Returns
    -------
    tuple[set[str], set[str]]
        Set of all named arguments and the subset that expose defaults.

    Raises
    ------
    JSONSettingsValidationError
        If the callable has no finite named-parameter interface or its generated
        signature cannot be parsed.
    """

    try:
        parameters = signature(factory).parameters.values()
    except (TypeError, ValueError):
        # Pybind functions do not always expose an inspectable signature. Their
        # generated first signature is then parsed from the docstring instead.
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
    """Compare an arbitrary JSON contract with an importable Tudatpy module.

    Parameters
    ----------
    contract_path : str or pathlib.Path
        Contract to validate.
    factory_module_name : str or None, optional
        Fully qualified module exposing the contracted factories. When omitted,
        the module is derived from the contract's location below the packaged
        contract root.
    type_module_names : collections.abc.Sequence[str], optional
        Fully qualified root modules searched for contracted enum types. By
        default, the public dynamics, mathematics, and astrodynamics modules are
        searched.

    Returns
    -------
    int
        Number of factory entries successfully validated.

    Raises
    ------
    JSONSettingsValidationError
        If contract properties or optional markers differ from the first exposed
        factory signature, or structural contract validation fails.
    ModuleNotFoundError
        If a requested module cannot be imported.
    """

    factory_module_name = factory_module_name or module_name_from_contract(contract_path)
    factory_module = import_module(factory_module_name)
    default_type_modules = DEFAULT_TYPE_MODULE_NAMES
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
    """Validate every JSON contract found below a contract root.

    Parameters
    ----------
    contract_root : str or pathlib.Path, optional
        Directory searched recursively for ``*.json`` contracts. The packaged
        contract directory is used by default.

    Returns
    -------
    int
        Total number of validated factory entries across all discovered files.

    Raises
    ------
    JSONSettingsValidationError
        If any contract is invalid or differs from its corresponding Tudat API.
    ModuleNotFoundError
        If a corresponding factory module cannot be imported.
    """

    total = 0
    for contract_path in sorted(Path(contract_root).rglob("*.json")):
        total += validate_contract_against_module(contract_path)
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
    """Validate one settings definition and invoke its contracted Tudat factory.

    Parameters
    ----------
    definition : Any
        JSON object containing exactly one factory name mapped to its arguments,
        or a bare factory-name string for a zero-argument call. A factory mapped
        to JSON null is also treated as a zero-argument call.
    contract : dict[str, ContractEntry]
        Validated contract entries available at this location in the settings
        tree.
    factory_module : types.ModuleType
        Module exposing the selected Tudat factory.
    type_modules : tuple[types.ModuleType, ...]
        Modules searched for enum members during argument conversion.
    path : str
        Logical settings location included in validation errors.
    contract_root : pathlib.Path, optional
        Root directory used to locate contracts for nested settings properties.
    factory_context : Mapping[str, Mapping[str, Any]] or None, optional
        Runtime-only keyword arguments keyed by qualified or unqualified factory
        name. Context arguments cannot also appear as JSON properties.

    Returns
    -------
    Any
        Tudat settings object returned by the selected factory. Generic test
        factories may return another Python value.

    Raises
    ------
    JSONSettingsValidationError
        If the definition is malformed, names an unknown factory, omits required
        properties, supplies undeclared or unsupported properties, conflicts with
        runtime context, fails conversion, or is rejected by Tudat.
    """

    # A bare factory name is shorthand for a zero-argument invocation.
    if isinstance(definition, str):
        definition = {definition: {}}
    definition = expect_object(definition, path)
    if len(definition) != 1:
        raise JSONSettingsValidationError(f"{path} must contain exactly one settings factory")

    factory_name, arguments = next(iter(definition.items()))
    if factory_name not in contract:
        raise JSONSettingsValidationError(
            f"{path} contains unknown settings factory {factory_name!r}"
        )

    # JSON null is the valid explicit no-arguments spelling. A dangling colon
    # without a value is not JSON and is rejected by the JSON parser itself.
    if arguments is None:
        arguments = {}
    arguments = expect_object(arguments, f"{path}.{factory_name}")
    entry = contract[factory_name]

    # A required unsupported property makes the factory unavailable. An
    # optional one is accepted only when omitted, so Tudat uses its default.
    supplied_unsupported = arguments.keys() & entry.unsupported
    required_unsupported = entry.unsupported - entry.optional.keys()
    unsupported = supplied_unsupported | required_unsupported
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

    # Contract defaults are operational: supplied values override them before
    # conversion. Null and unsupported defaults are omitted from the call and
    # consequently delegate to the Tudat factory default.
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
        if name not in entry.unsupported
        and (value is not None or name not in entry.optional or entry.optional[name] is not None)
    }

    factory = getattr(factory_module, factory_name)
    qualified_name = f"{factory_module.__name__}.{factory_name}"
    # Runtime-only objects, such as SystemOfBodies, are injected separately and
    # cannot be overridden by an identically named JSON property.
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
    """Derive the owning Tudatpy factory module from a packaged contract path.

    Parameters
    ----------
    contract_path : str or pathlib.Path
        Contract located below ``CONTRACT_ROOT``. Its directory components mirror
        Tudatpy modules; a final settings-family filename need not be a module.

    Returns
    -------
    str
        Longest importable fully qualified module name represented by the path.

    Raises
    ------
    ValueError
        If the contract is outside the packaged contract root.
    JSONSettingsValidationError
        If no path prefix corresponds to an importable Tudatpy module.
    ModuleNotFoundError
        If importing a candidate fails because one of its own dependencies is
        missing rather than because the candidate module does not exist.
    """

    relative = Path(contract_path).resolve().relative_to(CONTRACT_ROOT.resolve())
    parts = relative.with_suffix("").parts
    prefix = (
        "tudatpy.dynamics" if parts[0] in {"environment_setup", "propagation_setup"} else "tudatpy"
    )
    # Contract filenames may name a settings family below the module that owns
    # its factories. Choose the longest importable module prefix.
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
    """Walk an outer settings tree and convert all recognized factory leaves.

    Parameters
    ----------
    value : Any
        Current JSON subtree. Structural dictionaries are preserved while lists
        and one-key factory objects are converted recursively. A recognized bare
        factory-name string represents a zero-argument factory call.
    contract : dict[str, ContractEntry]
        Factory entries recognized at leaves of this tree.
    factory_module : types.ModuleType
        Module exposing the contracted factories.
    type_modules : tuple[types.ModuleType, ...]
        Modules searched for enum values during conversion.
    path : str
        Logical location of ``value`` for validation errors.
    factory_expected : bool, optional
        Whether ``value`` must be interpreted as a factory definition. List items
        set this flag because settings lists contain factory objects.
    factory_context : Mapping[str, Mapping[str, Any]] or None, optional
        Runtime-only keyword arguments available to factory calls.

    Returns
    -------
    Any
        Tree with structural mappings retained and every contracted factory leaf
        replaced by its returned Tudat settings object.

    Raises
    ------
    JSONSettingsValidationError
        If a required factory definition or any nested value is invalid.
    """

    if isinstance(value, list):
        # Lists in settings trees contain factory definitions, such as the
        # accelerations exerted by one body or dependent variables to save.
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
    if factory_expected or isinstance(value, str):
        return create_settings_object(
            value,
            contract,
            factory_module,
            type_modules,
            path,
            factory_context=factory_context,
        )
    if not isinstance(value, dict):
        return value
    if len(value) == 1 and next(iter(value)) in contract:
        # A recognized one-key object is a factory invocation. Other mappings
        # are structural dictionaries whose values are searched recursively.
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
    """Resolve a JSON document and construct settings using an arbitrary contract.

    Parameters
    ----------
    settings_path : str or pathlib.Path
        Root settings JSON document. Relative ``$ref`` objects are resolved before
        contract conversion begins.
    contract_path : str or pathlib.Path
        Contract defining the permitted factory names, arguments, types, defaults,
        and unsupported inputs.
    module_name : str or None, optional
        Fully qualified module exposing the contracted factories. When omitted,
        it is derived from ``contract_path``.
    factory_context : Mapping[str, Mapping[str, Any]] or None, optional
        Runtime-only keyword arguments to inject into selected factories.

    Returns
    -------
    Any
        Converted settings object or structural tree of settings objects matching
        the outer shape of the input document.

    Raises
    ------
    JSONSettingsValidationError
        If reading, reference resolution, contract loading, value conversion, or
        Tudat factory invocation fails.
    ModuleNotFoundError
        If the selected factory module cannot be imported.
    """

    factory_module = import_module(module_name or module_name_from_contract(contract_path))
    type_modules = tuple(import_module(name) for name in DEFAULT_TYPE_MODULE_NAMES)
    contract = load_contract(contract_path, factory_module, type_modules)
    # Unlike contracts, settings documents need not have an object root. A
    # bare factory-name string is the universal shorthand for invoking a
    # factory without explicitly supplied arguments.
    document = read_json_document(settings_path)
    return _convert_settings_tree(
        document,
        contract,
        factory_module,
        type_modules,
        "settings",
        factory_context=factory_context,
    )
