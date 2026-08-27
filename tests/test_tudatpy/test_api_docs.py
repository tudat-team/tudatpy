import importlib.util
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def api_docs_conf():
    conf_path = Path(__file__).resolve().parents[2] / "docs" / "tudatpy" / "source" / "conf.py"
    spec = importlib.util.spec_from_file_location("tudatpy_api_docs_conf", conf_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("annotation", "expected"),
    [
        ("typing.SupportsInt | typing.SupportsIndex", "int"),
        ("int | SupportsIndex", "int"),
        ("typing.SupportsFloat | typing.SupportsIndex", "float"),
        ("float | SupportsIndex", "float"),
        (
            "typing.SupportsComplex | typing.SupportsFloat | typing.SupportsIndex",
            "complex",
        ),
        ("SupportsComplex | float | SupportsIndex", "complex"),
        (
            "collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]",
            "list[float]",
        ),
        (
            "collections.abc.Mapping[str, typing.SupportsInt | typing.SupportsIndex]",
            "dict[str, int]",
        ),
        (
            'typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"]',
            "numpy.ndarray[numpy.float64[3, 1]]",
        ),
        (
            "typing.Annotated[numpy.typing.ArrayLike, numpy.float64]",
            "numpy.ndarray[numpy.float64]",
        ),
        (
            "typing.Annotated[numpy.typing.NDArray[numpy.float64], '[m, n]']",
            "numpy.ndarray[numpy.float64[m, n]]",
        ),
        (
            "~tudatpy.kernel.dynamics.environment.SystemOfBodies",
            "tudatpy.kernel.dynamics.environment.SystemOfBodies",
        ),
    ],
)
def test_simplify_type_annotations(api_docs_conf, annotation, expected):
    assert api_docs_conf.simplify_type_annotations(annotation) == expected


def test_overload_signatures_are_simplified(api_docs_conf):
    lines = [
        "Overloaded function.",
        "",
        "1. select(values: collections.abc.Sequence[typing.SupportsFloat | "
        "typing.SupportsIndex]) -> typing.SupportsInt | typing.SupportsIndex",
    ]

    assert api_docs_conf._rewrite_overloaded_list_items(lines) == [
        "Overloaded function.",
        "",
        "Overload 1:",
        "``select(values: list[float]) -> int``",
    ]


def test_readonly_property_is_marked_once(api_docs_conf):
    readonly_property = property(lambda instance: instance.value)
    lines = ["Property description."]

    api_docs_conf.mark_property_mutability(
        None, "property", "Example.value", readonly_property, None, lines
    )
    api_docs_conf.mark_property_mutability(
        None, "property", "Example.value", readonly_property, None, lines
    )

    assert lines == ["**read-only**", "", "Property description."]


def test_writable_property_is_not_marked_readonly(api_docs_conf):
    writable_property = property(
        lambda instance: instance.value,
        lambda instance, value: setattr(instance, "value", value),
    )
    lines = ["**read-only**", "", "Property description."]

    api_docs_conf.mark_property_mutability(
        None, "property", "Example.value", writable_property, None, lines
    )

    assert "**read-only**" not in lines
