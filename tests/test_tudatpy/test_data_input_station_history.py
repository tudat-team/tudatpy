"""Binding/import smoke only; parser and scientific assertions are in C++."""

from pathlib import Path
from tudatpy.data_input.environment_data import station_history

FIXTURES = Path(__file__).parents[1] / "test_tudat/reference_files/station_history"


def test_station_history_namespace_keys_and_array_binding():
    history = station_history.read_sinex_station_history(str(FIXTURES / "ab07_solutions.snx"))
    key = station_history.SinexSolutionKey("AB07", "A", "1")
    solution = station_history.get_station_coordinate_solution(history, key)
    assert isinstance(history.sites, list)
    assert isinstance(history.parameters, dict)
    assert len(solution.position_metres) == 3
    assert key.site.site_code == "AB07"
    duplicate_history = station_history.read_sinex_station_history(
        str(FIXTURES / "duplicate_sites.snx")
    )
    assert isinstance(duplicate_history.sites, list)
    assert len([site for site in duplicate_history.sites if site.key.site_code == "GALA"]) == 2
    duplicate_history.release = "explicit source release"
    assert duplicate_history.release == "explicit source release"
    station_history.validate_sinex_station_history(duplicate_history)


def test_unresolved_optional_binding_is_none():
    history = station_history.read_sinex_station_history(str(FIXTURES / "ab07_psd.snx"))
    relaxation = next(row for row in history.parameters.values() if row.type.startswith("TEXP"))
    assert relaxation.value_si is None
    assert relaxation.standard_deviation_si is None
