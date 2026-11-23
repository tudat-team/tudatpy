"""Exercise the ordinary Python factory and native solver, without Earth resources."""

import importlib.util
import os
from pathlib import Path
import pickle

import numpy as np
import pytest

import tudatpy
from tudatpy.dynamics import propagation_setup
from tudatpy.interface import spice

acceleration = propagation_setup.acceleration


def test_backend_information():
    available = tudatpy.is_openmp_available()
    assert type(available) is bool
    expected = os.environ.get("TUDAT_EXPECT_OPENMP")
    if expected is not None:
        assert expected in {"0", "1"}
        assert available == (expected == "1")


def test_sequential_default_and_existing_positional_arguments():
    assert acceleration.spherical_harmonic_gravity(20, 20).number_of_threads == 1
    assert acceleration.spherical_harmonic_gravity(20, 20, True).number_of_threads == 1
    assert (
        acceleration.spherical_harmonic_gravity(
            maximum_degree=20, maximum_order=20, number_of_threads=1
        ).number_of_threads
        == 1
    )


@pytest.mark.parametrize("threads", [0, -1, -20])
def test_reject_invalid_counts(threads):
    with pytest.raises(ValueError, match="positive"):
        acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=threads)


def test_settings_are_local_and_readonly():
    if not tudatpy.is_openmp_available():
        with pytest.raises(RuntimeError, match="TUDAT_BUILD_WITH_OPENMP"):
            acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=2)
        return
    parallel = acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=4)
    sequential = acceleration.spherical_harmonic_gravity(20, 20)
    assert parallel.number_of_threads == 4
    assert sequential.number_of_threads == 1
    with pytest.raises(AttributeError):
        parallel.number_of_threads = 2


@pytest.mark.requires_serialization
def test_pickle_restores_sequential_policy():
    threads = 4 if tudatpy.is_openmp_available() else 1
    settings = acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=threads)
    restored = pickle.loads(pickle.dumps(settings))
    assert type(restored) is type(settings)
    assert restored.number_of_threads == 1
    assert restored == settings


@pytest.fixture(scope="module")
def benchmark():
    path = (
        Path(__file__).resolve().parents[2] / "benchmarks/benchmark_python_spherical_harmonics.py"
    )
    spec = importlib.util.spec_from_file_location("gravity_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    spice.load_standard_kernels()
    return module


@pytest.mark.parametrize("mode", ["dynamics", "variational"])
def test_full_solver_history_parity(benchmark, mode):
    threads = 4 if tudatpy.is_openmp_available() else 1
    reference_case = benchmark.make_case(1, 20, 300.0, mode)
    reference = benchmark.collect(benchmark.execute(reference_case, mode), mode)
    actual_case = benchmark.make_case(threads, 20, 300.0, mode)
    actual = benchmark.collect(benchmark.execute(actual_case, mode), mode)
    maximum = benchmark.compare(reference, actual, 300.0)
    assert np.isfinite(maximum)
    # The native model/cache operation order is preserved in this build.
    assert maximum == 0.0


def test_python_rotation_callback_stays_on_caller_and_remains_alive(benchmark):
    import gc
    import threading
    import weakref
    from tudatpy.dynamics import environment_setup

    threads = 4 if tudatpy.is_openmp_available() else 1
    case = benchmark.make_case(threads, 20, 300.0, "dynamics")
    caller = threading.get_ident()
    calls = []

    class RotationCallback:
        def __call__(self, epoch):
            calls.append(threading.get_ident())
            return np.eye(3)

    callback = RotationCallback()
    reference = weakref.ref(callback)
    rotation = environment_setup.rotation_model.custom_rotation_model(
        "J2000", "IAU_MARS", callback, 1.0
    )
    environment_setup.add_rotation_model(case[0], "Mars", rotation)
    del callback, rotation
    gc.collect()
    assert reference() is not None
    results = benchmark.collect(benchmark.execute(case, "dynamics"), "dynamics")
    benchmark.compare(results, results, 300.0)
    assert calls
    assert set(calls) == {caller}
