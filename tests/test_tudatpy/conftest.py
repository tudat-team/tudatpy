import os
import socket
from urllib.error import HTTPError, URLError

import pytest
import requests

_REMOTE_SERVICE_UNAVAILABLE_STATUSES = (500, 502, 503, 504)
_JPL_HORIZONS_OUTAGE_SIGNATURES = (
    "wldini(): missing required file LTKERNL",
    "ERROR in VLRDC: Var not declared: IP_ADDR",
)


def _is_connectivity_failure(exception):
    """Return whether an exception represents an unavailable remote service."""
    # Servers and gateways respond with these statuses when the service is failing or unavailable.
    if isinstance(exception, requests.exceptions.HTTPError):
        return (
            exception.response is not None
            and exception.response.status_code in _REMOTE_SERVICE_UNAVAILABLE_STATUSES
        )
    if isinstance(exception, HTTPError):
        return exception.code in _REMOTE_SERVICE_UNAVAILABLE_STATUSES

    if isinstance(
        exception,
        (
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            ConnectionRefusedError,
            ConnectionResetError,
            socket.gaierror,
            TimeoutError,
        ),
    ):
        return not isinstance(
            exception,
            (requests.exceptions.SSLError, requests.exceptions.InvalidURL),
        )

    if isinstance(exception, URLError) and not isinstance(exception, HTTPError):
        return True

    # Astroquery raises ValueError when the Horizons API itself is reachable but
    # cannot serve requests because required backend resources are unavailable.
    exception_message = str(exception)
    return isinstance(exception, ValueError) and (
        "API SOURCE: NASA/JPL Horizons API" in exception_message
        and any(signature in exception_message for signature in _JPL_HORIZONS_OUTAGE_SIGNATURES)
    )


def pytest_addoption(parser):
    parser.addoption(
        "--remote-data",
        action="store_true",
        dest="remote_data",
        default=True,
        help="run tests that require a network connection (default)",
    )
    parser.addoption(
        "--no-remote-data",
        action="store_false",
        dest="remote_data",
        help="skip tests that require a network connection",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_serialization: requires a kernel built with TUDAT_BUILD_WITH_SERIALIZATION=ON",
    )
    config.addinivalue_line(
        "markers",
        "remote_data(required_env=(), service='remote-data'): mark a test as requiring "
        "internet access and, optionally, environment variables containing credentials",
    )


def _serialization_enabled(kernel):
    """Use the compiled option; an unknown/old binary must fail, not skip."""
    value = getattr(kernel, "_build_with_serialization", None)
    if type(value) is not bool:
        raise pytest.UsageError(
            "The loaded tudatpy.kernel does not expose a boolean "
            "_build_with_serialization flag. Rebuild the kernel and check its "
            "import path; serialization availability must not be inferred "
            "from missing methods."
        )
    return value


def pytest_collection_modifyitems(config, items):
    serialization_enabled = None
    if any(item.get_closest_marker("requires_serialization") for item in items):
        from tudatpy import kernel

        serialization_enabled = _serialization_enabled(kernel)
    run_remote_tests = config.getoption("remote_data")
    missing_credentials_by_service = {}

    for item in items:
        if item.get_closest_marker("requires_serialization") and not serialization_enabled:
            item.add_marker(
                pytest.mark.skip(reason="kernel built with TUDAT_BUILD_WITH_SERIALIZATION=OFF")
            )
        remote_data_marker = item.get_closest_marker("remote_data")
        if remote_data_marker is None:
            continue

        if not run_remote_tests:
            item.add_marker(
                pytest.mark.skip(reason="remote-data tests disabled with --no-remote-data")
            )
            continue

        required_env = remote_data_marker.kwargs.get("required_env", ())
        missing_env = [variable for variable in required_env if not os.getenv(variable)]
        if not missing_env:
            continue

        item.add_marker(
            pytest.mark.skip(
                reason="missing required environment variable(s): " + ", ".join(missing_env)
            )
        )
        service = remote_data_marker.kwargs.get("service", "Remote-data")
        service_missing_env = missing_credentials_by_service.setdefault(service, [])
        service_missing_env.extend(
            variable for variable in missing_env if variable not in service_missing_env
        )

    terminal_reporter = config.pluginmanager.get_plugin("terminalreporter")
    if terminal_reporter is None:
        return

    for service, missing_env in missing_credentials_by_service.items():
        terminal_reporter.write_line(
            f"NOTICE: {service} remote-data tests will be skipped; "
            "missing environment variable(s): " + ", ".join(missing_env),
            yellow=True,
        )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Treat unavailable remote services as skips, not test failures."""
    outcome = yield
    report = outcome.get_result()

    if (
        report.when not in ("setup", "call")
        or report.outcome != "failed"
        or item.get_closest_marker("remote_data") is None
        or call.excinfo is None
        or not _is_connectivity_failure(call.excinfo.value)
    ):
        return

    report.outcome = "skipped"
    report.longrepr = (
        str(item.path),
        call.excinfo.traceback[-1].lineno,
        f"Remote service unavailable: {call.excinfo.value}",
    )


@pytest.fixture
def serialization_simulator():
    """Fresh short propagation for real result-serialization tests."""

    def create(*, variational=False, hybrid=False):
        from tudatpy.astro import element_conversion
        from tudatpy.astro.time_representation import DateTime
        from tudatpy.data_input.environment_data import spice
        from tudatpy.dynamics import (
            environment_setup,
            propagation_setup,
            parameters_setup,
            simulator,
        )

        spice.load_standard_kernels()
        settings = environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
        settings.add_empty_settings("Delfi-C3")
        settings.get("Delfi-C3").constant_mass = 2.2
        bodies = environment_setup.create_system_of_bodies(settings)
        start = DateTime(2008, 4, 28).to_epoch()
        initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )
        termination = propagation_setup.propagator.time_termination(start + 600.0)
        if hybrid:
            termination = propagation_setup.propagator.hybrid_termination(
                [termination, propagation_setup.propagator.cpu_time_termination(600.0)],
                fulfill_single_condition=True,
            )
        accelerations = propagation_setup.create_acceleration_models(
            bodies,
            {"Delfi-C3": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
            ["Delfi-C3"],
            ["Earth"],
        )
        propagation = propagation_setup.propagator.translational(
            ["Earth"],
            accelerations,
            ["Delfi-C3"],
            initial_state,
            start,
            propagation_setup.integrator.runge_kutta_fixed_step(
                30.0,
                coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4,
            ),
            termination,
            output_variables=[
                propagation_setup.dependent_variable.relative_position("Delfi-C3", "Earth")
            ],
        )
        if variational:
            parameter_settings = parameters_setup.initial_states(propagation, bodies)
            parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
            parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
            return simulator.create_variational_equations_solver(bodies, propagation, parameters)
        return simulator.create_dynamics_simulator(bodies, propagation)

    return create
