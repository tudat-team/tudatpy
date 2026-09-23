import os
import socket
from urllib.error import HTTPError, URLError

import pytest
import requests

_REMOTE_SERVICE_UNAVAILABLE_STATUSES = (502, 503, 504)
_JPL_HORIZONS_OUTAGE_SIGNATURES = (
    "wldini(): missing required file LTKERNL",
    "ERROR in VLRDC: Var not declared: IP_ADDR",
)


def _implicit_decimal_field(value, integer_width, fraction_width=4, signed=False):
    """Return a fixed-width MPC field with an implied decimal point."""
    if value is None:
        return " " * (integer_width + fraction_width)
    width = integer_width + fraction_width + 1
    return f"{value:{'+' if signed else ''}0{width}.{fraction_width}f}".replace(".", "")


def _mpc_radar_pair(
    number="00433",
    date="1990 07 15.326389",
    delay_us=None,
    delay_sigma_us=None,
    doppler_hz=None,
    doppler_sigma_hz=None,
    frequency_mhz=2380.0,
    transmitter="251",
    receiver="251",
    bounce_point="C",
):
    """Return two 80-column MPC radar records."""
    head = f"{number:<5}{'':7}  "
    tail = f"{transmitter:>3}{'':6}{receiver:>3}"
    first = (
        head
        + "R"
        + date
        + _implicit_decimal_field(delay_us, 11)
        + _implicit_decimal_field(doppler_hz, 11, signed=True)
        + _implicit_decimal_field(frequency_mhz, 5, 1)
        + tail
    )
    second = (
        head
        + "r"
        + date
        + bounce_point
        + _implicit_decimal_field(delay_sigma_us, 10)
        + _implicit_decimal_field(doppler_sigma_hz, 11)
        + " " * 6
        + tail
    )
    assert len(first) == len(second) == 80
    return [first, second]


@pytest.fixture
def mpc_radar_pair():
    """Return a factory for MPC radar record pairs."""
    return _mpc_radar_pair


def _is_connectivity_failure(exception):
    """Return whether an exception represents an unavailable remote service."""
    # Gateways can respond even when the upstream service is unavailable.
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
        "remote_data(required_env=(), service='remote-data'): mark a test as requiring "
        "internet access and, optionally, environment variables containing credentials",
    )


def pytest_collection_modifyitems(config, items):
    run_remote_tests = config.getoption("remote_data")
    missing_credentials_by_service = {}

    for item in items:
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
