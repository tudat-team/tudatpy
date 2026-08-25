import os

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--remote-data",
        action="store_true",
        dest="remote_data",
        default=False,
        help="run tests that require a network connection (disabled by default)",
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
                pytest.mark.skip(reason="remote-data tests disabled; use --remote-data to enable")
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
