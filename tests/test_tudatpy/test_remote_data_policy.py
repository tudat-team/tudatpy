import importlib.util
import socket
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError, URLError

import pytest
import requests

_CONFTEST_SPEC = importlib.util.spec_from_file_location(
    "remote_data_conftest", Path(__file__).with_name("conftest.py")
)
conftest = importlib.util.module_from_spec(_CONFTEST_SPEC)
_CONFTEST_SPEC.loader.exec_module(conftest)


class _HookOutcome:
    def __init__(self, report):
        self.report = report

    def get_result(self):
        return self.report


def _run_report_hook(exception, remote_data=True):
    try:
        raise exception
    except Exception:
        excinfo = pytest.ExceptionInfo.from_current()

    item = SimpleNamespace(
        path="test_remote.py",
        get_closest_marker=lambda marker: (
            object() if marker == "remote_data" and remote_data else None
        ),
    )
    call = SimpleNamespace(excinfo=excinfo)
    report = SimpleNamespace(when="call", outcome="failed", longrepr=None)
    hook = conftest.pytest_runtest_makereport(item, call)
    next(hook)
    with pytest.raises(StopIteration):
        hook.send(_HookOutcome(report))
    return report


@pytest.mark.parametrize(
    "exception",
    (
        requests.ConnectionError("offline"),
        requests.ConnectTimeout("connect timed out"),
        requests.ReadTimeout("read timed out"),
        ConnectionRefusedError("connection refused"),
        ConnectionResetError("connection reset"),
        socket.gaierror("name resolution failed"),
        URLError("offline"),
    ),
)
def test_remote_connectivity_failure_is_skipped(exception):
    report = _run_report_hook(exception)

    assert report.outcome == "skipped"
    assert "Remote service unavailable:" in report.longrepr[2]


def _http_error(status_code, client):
    if client == "urllib":
        return HTTPError("https://example.invalid", status_code, "HTTP error", None, None)
    response = requests.Response()
    response.status_code = status_code
    return requests.HTTPError(f"HTTP {status_code}", response=response)


@pytest.mark.parametrize("client", ("requests", "urllib"))
@pytest.mark.parametrize("status_code", (502, 503, 504))
@pytest.mark.parametrize("remote_data, expected_outcome", ((True, "skipped"), (False, "failed")))
def test_http_service_outage_requires_remote_data_marker(
    client, status_code, remote_data, expected_outcome
):
    report = _run_report_hook(_http_error(status_code, client), remote_data=remote_data)

    assert report.outcome == expected_outcome
    if remote_data:
        assert "Remote service unavailable:" in report.longrepr[2]


@pytest.mark.parametrize("client", ("requests", "urllib"))
@pytest.mark.parametrize("status_code", (400, 401, 403, 404, 429, 500, 501))
def test_other_http_errors_still_fail(client, status_code):
    report = _run_report_hook(_http_error(status_code, client))

    assert report.outcome == "failed"


@pytest.mark.parametrize(
    "outage_message",
    (
        "wldini(): missing required file LTKERNL",
        "ERROR in VLRDC: Var not declared: IP_ADDR",
    ),
)
def test_jpl_horizons_backend_outage_is_skipped(outage_message):
    report = _run_report_hook(
        ValueError(
            "Query failed without known error message; received the following "
            "response:\nAPI SOURCE: NASA/JPL Horizons API\n" + outage_message
        )
    )

    assert report.outcome == "skipped"
    assert "Remote service unavailable:" in report.longrepr[2]


@pytest.mark.parametrize(
    "exception",
    (
        AssertionError("incorrect result"),
        ValueError("invalid test data"),
        requests.HTTPError("HTTP 500"),
        requests.exceptions.InvalidURL("invalid URL"),
        requests.exceptions.InvalidJSONError("invalid JSON"),
        requests.exceptions.SSLError("certificate verification failed"),
        HTTPError("https://example.invalid", 500, "server error", None, None),
    ),
)
def test_non_connectivity_failure_still_fails(exception):
    report = _run_report_hook(exception)

    assert report.outcome == "failed"
