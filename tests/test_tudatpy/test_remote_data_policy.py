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


def _run_report_hook(exception):
    try:
        raise exception
    except Exception:
        excinfo = pytest.ExceptionInfo.from_current()

    item = SimpleNamespace(
        path="test_remote.py",
        get_closest_marker=lambda marker: object() if marker == "remote_data" else None,
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


@pytest.mark.parametrize(
    "exception",
    (
        AssertionError("incorrect result"),
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
