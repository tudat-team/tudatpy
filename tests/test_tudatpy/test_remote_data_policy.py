import importlib.util
from pathlib import Path
from types import SimpleNamespace

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


def test_remote_connection_failure_is_skipped():
    report = _run_report_hook(requests.ConnectionError("offline"))

    assert report.outcome == "skipped"
    assert "Remote service unavailable: offline" in report.longrepr[2]


def test_remote_assertion_failure_still_fails():
    report = _run_report_hook(AssertionError("incorrect result"))

    assert report.outcome == "failed"
