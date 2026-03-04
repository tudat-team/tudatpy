import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--remote-data", action="store_true", default=False, help="run tests that require network connection"
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "remote_data: mark test as requiring internet access")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--remote-data"):
        # --remote-data given: do not skip remote tests
        return
    skip_remote = pytest.mark.skip(reason="need --remote-data option to run")
    for item in items:
        if "remote_data" in item.keywords:
            item.add_marker(skip_remote)