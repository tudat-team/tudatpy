from tudatpy.data_input import data_retrieval
from tudatpy.data_input.data_retrieval import missions
from tudatpy.data_input.data_retrieval.missions import LoadPDS


def test_mission_downloader_public_imports():
    assert data_retrieval.missions is missions
    assert missions.LoadPDS is LoadPDS


def test_load_pds_includes_mission_mixins():
    loader = LoadPDS()

    for method_name in (
        "get_mission_files",
        "get_cassini_flyby_files",
        "get_grail_a_files",
        "get_grail_b_files",
        "get_juice_files",
        "get_mex_files",
        "get_mro_files",
        "get_ro_files",
    ):
        assert callable(getattr(loader, method_name))
