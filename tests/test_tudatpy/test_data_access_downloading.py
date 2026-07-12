from tudatpy.data_access import downloading
from tudatpy.data_access.downloading import missions
from tudatpy.data_access.downloading.missions import LoadPDS


def test_mission_downloader_public_imports():
    assert downloading.missions is missions
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
