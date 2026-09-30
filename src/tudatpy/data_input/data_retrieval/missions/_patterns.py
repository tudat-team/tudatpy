"""Regex filename patterns and static mission reference data for LoadPDS."""

SUPPORTED_PATTERNS = {
    "juice": {
        "ck": r"^(?P<mission>juice)_(?P<instrument>(sc|mga|sa))_(?P<data_type>(meas|crema_\w+_baseline))?_?(?P<reference_file>[a-zA-Z]+)?_?(?P<prod_ID>[a-zA-Z]+)?_?(?P<start_date_file>(\d{6}|\d{8}))?_?(?P<end_date_file>\d{6})?_?(?P<sclk>(s|t|f)\d{6})?_?(?P<version>v\d{2})(?P<extension>\.bc)$",
        "spk": r"^(?P<mission>juice)_(?P<data_type>(orbc|crema))_(?P<reference_file_prod_ID>(\d{6}|\w+_plan))?_?(?P<start_date_file>\d{6})?_?(?P<end_date_file>\d{6})?_?(?P<version>v\d{2})(?P<extension>\.bsp)$",
        "fk": r"^(?P<mission>juice)_(?P<data_type>(dsk_surfaces|events_\w+|ops|roi|sci|stations_topo))_(?P<version>v\d{2})(?P<extension>\.tf)$",
    },
    "mro": {
        "ck": r"^(?P<mission>mro)_(?P<instrument>(sc|hga|sa))_(?P<mission_phase>(cru|ab|psp|rel|psp))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<version>v\d{2})?(?P<extension>\.bc)$",
        "spk": r"^(?P<mission>mro)_(?P<instrument>(sc|hga|sa|struct))_(?P<mission_phase>(cru|ab|psp|rel|psp))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<version>v\d{2})?(?P<extension>\.bsp)$",
        "odf": r"^(?P<mission>mro)(?P<dataset>magr)(?P<date_file>\d{4}_\d{3}_\d{4})(?P<uplink>x|n)(?P<station>nn|mm|[0-9]{2})(?P<downlink>[123m])(?P<version>v\d+)(?P<extension>\.odf)$",
        "tnf": r"^(?P<mission>mro)(?P<dataset>magr)(?P<date_file>\d{4}_\d{3}_\d{4})(?P<uplink>x|n)(?P<station>nn|mm|[0-9]{2})(?P<downlink>[123m])(?P<version>v\d+)(?P<extension>\.tnf)$",
        "tro": r"^(?P<mission>mro)(?P<dataset>magr)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.tro)$",
        "ion": r"^(?P<mission>mro)(?P<dataset>magr)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.ion)$",
    },
    "lro": {
        "ck": r"^(?P<mission>lro)(?P<instrument>(sc|hg|sa|dv|))_(?P<start_date_file>\d{7})_(?P<end_date_file>\d{7})_(?P<version>v\d{2})(?P<extension>\.bc)$"
    },
    "cassini": {
        "ck": r"^(?P<start_date_file>\d{5})_(?P<end_date_file>\d{5})(?P<type>(p|r))(?P<version>[a-z])(?P<freeform>\w+)?(?P<extension>\.bc)$"
    },
    "insight": {
        "ck": r"^(?P<mission>insight)_(?P<purpose>(cru_rec|edl_rec|surf_ops|))_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})_(?P<version>v\d{2})(?P<extension>\.bc)"
    },
    "mex": {
        "ck": r"^(?P<mission>MEX)?_?(?P<data>ATNM)?_?(?P<purpose>(MEASURED|T6|SA))?_?(?P<start_date_file>(\d{4}|\d{6}|P\d{12}))?_?(?P<end_date_file>\d{6})?_?(?P<sclk>S\d{6})?_?(?P<version>(V\d+|\d+))?(?P<extension>\.BC)?$",
        "spk": r"^(?P<data>(ORMM|ORMF))_(?P<SPK_type>T19)_(?P<start_date_file>\d{6})_?(?P<end_date_file>\d{6})_(?P<version>\d{5})(?P<extension>\.BSP)?$",
        "ifms": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)",
        "dp2": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)",
        "dpx": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)",
        "dps": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.tab$)",
    },
    "grail-a": {
        "ck": r"^(?P<mission>gra)_(?P<instrument>rec)_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<extension>\.bc)$",
        "spk": r"^(?P<mission>grail)_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})_(?P<instrument>(sci|nav))_(?P<version>v\d{2})?(?P<extension>\.bsp)$",
        "odf": r"^(?P<mission>(gra))(?P<experiment>lugf)(?P<date_file>\d{4}_\d{3}_\d{4})(?P<pre_version>smmm)(?P<version>v\d{1,2})(?P<extension>\.odf)$",
        "tro": r"^(?P<mission>grx)(?P<experiment>lugf)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.tro)$",
        "ion": r"^(?P<mission>gra)(?P<experiment>lugf)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.ion)$",
    },
    "grail-b": {
        "ck": r"^(?P<mission>grb)_(?P<instrument>rec)_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})(?P<extension>\.bc)$",
        "spk": r"^(?P<mission>grail)_(?P<start_date_file>\d{6})_(?P<end_date_file>\d{6})_(?P<instrument>(sci|nav))_(?P<version>v\d{2})?(?P<extension>\.bsp)$",
        "odf": r"^(?P<mission>(grb))(?P<experiment>lugf)(?P<date_file>\d{4}_\d{3}_\d{4})(?P<pre_version>smmm)(?P<version>v\d{1,2})(?P<extension>\.odf)$",
        "tro": r"^(?P<mission>grx)(?P<experiment>lugf)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.tro)$",
        "ion": r"^(?P<mission>grb)(?P<experiment>lugf)(?P<start_date_file>\d{4}_\d{3})_(?P<end_date_file>\d{4}_\d{3})(?P<extension>\.ion)$",
    },
    "ro": {
        # "DP2": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.TAB$)",
        "dp2": r"^(?P<mission>[a-zA-Z0-9]+)_(?P<band>[a-zA-Z0-9]+)_(?P<date_file>[0-9]{9})_(?P<version>[0-9]{2})(?P<extension>\.TAB$)",
        "ck": r"^(?P<mission>ROS)?_?(?P<data>(ATNM|CATT|LATT|RATM|ROTT))?(?:_+(?P<token1>[A-Z0-9]+))?(?:_+(?P<token2>[A-Z0-9]+))?(?:_+(?P<token3>[A-Z0-9]+))?(?:_+(?P<token4>[A-Z0-9]+))?$",
        "mes": r"()",
        "spk": r"()",
    },
}

CASSINI_TITAN_FLYBY_DICT = {
    "T011": {
        "experiment": "tigr3",
        "pds_repo": "cors_0133",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "27.02.",
        "doy": 58,
        "year": 2006,
    },
    "T022": {
        "experiment": "tigr6",
        "pds_repo": "cors_0168",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "28.12.",
        "doy": 362,
        "year": 2006,
    },
    "T033": {
        "experiment": "tigr8",
        "pds_repo": "cors_0176",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "29.06.",
        "doy": 180,
        "year": 2007,
    },
    "T045": {
        "experiment": "tigr11",
        "pds_repo": "cors_0239",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "30.07.",
        "doy": 212,
        "year": 2008,
    },
    "T068": {
        "experiment": "tigr15",
        "pds_repo": "cors_0320",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "19.05.",
        "doy": 139,
        "year": 2010,
    },
    "T074": {
        "experiment": "tigr16",
        "pds_repo": "cors_0349, cors_0350",
        "ancillary_repo": "cors_0349",
        "cumindex_repo": "cors_0350",
        "date": "18.02.",
        "doy": 49,
        "year": 2011,
    },
    "T089": {
        "experiment": "tigr17",
        "pds_repo": "cors_0394, cors_0395",
        "ancillary_repo": "cors_0394",
        "cumindex_repo": "cors_0395",
        "date": "16.02.",
        "doy": 47,
        "year": 2013,
    },
    "T099": {
        "experiment": "tigr18",
        "pds_repo": "cors_0487, cors_0488",
        "ancillary_repo": "cors_0487",
        "cumindex_repo": "cors_0488",
        "date": "06.03.",
        "doy": 65,
        "year": 2014,
    },
    "T110": {
        "experiment": "tigr19",
        "pds_repo": "cors_0525",
        "ancillary_repo": None,
        "cumindex_repo": None,
        "date": "16.03.",
        "doy": 75,
        "year": 2015,
    },
    "T122": {
        "experiment": "tigr20",
        "pds_repo": "cors_0566, cors_0567",
        "ancillary_repo": "cors_0566",
        "cumindex_repo": "cors_0567",
        "date": "09.08.",
        "doy": 222,
        "year": 2016,
    },
}
