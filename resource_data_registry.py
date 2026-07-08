# Defines the structure and provenance of resource files.
TUDAT_RESOURCES_VERSION = "2.4"
TUDAT_RESOURCES_URL = f"https://raw.githubusercontent.com/tudat-team/tudat-resources/refs/tags/{TUDAT_RESOURCES_VERSION}/resource/"
DATA_REGISTRY = {
    "atmosphere_tables": {
        "dtm_mars.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/dtm_mars.dat",
        "MCDMeanAtmosphere.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphere.dat",
        "nrlmsise00_validation_data.pkl": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/nrlmsise00_validation_data.pkl",
        "USSA1976Until100kmPer100mUntil1000kmPer1000m_wHR_GC.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/USSA1976Until100kmPer100mUntil1000kmPer1000m_wHR_GC.dat",
        "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat",
        "USSA1976Until86kmPer100m.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/USSA1976Until86kmPer100m.dat",
        "MCDMeanAtmosphereTimeAverage": {
            "density.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphereTimeAverage/density.dat",
            "gasConstant.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphereTimeAverage/gasConstant.dat",
            "pressure.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphereTimeAverage/pressure.dat",
            "specificHeatRatio.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat",
            "temperature.dat": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/MCDMeanAtmosphereTimeAverage/temperature.dat"
        },
        "vmf3": {
            "anm_bh.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/anm_bh.txt",
            "anm_bw.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/anm_bw.txt",
            "anm_ch.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/anm_ch.txt",
            "anm_cw.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/anm_cw.txt",
            "bnm_bh.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/bnm_bh.txt",
            "bnm_bw.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/bnm_bw.txt",
            "bnm_ch.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/bnm_ch.txt",
            "bnm_cw.txt": f"{TUDAT_RESOURCES_URL}/atmosphere_tables/vmf3/bnm_cw.txt"
        }
    },
    "earth_deformation": {
        "diurnalDisplacementFrequencyDependence.txt": f"{TUDAT_RESOURCES_URL}/earth_deformation/diurnalDisplacementFrequencyDependence.txt",
        "diurnalDisplacementFrequencyDependence2.txt": f"{TUDAT_RESOURCES_URL}/earth_deformation/diurnalDisplacementFrequencyDependence2.txt",
        "ilrsStations.blq": f"{TUDAT_RESOURCES_URL}/earth_deformation/ilrsStations.blq",
        "longPeriodDisplacementFrequencyDependence.txt": f"{TUDAT_RESOURCES_URL}/earth_deformation/longPeriodDisplacementFrequencyDependence.txt",
        "oceanTidetest.blq": f"{TUDAT_RESOURCES_URL}/earth_deformation/oceanTidetest.blq",
    },
    "earth_orientation": {
        "eopc04_08_IAU2000.62-now.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/eopc04_08_IAU2000.62-now.txt",
        "eopc04_14_IAU2000.62-now.txt": "https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt",
        "historicalDeltaT.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/historicalDeltaT.txt",
        "polarMotionLibrationAmplitudes.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionLibrationAmplitudes.txt",
        "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
        "polarMotionLibrationFundamentalArgumentMultipliers.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionLibrationFundamentalArgumentMultipliers.txt",
        "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
        "polarMotionOceanTidesAmplitudes.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionOceanTidesAmplitudes.txt",
        "polarMotionOceanTidesFundamentalArgumentMultipliers.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/polarMotionOceanTidesFundamentalArgumentMultipliers.txt",
        "utcLibrationAmplitudes.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/utcLibrationAmplitudes.txt",
        "utcLibrationFundamentalArgumentMultipliers.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/utcLibrationFundamentalArgumentMultipliers.txt",
        "utcOceanTidesAmplitudes.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/utcOceanTidesAmplitudes.txt",
        "utcOceanTidesFundamentalArgumentMultipliers.txt": f"{TUDAT_RESOURCES_URL}/earth_orientation/utcOceanTidesFundamentalArgumentMultipliers.txt",
    },
    "ephemeris": {
        "eros_obs.txt": f"{TUDAT_RESOURCES_URL}/ephemeris/eros_obs.txt",
        "p_elem_t1.txt": f"{TUDAT_RESOURCES_URL}/ephemeris/p_elem_t1.txt",
        "p_elem_t2.txt": f"{TUDAT_RESOURCES_URL}/ephemeris/p_elem_t2.txt"
    },
    "gravity_models": {
        "Earth": {
            "egm96.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Earth/egm96.txt",
            "ggm02c.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Earth/ggm02c.txt",
            "ggm02s.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Earth/ggm02s.txt",
            "GOCO05c.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Earth/GOCO05c.txt"
        },
        "Mars": {
            "jgmro120d.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Mars/jgmro120d.txt"
        },
        "Mercury": {
            "jgmess_160a_sha.tab": "https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.tab"
        },
        "Moon": {
            "gggrx_1200l_sha.tab": "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200l_sha.tab",
            "glgm3150.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Moon/glgm3150.txt",
            "lpe200.txt": f"{TUDAT_RESOURCES_URL}/gravity_models/Moon/lpe200.txt"
        },
        "Venus": {
            "shgj180u.a01": "https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/shgj180u.a01"
        }
    },
    "quadrature": {
        "gaussianNodes.txt": f"{TUDAT_RESOURCES_URL}/quadrature/gaussianNodes.txt",
        "gaussianWeights.txt": f"{TUDAT_RESOURCES_URL}/quadrature/gaussianWeights.txt"
    },
    "space_weather": {
        "sw19571001.txt": "https://celestrak.org/SpaceData/sw19571001.txt"
    },
    "spice_kernels": {
        "celestrak.json": f"{TUDAT_RESOURCES_URL}/spice_kernels/celestrak.json",
        "codes_300ast_20100725.bsp": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.bsp",
        "codes_300ast_20100725.cmt": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.cmt",
        "codes_300ast_20100725.tf": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.tf",
        "de430_mar097_small.bsp": f"{TUDAT_RESOURCES_URL}/spice_kernels/de430_mar097_small.bsp",
        "de440.bsp": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
        "earth_200101_990825_predict.bpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_200101_990825_predict.bpc",
        "earth_200101_990825_predict.cmt": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_200101_990825_predict.cmt",
        "earth_720101_230601.bpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_720101_230601.bpc",
        "earth_720101_230601.cmt": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_720101_230601.cmt",
        "earth_fixed.tf": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_fixed.tf",
        "earth_latest_high_prec.bpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc",
        "earth_latest_high_prec.cmt": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.cmt",
        "gm_de431.tpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc",
        "inpop19a_TCB_m100_p100_asc": {
            "__ALL__": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TCB_m100_p100_asc.tar.gz"
            "inpop19a_TCB_m100_p100_asc_header.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Ear.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_EMB.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Jup.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Lib.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Mar.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Mer.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Moo.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Nep.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Plu.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Sat.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Sun.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_TCG.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Ura.asc": "",
            "inpop19a_TCB_m100_p100_asc_pos_Ven.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Ear.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_EMB.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Jup.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Lib.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Mar.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Mer.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Moo.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Nep.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Plu.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Sat.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Sun.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Ura.asc": "",
            "inpop19a_TCB_m100_p100_asc_vel_Ven.asc": ""
        },
        "inpop19a_TDB_m100_p100_asc": {
            "__ALL__": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TDB_m100_p100_asc.tar.gz"
            "inpop19a_TDB_m100_p100_asc_header.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Ear.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_EMB.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Jup.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Lib.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Mar.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Mer.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Moo.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Nep.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Plu.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Sat.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Sun.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_TT.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Ura.asc": "",
            "inpop19a_TDB_m100_p100_asc_pos_Ven.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Ear.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_EMB.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Jup.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Lib.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Mar.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Mer.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Moo.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Nep.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Plu.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Sat.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Sun.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Ura.asc": "",
            "inpop19a_TDB_m100_p100_asc_vel_Ven.asc": ""
        },
        "inpop19a_TDB_m100_p100_spice": {
            "__ALL__": "https://ftp.imcce.fr/pub/ephem/planets/inpop19a/inpop19a_TDB_m100_p100_spice.tar.gz",
            "inpop19a_TDB_m100_p100_spice.bpc": "",
            "inpop19a_TDB_m100_p100_spice.bsp": "",
            "inpop19a_TDB_m100_p100_spice.tpc": ""
        }
        "juice_mat_crema_4_0_20220601_20330626_v01.bsp": "https://spiftp.esac.esa.int/data/SPICE/JUICE/kernels/spk/juice_mat_crema_4_0_20220601_20330626_v01.bsp",
        "mars_iau2000_v1.tpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/mars_iau2000_v1.tpc",
        "moon_080317.tf": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf",
        "moon_assoc_pa.tf": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_assoc_pa.tf",
        "moon_de440_200625.tf": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/a_old_versions/moon_de440_200625.tf",
        "moon_pa_de440_200625.bpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc",
        "moon_pa_de440_200625.cmt": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.cmt",
        "naif.json": f"{TUDAT_RESOURCES_URL}/spice_kernels/naif.json",
        "naif0012.tls": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
        "NOE-4-2020.bsp": "https://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/NOE-4-2020.bsp",
        "NOE-4-2020.tpc": "https://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/NOE-4-2020.tpc",
        "NOE-5-2021.bsp": "https://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/NOE-5-2021.bsp",
        "NOE-5-2021.tpc": "https://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/NOE-5-2021.tpc",
        "NOE-6-2018-MAIN-v2.bsp": "https://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/NOE-6-2018-MAIN-v2.bsp",
        "NOE-6-2018-MAIN-v2.tpc": "https://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/NOE-6-2018-MAIN-v2.tpc",
        "pck00010.tpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
        "pck00011.tpc": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc",
        "tudat_merged_spk_kernel.bsp": f"{TUDAT_RESOURCES_URL}/spice_kernels/tudat_merged_spk_kernel.bsp",
        "tudat_merged_spk_kernel.inp": f"{TUDAT_RESOURCES_URL}/spice_kernels/tudat_merged_spk_kernel.inp",
        "tudat_space.json": f"{TUDAT_RESOURCES_URL}/spice_kernels/tudat_space.json"
    },
    "star_catalog_biases": {
        "debias_2018": {
            "bias.dat": f"{TUDAT_RESOURCES_URL}/star_catalog_biases/debias_2018/bias.dat",
            "README.txt": f"{TUDAT_RESOURCES_URL}/star_catalog_biases/debias_2018/README.txt",
            "tiles.dat": f"{TUDAT_RESOURCES_URL}/star_catalog_biases/debias_2018/tiles.dat"
        }
    },
    "station_locations": {
        "glo.sit": f"{TUDAT_RESOURCES_URL}/station_locations/glo.sit",
        "glo.vel": f"{TUDAT_RESOURCES_URL}/station_locations/glo.vel",
        "mpc_codes.dat": f"{TUDAT_RESOURCES_URL}/station_locations/mpc_codes.dat",
        "mpc.sit": f"{TUDAT_RESOURCES_URL}/station_locations/mpc.sit",
        "mpc.vel": f"{TUDAT_RESOURCES_URL}/station_locations/mpc.vel",
        "ns_codes.dat": f"{TUDAT_RESOURCES_URL}/station_locations/ns_codes.dat",
    }
}
