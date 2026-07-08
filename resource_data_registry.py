# Defines the structure and provenance of resource files.
TUDAT_RESOURCES_VERSION = "2.4"
TUDAT_RESOURCES_URL = f"https://raw.githubusercontent.com/tudat-team/tudat-resources/refs/tags/{TUDAT_RESOURCES_VERSION}/resource/"
DATA_REGISTRY = {
    "atmosphere_tables": {
        "dtm_mars.dat": "",
        "MCDMeanAtmosphere.dat": "",
        "nrlmsise00_validation_data.pkl": "",
        "USSA1976Until100kmPer100mUntil1000kmPer1000m_wHR_GC.dat": "",
        "USSA1976Until100kmPer100mUntil1000kmPer1000m.dat": "",
        "USSA1976Until86kmPer100m.dat": "",
        "MCDMeanAtmosphereTimeAverage": {
            "density.dat": "",
            "gasConstant.dat": "",
            "pressure.dat": "",
            "specificHeatRatio.dat": "",
            "temperature.dat": ""
        },
        "vmf3": {
            "anm_bh.txt": "",
            "anm_bw.txt": "",
            "anm_ch.txt": "",
            "anm_cw.txt": "",
            "bnm_bh.txt": "",
            "bnm_bw.txt": "",
            "bnm_ch.txt": "",
            "bnm_cw.txt": ""
        }
    },
    "earth_deformation": {
        "diurnalDisplacementFrequencyDependence.txt": "",
        "diurnalDisplacementFrequencyDependence2.txt": "",
        "ilrsStations.blq": "",
        "longPeriodDisplacementFrequencyDependence.txt": "",
        "oceanTidetest.blq": "",
    },
    "earth_orientation": {
        "eopc04_08_IAU2000.62-now.txt": "",
        "eopc04_14_IAU2000.62-now.txt": "",
        "historicalDeltaT.txt": "",
        "polarMotionLibrationAmplitudes.txt": "",
        "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt": "",
        "polarMotionLibrationFundamentalArgumentMultipliers.txt": "",
        "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt": "",
        "polarMotionOceanTidesAmplitudes.txt": "",
        "polarMotionOceanTidesFundamentalArgumentMultipliers.txt": "",
        "utcLibrationAmplitudes.txt": "",
        "utcLibrationFundamentalArgumentMultipliers.txt": "",
        "utcOceanTidesAmplitudes.txt": "",
        "utcOceanTidesFundamentalArgumentMultipliers.txt": "",
    },
    "ephemeris": {
        "eros_obs.tx": "",
        "p_elem_t1.txt": "",
        "p_elem_t2.txt": ""
    },
    "gravity_models": {
        "Earth": {
            "egm96.txt": "",
            "ggm02c.txt": "",
            "ggm02s.txt": "",
            "GOCO05c.txt": ""
        },
        "Mars": {
            "jgmro120d.txt": ""
        },
        "Mercury": {
            "jgmess_160a_sha.tab": ""
        },
        "Moon": {
            "gggrx_1200l_sha.tab": "",
            "glgm3150.txt": "",
            "lpe200.txt": ""
        },
        "Venus": {
            "shgj180u.a01": ""
        }
    },
    "quadrature": {
        "gaussianNodes.txt": "",
        "gaussianWeights.txt": ""
    },
    "space_weather": {
        "sw19571001.txt": ""
    },
    "spice_kernels": {
        "celestrak.json": "",
        "codes_300ast_20100725.bsp": "",
        "codes_300ast_20100725.cmt": "",
        "codes_300ast_20100725.tf": "",
        "de430_mar097_small.bsp": "",
        "de440.bsp": "",
        "earth_200101_990825_predict.bpc": "",
        "earth_200101_990825_predict.cmt": "",
        "earth_720101_230601.bpc": "",
        "earth_720101_230601.cmt": "",
        "earth_fixed.tf": "",
        "earth_latest_high_prec.bpc": "",
        "earth_latest_high_prec.cmt": "",
        "gm_de431.tpc": "",
        "inpop19a_TCB_m100_p100_asc": {
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
        "inpop19a_TDB_m100_p100_spice.bpc": "",
        "inpop19a_TDB_m100_p100_spice.bsp": "",
        "inpop19a_TDB_m100_p100_spice.tpc": "",
        "juice_mat_crema_4_0_20220601_20330626_v01.bsp": "",
        "mars_iau2000_v1.tpc": "",
        "moon_080317.tf": "",
        "moon_assoc_pa.tf": "",
        "moon_de440_200625.tf": "",
        "moon_pa_de440_200625.bpc": "",
        "moon_pa_de440_200625.cmt": "",
        "naif.json": "",
        "naif0012.tls": "",
        "NOE-4-2020.bsp": "",
        "NOE-4-2020.tpc": "",
        "NOE-5-2021.bsp": "",
        "NOE-5-2021.tpc": "",
        "NOE-6-2018-MAIN-v2.bsp": "",
        "NOE-6-2018-MAIN-v2.tpc": "",
        "pck00010.tpc": "",
        "pck00011.tpc": "",
        "tudat_merged_spk_kernel.bsp": "",
        "tudat_merged_spk_kernel.inp": "",
        "tudat_space.json": ""
    },
    "star_catalog_biases": {
        "debias_2018": {
            "bias.dat": "",
            "README.txt": "",
            "tiles.dat": ""
        }
    },
    "station_locations": {
        "glo.sit": "",
        "glo.vel": "",
        "mpc_codes.dat": "",
        "mpc.sit": "",
        "mpc.vel": "",
        "ns_codes.dat": "",
    }
}