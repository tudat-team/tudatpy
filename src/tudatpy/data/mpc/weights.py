import pandas as pd
import numpy as np
from astropy.units import Quantity
import astropy.units as u
from astroquery.mpc import MPC


def get_weights_VFCC17(
    mpc_table: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """Retrieves observation weights using the weighting scheme presented in
    "Statistical analysis of astrometric errors for the most productive
    asteroid surveys" by Veres et al. (2017).

    Observation types: "x", "X", "V", "v", "W", "w", "R", "r", "Q", "q", "O",
    are not described by the paper and receive a placeholder weight of 1/100
    if provided.

    Note that VFCC17 produces a single positional weight per observation, not
    separate RA/DEC weights - the two returned arrays are therefore identical,
    provided as a pair only for interface symmetry with get_biases_EFCC18.

    Parameters
    ----------
    mpc_table : pd.DataFrame
        Table retrieved by calling the mpc.BatchMPC.table property. Must contain
        'number', 'epoch' (Julian Days), 'note2', 'observatory' and 'catalog' columns.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Right Ascension weights, Declination weights (identical values), one entry
        per row of `mpc_table`, in the same order.
    """

    table = mpc_table.copy()
    table["observatory"] = table["observatory"].astype(str).str.strip().str.zfill(3)
    table["number"] = table["number"].astype(str).str.strip()

    # NOTE 1000 is a placeholder. The following observation types are not processed and receive the placeholder value:
    # first_discoveries = ["x", "X"], roaming = ["V", "v", "W", "w"], radar = ["R", "r", "Q", "q"], offset = ["O"]
    table = table.assign(inv_w=lambda _: 1000)

    # get an approximate timezone based on the observatory code's longitude
    observatories_table = MPC.get_observatory_codes().to_pandas()
    observatories_table["Code"] = observatories_table["Code"].astype(str).str.strip().str.zfill(3)
    observatories_table = (
        observatories_table.assign(lon_wrapping=lambda x: (x.Longitude + 180) % 360 - 180)
        .assign(approx_tz=lambda x: ((x.lon_wrapping / 180) * 12))
        .assign(jd_tz=lambda x: (x.lon_wrapping / 360).fillna(0))
        .loc[:, ["Code", "approx_tz", "jd_tz"]]
    )

    # the reset_index + set_index is a cheeky way to keep the original table index
    # this prevents mismatches when adding weights to the tables downstream
    table = (
        pd.merge(
            how="left",
            left=table.reset_index(),
            right=observatories_table,
            left_on="observatory",
            right_on="Code",
        )
        .drop("Code", axis=1)
        .set_index("index")
    )

    # create column that modifies the JD with an approximate timezone.
    # the time JD.50 will then be the approximate midnight at that timezone.
    # if we then take the int, we can group by this number to get all observations that night.
    table = table.assign(
        epochJD_tz_int=lambda x: np.floor(x.epoch + x.jd_tz)
    )  # epoch is in Julian Days.

    # Below are the weights applied as described per table in:
    # https://www.sciencedirect.com/science/article/pii/S0019103517301987

    # NON-CCD
    # ###################
    # TABLE 5: Non-CCD residuals
    # Conditions for Table 5
    # 1890-1-1 and 1950-1-1
    pre_1890 = table.epoch <= 2411368.0  # check on Julian Period
    between_1890_1950 = (table.epoch > 2411368.0) & (table.epoch <= 2433282.0)
    after_1950 = table.epoch > 2433282.0

    photographic = table.note2.isin([np.nan, "P", "A", "N", "Z"])
    occultations = table.note2 == "E"
    hipparcos = table.note2 == "H"
    transit_circle = table.note2 == "T"
    encoder = table.note2 == "e"
    micrometer = table.note2 == "M"
    satellite = table.note2.isin(["S", "s"])
    multinormal_place = table.note2 == "n"

    # Apply Table 5
    table.inv_w = table.inv_w.mask((photographic & pre_1890), 10.0)
    table.inv_w = table.inv_w.mask((photographic & between_1890_1950), 5.0)
    table.inv_w = table.inv_w.mask((photographic & after_1950), 2.5)

    table.inv_w = table.inv_w.mask((occultations), 0.2)
    table.inv_w = table.inv_w.mask((hipparcos), 0.2)
    table.inv_w = table.inv_w.mask((transit_circle), 0.5)
    table.inv_w = table.inv_w.mask((encoder), 0.75)
    table.inv_w = table.inv_w.mask((micrometer), 2.0)
    table.inv_w = table.inv_w.mask((satellite), 1.5)
    table.inv_w = table.inv_w.mask((multinormal_place), 1.0)
    # ###################

    # CCD
    # ###################
    # TABLE 3: astroid observers
    # Table 3 conditions:
    # NOTE for now CCD will receive the same base weighting as CCD.
    # CMOS sensors for astronomy is a new development.
    ccd = table.note2.isin(["C", "c", "D"])
    cmos = table.note2.isin(["B"])
    ccd = ccd | cmos

    tab3_no_catalog = table.catalog.isin(["unknown", np.nan])

    # Apply Table 3:
    table.inv_w = table.inv_w.mask((ccd & ~tab3_no_catalog), 1.0)
    table.inv_w = table.inv_w.mask((ccd & tab3_no_catalog), 1.5)

    table.inv_w = table.inv_w.mask((table.observatory == "704"), 1.0)
    table.inv_w = table.inv_w.mask((table.observatory == "G96"), 0.5)
    table.inv_w = table.inv_w.mask((table.observatory == "F51"), 0.2)
    table.inv_w = table.inv_w.mask((table.observatory == "G45"), 0.6)
    table.inv_w = table.inv_w.mask((table.observatory == "699"), 0.8)
    table.inv_w = table.inv_w.mask((table.observatory == "D29"), 0.75)

    table.inv_w = table.inv_w.mask((table.observatory == "C51"), 1.0)
    table.inv_w = table.inv_w.mask((table.observatory == "E12"), 0.75)
    table.inv_w = table.inv_w.mask((table.observatory == "608"), 0.6)
    table.inv_w = table.inv_w.mask((table.observatory == "J75"), 1.0)
    # ###################

    # ###################
    # TABLE 2: epoch dependant residuals
    # Table 2 conditions:
    tab2_703_epoch = table.epoch < 2456658.0  # 2014-1-1
    tab2_691_epoch = table.epoch < 2452640.0  # 2003-1-1
    tab2_644_epoch = table.epoch < 2452883.0  # 2003-9-1

    # Apply Table 2:
    table.inv_w = table.inv_w.mask((tab2_703_epoch & (table.observatory == "703")), 1.0)
    table.inv_w = table.inv_w.mask((~tab2_703_epoch & (table.observatory == "703")), 0.8)
    table.inv_w = table.inv_w.mask((tab2_691_epoch & (table.observatory == "691")), 0.6)
    table.inv_w = table.inv_w.mask((~tab2_691_epoch & (table.observatory == "691")), 0.5)
    table.inv_w = table.inv_w.mask((tab2_644_epoch & (table.observatory == "644")), 0.6)
    table.inv_w = table.inv_w.mask((~tab2_644_epoch & (table.observatory == "644")), 0.4)
    # ###################

    # ###################
    # TABLE 4: Catalog dependant + NEO observers
    # Table 4 conditions:
    # NOTE these are the LCO observatories described in the paper
    LCO_original = [
        "K92",
        "K93",
        "Q63",
        "Q64",
        "V37",
        "W84",
        "W85",
        "W86",
        "W87",
        "K91",
        "E10",
        "F65",
    ]
    # NOTE these are additional LCO observatories, mostly online after publication
    # Aqawan and Clamshell (0.4m) style observatories have not been included due to their previous omssion in the paper.
    # Those thus assume the default inv_weight of 1.0
    # Since remaining 1m telescopes are duplicates of existing observatories, they are included.
    # Dates retrieved from: https://sbnmpc.astro.umd.edu/mpecwatch/obs.html
    # See also: https://lco.global/observatory/sites/mpccodes/
    LCO_new = [
        "V39",  # online 2019 1m
        "Z24",  # online 2021 1m
        "Z31",  # online 2021 1m
    ]

    LCO_obs = LCO_original + LCO_new
    MAUNAKEA_obs = ["T09", "T12", "T14"]
    tab4_LCO_observatories = table.observatory.isin(LCO_obs)
    # NOTE see catalog codes info: https://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
    tab4_Catalog_UCAC4 = table.catalog == "q"
    tab4_Catalog_PPMXL = table.catalog == "t"
    tab4_Catalog_GAIA = table.catalog.isin(["U", "V", "W", "X", "3", "6"])
    tab4_Catalog_USNOB12 = table.catalog.isin(["o", "s"])

    tab4_G83_UCAC4_PPMXL = (table.observatory == "G83") & (tab4_Catalog_UCAC4 | tab4_Catalog_PPMXL)
    tab4_G83_GAIA = (table.observatory == "G83") & tab4_Catalog_GAIA

    tab4_Y28_GAIA_PPMXL = (table.observatory == "Y28") & (tab4_Catalog_PPMXL & tab4_Catalog_GAIA)
    tab4_568_USNOB = (table.observatory == "568") & tab4_Catalog_USNOB12
    tab4_568_GAIA = (table.observatory == "568") & tab4_Catalog_GAIA
    tab4_568_PPMXL = (table.observatory == "568") & tab4_Catalog_PPMXL
    tab4_T09_T12_T14_GAIA = (table.observatory.isin(MAUNAKEA_obs)) & tab4_Catalog_GAIA
    tab4_309_UCAC4_PPMXL = (table.observatory == "309") & (tab4_Catalog_UCAC4 | tab4_Catalog_PPMXL)
    tab4_309_GAIA = (table.observatory == "309") & tab4_Catalog_GAIA

    # Apply Table 4:
    table.inv_w = table.inv_w.mask(table.observatory == "645", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "673", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "689", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "950", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "H01", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "J04", 0.4)
    table.inv_w = table.inv_w.mask(tab4_G83_UCAC4_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_G83_GAIA, 0.2)
    table.inv_w = table.inv_w.mask(tab4_LCO_observatories, 0.4)
    table.inv_w = table.inv_w.mask(table.observatory == "W84", 0.5)  # LCO includes W84?

    table.inv_w = table.inv_w.mask(tab4_Y28_GAIA_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_568_USNOB, 0.5)
    table.inv_w = table.inv_w.mask(tab4_568_GAIA, 0.1)
    table.inv_w = table.inv_w.mask(tab4_568_PPMXL, 0.2)
    table.inv_w = table.inv_w.mask(tab4_T09_T12_T14_GAIA, 0.1)
    table.inv_w = table.inv_w.mask(tab4_309_UCAC4_PPMXL, 0.3)
    table.inv_w = table.inv_w.mask(tab4_309_GAIA, 0.2)
    # ###################

    # Transform residual into weight:
    table = table.assign(
        weight_pre=lambda x: 1 / np.square(Quantity(x.inv_w, unit=u.arcsec).to(u.rad).value)
    )

    # Reduce weight if there are more than 4 observations that night:
    # NOTE For satellites, this is done per julian day.
    # calculate number of observations on one night at one observatory for one target:
    table = table.assign(
        observations_on_epoch=lambda x: x.groupby(
            ["epochJD_tz_int", "observatory", "number"]
        ).epoch.transform("count")
    )

    # Veres states that, for N > 4,  "uncertainties have to be inflated by a factor of sqrt(N)".
    # This means sigma_new = sigma_old * sqrt(N)
    # Which also means: weight_new = 1/(sigma_new^2) = W_old/N
    # Since we want the weight to be 1 when N = 4, we end up dividing by 4.
    # How this is done in practice: we divide the weight by N if N > 4, else 1.
    table = table.assign(mult_obs_deweight=lambda x: np.maximum(x.observations_on_epoch / 4, 1.0))

    table = table.assign(weight=lambda x: x.weight_pre / x.mult_obs_deweight)

    weight = table.weight.to_numpy()
    RA_weight = weight
    DEC_weight = weight

    return RA_weight, DEC_weight
