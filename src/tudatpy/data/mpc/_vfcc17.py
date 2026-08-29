from __future__ import annotations

import warnings

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

_MIGRATION_GUIDE_URL = (
    "https://docs.tudat.space/en/latest/user-guide/project-updates/migration-guide.html"
)


def get_weights_VFCC17(
    MPC_codes: pd.Series | list | np.ndarray | None = None,
    epoch: pd.Series | list | np.ndarray | None = None,
    observation_type: pd.Series | list | np.ndarray | None = None,
    observatory: pd.Series | list | np.ndarray | None = None,
    star_catalog: pd.Series | list | np.ndarray | None = None,
    mpc_table: pd.DataFrame | None = None,
    return_full_table=False,
) -> np.ndarray | pd.DataFrame:
    """Return MPC astrometry weights using the Veres et al. (2017) scheme.

    This is the deprecated compatibility implementation that used to live in
    :mod:`tudatpy.data.mpc`. There is no one-to-one replacement in the new
    data-input workflow that returns the same direct weight array or
    intermediate table.
    """
    import numpy as np
    import pandas as pd
    from astropy.units import Quantity
    import astropy.units as u
    from astroquery.mpc import MPC

    warnings.warn(
        (
            "tudatpy.data.mpc.get_weights_VFCC17 is deprecated. "
            "There is no one-to-one equivalent in the new data_input workflow. "
            "The old function remains supported during the deprecation period. "
            "In the new setup, request VFCC17 weighting with add_weights=True "
            "when creating optical TrackingData from MPC, 80-column, pandas, "
            "or astropy data; direct weight arrays and intermediate weighting "
            "tables are no longer exposed as a standalone public interface. "
            f"See the TudatPy migration guide: {_MIGRATION_GUIDE_URL}"
        ),
        DeprecationWarning,
        stacklevel=2,
    )

    if (
        (mpc_table is None)
        and (epoch is not None)
        and (observation_type is not None)
        and (observatory is not None)
        and (star_catalog is not None)
    ):
        if not (len(epoch) == len(observation_type) == len(observatory) == len(star_catalog)):
            raise ValueError("All inputs must have same size")

        table = pd.DataFrame.from_dict(
            {
                "number": MPC_codes,
                "epoch": epoch,
                "note2": observation_type,
                "observatory": observatory,
                "catalog": star_catalog,
            }
        )
    elif (
        (mpc_table is not None)
        and (epoch is None)
        and (observation_type is None)
        and (observatory is None)
        and (star_catalog is None)
    ):
        table = mpc_table.copy()
    else:
        raise ValueError(
            "Must provide either parameters: `epoch`, `observation_type`, "
            "`observatory` and `star_catalog` OR `mpc_table`."
        )

    table["observatory"] = table["observatory"].astype(str).str.strip().str.zfill(3)
    table["number"] = table["number"].astype(str).str.strip()
    table = table.assign(inv_w=lambda _: 1000)

    observatories_table = MPC.get_observatory_codes().to_pandas()
    observatories_table["Code"] = observatories_table["Code"].astype(str).str.strip().str.zfill(3)
    observatories_table = (
        observatories_table.assign(lon_wrapping=lambda x: (x.Longitude + 180) % 360 - 180)
        .assign(approx_tz=lambda x: ((x.lon_wrapping / 180) * 12))
        .assign(jd_tz=lambda x: (x.lon_wrapping / 360).fillna(0))
        .loc[:, ["Code", "approx_tz", "jd_tz"]]
    )

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
    table = table.assign(epochJD_tz_int=lambda x: np.floor(x.epoch + x.jd_tz))

    pre_1890 = table.epoch <= 2411368.0
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

    ccd = table.note2.isin(["C", "c", "D"]) | table.note2.isin(["B"])
    tab3_no_catalog = table.catalog.isin(["unknown", np.nan])

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

    tab2_703_epoch = table.epoch < 2456658.0
    tab2_691_epoch = table.epoch < 2452640.0
    tab2_644_epoch = table.epoch < 2452883.0

    table.inv_w = table.inv_w.mask((tab2_703_epoch & (table.observatory == "703")), 1.0)
    table.inv_w = table.inv_w.mask((~tab2_703_epoch & (table.observatory == "703")), 0.8)
    table.inv_w = table.inv_w.mask((tab2_691_epoch & (table.observatory == "691")), 0.6)
    table.inv_w = table.inv_w.mask((~tab2_691_epoch & (table.observatory == "691")), 0.5)
    table.inv_w = table.inv_w.mask((tab2_644_epoch & (table.observatory == "644")), 0.6)
    table.inv_w = table.inv_w.mask((~tab2_644_epoch & (table.observatory == "644")), 0.4)

    lco_original = [
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
    lco_new = ["V39", "Z24", "Z31"]

    lco_obs = lco_original + lco_new
    maunakea_obs = ["T09", "T12", "T14"]
    tab4_lco_observatories = table.observatory.isin(lco_obs)
    tab4_catalog_ucac4 = table.catalog == "q"
    tab4_catalog_ppmxl = table.catalog == "t"
    tab4_catalog_gaia = table.catalog.isin(["U", "V", "W", "X", "3", "6"])
    tab4_catalog_usnob12 = table.catalog.isin(["o", "s"])

    tab4_g83_ucac4_ppmxl = (table.observatory == "G83") & (tab4_catalog_ucac4 | tab4_catalog_ppmxl)
    tab4_g83_gaia = (table.observatory == "G83") & tab4_catalog_gaia

    tab4_y28_gaia_ppmxl = (table.observatory == "Y28") & (tab4_catalog_ppmxl & tab4_catalog_gaia)
    tab4_568_usnob = (table.observatory == "568") & tab4_catalog_usnob12
    tab4_568_gaia = (table.observatory == "568") & tab4_catalog_gaia
    tab4_568_ppmxl = (table.observatory == "568") & tab4_catalog_ppmxl
    tab4_t09_t12_t14_gaia = (table.observatory.isin(maunakea_obs)) & tab4_catalog_gaia
    tab4_309_ucac4_ppmxl = (table.observatory == "309") & (tab4_catalog_ucac4 | tab4_catalog_ppmxl)
    tab4_309_gaia = (table.observatory == "309") & tab4_catalog_gaia

    table.inv_w = table.inv_w.mask(table.observatory == "645", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "673", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "689", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "950", 0.5)
    table.inv_w = table.inv_w.mask(table.observatory == "H01", 0.3)
    table.inv_w = table.inv_w.mask(table.observatory == "J04", 0.4)
    table.inv_w = table.inv_w.mask(tab4_g83_ucac4_ppmxl, 0.3)
    table.inv_w = table.inv_w.mask(tab4_g83_gaia, 0.2)
    table.inv_w = table.inv_w.mask(tab4_lco_observatories, 0.4)
    table.inv_w = table.inv_w.mask(table.observatory == "W84", 0.5)

    table.inv_w = table.inv_w.mask(tab4_y28_gaia_ppmxl, 0.3)
    table.inv_w = table.inv_w.mask(tab4_568_usnob, 0.5)
    table.inv_w = table.inv_w.mask(tab4_568_gaia, 0.1)
    table.inv_w = table.inv_w.mask(tab4_568_ppmxl, 0.2)
    table.inv_w = table.inv_w.mask(tab4_t09_t12_t14_gaia, 0.1)
    table.inv_w = table.inv_w.mask(tab4_309_ucac4_ppmxl, 0.3)
    table.inv_w = table.inv_w.mask(tab4_309_gaia, 0.2)

    table = table.assign(
        weight_pre=lambda x: 1 / np.square(Quantity(x.inv_w, unit=u.arcsec).to(u.rad).value)
    )
    table = table.assign(
        observations_on_epoch=lambda x: x.groupby(
            ["epochJD_tz_int", "observatory", "number"]
        ).epoch.transform("count")
    )
    table = table.assign(mult_obs_deweight=lambda x: np.maximum(x.observations_on_epoch / 4, 1.0))
    table = table.assign(weight=lambda x: x.weight_pre / x.mult_obs_deweight)

    if return_full_table:
        return table
    return table.weight.to_numpy()
