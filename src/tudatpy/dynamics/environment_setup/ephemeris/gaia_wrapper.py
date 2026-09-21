"""Environment-model helpers for Gaia archive data."""

from tudatpy.dynamics.environment_setup import ephemeris


def gaia_from_astrometry(gaia_astrometry, geocentric: bool = True):
    """Create tabulated ephemeris settings from loaded Gaia astrometry.

    Parameters
    ----------
    gaia_astrometry : tudatpy.data_input.tracking_data.gaia.GaiaAstrometry
        Loaded Gaia astrometry containing spacecraft state vectors.
    geocentric : bool, default True
        Use the geocentric state vectors and Earth as frame origin. If false,
        use barycentric vectors and the solar-system barycentre.
    """
    return ephemeris.tabulated(
        gaia_astrometry.get_gaia_state_history(geocentric),
        frame_origin="Earth" if geocentric else "SSB",
        frame_orientation="J2000",
    )
