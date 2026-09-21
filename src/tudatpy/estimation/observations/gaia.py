"""Conversion of loaded Gaia astrometry to observation datasets."""


def create_observation_dataset_from_gaia_astrometry(gaia_astrometry, bodies):
    """Create a correlated observation dataset from loaded Gaia astrometry."""
    return gaia_astrometry.to_observation_dataset(bodies)
