"""Read the two real loading grids and construct settings, without model evaluation."""

from pathlib import Path
from tudatpy.data_input.environment_data import atmospheric_loading
from tudatpy.dynamics import environment_setup

REFERENCE = Path("/home/dominic/Tudat/earth-station-model-work/references")
shape = environment_setup.shape_deformation

# Official Desai grid is already in CM. Native model uses the selected Earth rotation converter and 2018 secular pole.
ocean_pole = shape.ocean_pole_tidal(str(REFERENCE / "ocean-pole-tide/opoleloadcoefcmcor.txt"))

# The original RP03 CE grid uses millimetre coefficients; the reader converts them to metres.
metadata = atmospheric_loading.AtmosphericTideMetadata()
metadata.reference_origin = atmospheric_loading.AtmosphericLoadingReferenceOrigin.solid_earth_centre
metadata.source_profile = (
    atmospheric_loading.AtmosphericTideSourceProfile.iers2010_rp03_no_ib_farrell
)
metadata.latitude_convention = atmospheric_loading.AtmosphericTideLatitudeConvention.geocentric
metadata.displacement_basis = (
    atmospheric_loading.AtmosphericTideDisplacementBasis.geocentric_radial_north_east
)
metadata.expected_longitude_count = 361
metadata.expected_latitude_count = 181
grid = atmospheric_loading.read_atmospheric_tide_grid(
    str(REFERENCE / "atmospheric-tides/s1_s2_def_ce.dat"), metadata
)

# Native synthesis uses continuous UT1 and adds the identified IERS Table7.6 CE-to-CM translation once in Cartesian axes.
evaluation = shape.AtmosphericTideEvaluationSettings()
evaluation.lookup_latitude_convention = (
    atmospheric_loading.AtmosphericTideLatitudeConvention.geocentric
)
evaluation.displacement_basis = (
    atmospheric_loading.AtmosphericTideDisplacementBasis.geocentric_radial_north_east
)
evaluation.convert_to_centre_of_mass = True
atmospheric_tide = shape.atmospheric_tide_from_grid(grid, evaluation)

# Append these two settings once to the existing Earth body deformation list before creating the bodies.
# Existing body-deformation station motion then evaluates both stations at their own source coordinates.
print(type(ocean_pole).__name__, type(atmospheric_tide).__name__)
print(grid.metadata.model_name, grid.metadata.greens_function, grid.metadata.reference_origin)
