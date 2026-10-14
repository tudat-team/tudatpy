/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_shape_deformation.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/ephemerides/ephemeris.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/simulation/environment_setup/createBodyDeformationModel.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace shape_deformation
{

void expose_shape_deformation_setup( py::module& m )
{
    py::class_< tss::BodyDeformationSettings, std::shared_ptr< tss::BodyDeformationSettings > >( m,
                                                                                                 "BodyDeformationSettings",
                                                                                                 R"doc(

         Base class for providing settings for body shape deformation model.

         Functional (base) class for settings of body shape deformation models that require no information in addition to their type.
         Body shape deformation model settings requiring additional information must be defined using an object derived from this class.


      )doc" );

    py::class_< tss::BasicSolidBodyDeformationSettings,
                std::shared_ptr< tss::BasicSolidBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m,
                                                "BasicSolidBodyDeformationSettings",
                                                R"doc(

         Class for defining model settings for simple tidal solid-body shape deformation.

         `BodyDeformationSettings` derived class for simple tidal solid-body shape deformation.




      )doc" );

    namespace ba = tudat::basic_astrodynamics;
    py::enum_< ba::AtmosphericTideEndpointPolicy >( m, "AtmosphericTideEndpointPolicy" )
            .value( "require_consistency", ba::AtmosphericTideEndpointPolicy::requireConsistency )
            .value( "average_seam_and_polar_cartesian_approximation",
                    ba::AtmosphericTideEndpointPolicy::averageSeamAndPolarCartesianApproximation );
    py::enum_< ba::AtmosphericLoadingOriginPolicy >( m, "AtmosphericLoadingOriginPolicy" )
            .value( "preserve_source_origin", ba::AtmosphericLoadingOriginPolicy::preserveSourceOrigin )
            .value( "convert_identified_rp03_ce_to_cm", ba::AtmosphericLoadingOriginPolicy::convertIdentifiedRp03CeToCm );
    py::class_< ba::AtmosphericTideEvaluationSettings >(
            m,
            "AtmosphericTideEvaluationSettings",
            "Independent grid lookup latitude and displacement basis declarations. Geodetic choices require Earth's configured shape. "
            "Strict endpoint consistency and preserved source origin are defaults." )
            .def( py::init<>( ) )
            .def_readwrite( "lookup_latitude_convention", &ba::AtmosphericTideEvaluationSettings::lookupLatitudeConvention )
            .def_readwrite( "displacement_basis", &ba::AtmosphericTideEvaluationSettings::displacementBasis )
            .def_readwrite( "endpoint_policy", &ba::AtmosphericTideEvaluationSettings::endpointPolicy )
            .def_readwrite( "origin_policy", &ba::AtmosphericTideEvaluationSettings::originPolicy );
    py::class_< tss::AtmosphericTideBodyDeformationSettings,
                std::shared_ptr< tss::AtmosphericTideBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m, "AtmosphericTideBodyDeformationSettings" );
    m.def( "atmospheric_tide_from_grid",
           &tss::atmosphericTideBodyShapeDeformationFromGrid,
           py::arg( "grid" ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           R"doc(Create conventional S1/S2 loading from preloaded SI grid records.

The environment epoch is TDB seconds since J2000. Phase is measured from UT1
midnight using Earth's shared continuous time converter. Grid latitude lookup
and displacement component basis must each be declared explicitly; the provider
latitude label does not resolve either convention. Source origin is preserved.
The optional identified RP03 CE-to-CM policy adds Table 7.6 XYZ translation once
after local rotation. Actual RP03 seam/polar rounding discrepancies require an
explicit approximation policy, which records their maximum size and preserves
raw records. Position evaluation works through inserted leap seconds. State
and velocity are analytic derivatives with respect to TDB, chaining the model's
continuous UT1/TAI rate and geocentric TT/TDB rate. Missing or unsupported rate
stencils raise for state/rate evaluation without disabling valid positions.
)doc" );
    m.def( "atmospheric_tide_from_stations",
           &tss::atmosphericTideBodyShapeDeformationFromStations,
           py::arg( "stations" ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           "Create S1/S2 loading from exact station identifiers and precomputed SI coefficients. Explicit displacement basis is required; "
           "source coordinates are not used for lookup. Epochs are TDB, phase UT1; analytic d/dTDB state/rates require a valid selected "
           "UT1-rate stencil and geocentric TT/TDB clock." );
    m.def( "atmospheric_tide_from_grid_file",
           &tss::atmosphericTideBodyShapeDeformationFromGridFile,
           py::arg( "file_path" ),
           py::arg( "metadata" ) = tudat::input_output::AtmosphericTideMetadata( ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           "Read a grid once during construction; use the same explicit geometry, endpoint and origin policies as "
           "atmospheric_tide_from_grid." );
    m.def( "atmospheric_tide_from_station_file",
           &tss::atmosphericTideBodyShapeDeformationFromStationFile,
           py::arg( "file_path" ),
           py::arg( "metadata" ) = tudat::input_output::AtmosphericTideMetadata( ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           "Read ATL station coefficients once during construction; use the same basis/origin policies as "
           "atmospheric_tide_from_stations." );

    m.def( "basic_solid_body_tidal",
           &tss::basicTidalBodyShapeDeformation,
           py::arg( "tide_raising_bodies" ),
           py::arg( "displacement_love_numbers" ),
           py::arg( "reference_radius" ) = TUDAT_NAN,
           R"doc(

 Function for creating basic tidal solid-body shape deformation

 Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{m}` and :math:`l_{m}` (with only :math:`m=2,3` currently supported). This function implements equations (7.5) and (7.6) of the `IERS Conventions 2010 <https://iers-conventions.obspm.fr/conventions_material.php>`_.


 Parameters
 ----------
 tide_raising_bodies : list[ string ]
     List of bodies that raise a tide that induces the shape variation.
 displacement_love_numbers : dict[int, tuple[float, float]]
     Dictionary of pairs. The dictionary key the spherical harmonic degree :math:`l` of the tidal deformation (2 or 3 are currently supported). The dictionary value is comprised of a pair of floats representing the :math:`h_{2}` and :math:`l_{2}` deformation Love numbers
 reference_radius : float, default = NaN
     Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
 Returns
 -------
 BasicSolidBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





 Examples
 --------
 In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:

 .. code-block:: python

   # Create Love numbers
   love_numbers = dict()
   love_numbers[2] = (0.6, 0.08)

   # Create tide raising bodies
   tide_raising_bodies = ["Sun", "Moon"]

   # Append shape variation settings to existing (default is empty) list
   body_settings.get( "Earth" ).shape_deformation_settings.append(
       environment_setup.shape_deformation.basic_solid_body_tidal(
           tide_raising_bodies, love_numbers ) )


     )doc" );

    m.def( "degree_two_basic_solid_body_tidal",
           &tss::degreeTwoBasicTidalBodyShapeDeformation,
           py::arg( "tide_raising_bodies" ),
           py::arg( "love_number" ),
           py::arg( "shida_number" ),
           py::arg( "reference_radius" ) = TUDAT_NAN,
           R"doc(

 Function for creating degree 2 basic tidal solid-body shape deformation

 Function for creating basic tidal solid-body shape deformation, computing the tidal shape variation due to any number of bodies causing the deformation, and a tidal response define by the deformation Love and Shida numbers :math:`h_{2}` and :math:`l_{2}`. This function implements equations (7.5) of the IERS Conventions 2010, and provides a simplified interface (for degree 2 only) of :func:`basic_solid_body_tidal`.


 Parameters
 ----------
 tide_raising_bodies : list[ string ]
     List of bodies that raise a tide that induces the shape variation.
 love_number : float
     Value of :math:`h_{2}` deformation Love number`
 shida_number : float
     Value of :math:`l_{2}` deformation Shida number`
 reference_radius : float, default = NaN
     Spherical harmonic reference radius of the deformed body. If this value is left undefined (e.g at NaN), the reference radius of the existing spherical harmonic gravity field of the deformed body is used.
 Returns
 -------
 BasicSolidBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings` class





 Examples
 --------
 In this example, we create a settings for degree 2 tidal deformation of the Earth due to the Sun and Moon:

 .. code-block:: python

   # Create Love numbers
   h2_love_number = 0.6
   l2_shida_number = 0.08

   # Create tide raising bodies
   tide_raising_bodies = ["Sun", "Moon"]

   # Append shape variation settings to existing (default is empty) list
   body_settings.get( "Earth" ).shape_deformation_settings.append(
       environment_setup.shape_deformation.degree_two_basic_solid_body_tidal(
           tide_raising_bodies, h2_love_number, l2_shida_number ) )


     )doc" );

    m.def( "iers_2010_solid_body_tidal",
           &tss::iers2010TidalBodyShapeDeformation,
           R"doc(

 Function for creating full IERS 2010 solid Earth tide shape deformation model

 Function for creating full IERS 2010 solid Earth tide shape deformation model, computing the tidal shape variation due to the full model defined in
 Section 7.1.1 of the IERS Conventions 2010, implementing Eqs. (7.5)-(7.13), including all terms from Tables 7.3a and 7.3b. At present, none of the input parameters of the model can be varied.

 Returns
 -------
 BodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` defining the IERS 2010 settings

     )doc" );

    m.def( "pole_tidal", &tss::poleTideBodyShapeDeformation, R"doc(

Function for creating IERS 2010 pole tide deformation model

Function for creating IERS 2010  2010 pole tide deformation model, computing the shape variation as a result of polar-motion induced
rotational variation, as described ion Section 7.1.4 of the IERS Conventions 2010

Returns
-------
BodyDeformationSettings
    Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived class

)doc" );

    using OceanPolePolicy = tudat::basic_astrodynamics::OceanPoleTideLatitudePolicy;
    py::enum_< OceanPolePolicy >( m, "OceanPoleTideLatitudePolicy" )
            .value( "require_grid_coverage", OceanPolePolicy::requireGridCoverage )
            .value( "clamp_to_nearest_grid_latitude_approximation", OceanPolePolicy::clampToNearestGridLatitudeApproximation );
    m.def( "ocean_pole_tidal",
           &tss::oceanPoleTideBodyShapeDeformation,
           py::arg( "grid_file" ),
           py::arg( "latitude_policy" ) = OceanPolePolicy::requireGridCoverage,
           R"doc(
Create conventional ocean pole-tide loading (IERS 7.1.5, Eq. 29), with the 2018 secular pole and IERS 2010 constants.
Supply the decompressed IERS CM grid as a plain-text file; it is loaded only at model construction.
The official grid already includes its geocentre correction. Source origin declarations are retained without conversion.
Bilinear interpolation is periodic in longitude for a complete grid. The default rejects latitudes outside the supplied
coverage (the official grid ends at +/-89.75 degrees). The explicitly named clamp approximation uses the nearest grid latitude.
Cartesian output uses the geocentric radial/north/east basis; station ENU uses its actual topocentric basis.
Standalone epochs are TT seconds since J2000; the existing environment currently approximates its TDB argument as TT.
Requires Earth's GCRS/ITRS rotation model and its polar motion data. Grid coefficients are dimensionless; output is metres.
)doc" );
    m.def( "ocean_pole_tidal_from_data",
           &tss::oceanPoleTideBodyShapeDeformationFromData,
           py::arg( "grid" ),
           py::arg( "latitude_policy" ) = OceanPolePolicy::requireGridCoverage,
           "Create the same conventional loading model from a preloaded data_input.environment_data.ocean_loading.OceanPoleTideGrid." );

    py::enum_< tudat::basic_astrodynamics::OceanTideSynthesis >( m, "OceanTideSynthesis" )
            .value( "iers_hardisp", tudat::basic_astrodynamics::OceanTideSynthesis::iersHardisp )
            .value( "legacy_eleven_constituents", tudat::basic_astrodynamics::OceanTideSynthesis::legacyElevenConstituents );

    m.def( "ocean_tidal",
           &tss::oceanTideBodyShapeDeformation,
           py::arg( "blq_files" ),
           py::arg( "synthesis" ) = tudat::basic_astrodynamics::OceanTideSynthesis::iersHardisp,
           R"doc(

Function for creating IERS 2010 ocean tide deformation model

Creates the IERS ocean tidal loading model. The default ``iers_hardisp`` expands the eleven BLQ coefficients into 342 constituents
using the corrected IERS HARDISP admittance interpolation. ``legacy_eleven_constituents`` explicitly retains the earlier eleven-term
synthesis for reproducibility. This changes the previous default behavior.

BLQ files can be obtained from the `Ocean Loading Provider <https://barre.oso.chalmers.se/loading/l.php>`_. Files are read when the
model is constructed, and source headers, model names and CM/CE/unspecified origin declarations are retained. An origin declaration
describes the supplied loading coefficients; this model does not silently convert their reference origin.

Evaluation uses TT seconds since J2000 and corresponding UTC astronomical phases. Standalone model calls require TT;
existing environment evaluation currently passes its TDB time as TT, retaining the millisecond TDB-TT approximation. The local physical components use the spherical
radial/west/south basis; Cartesian output uses geocentric station latitude, including for an oblate station shape.

Parameters
----------
blq_files : list[ string ]
    List of BLQ files defining loading amplitudes and Greenwich phase lags for the requested stations.
synthesis : OceanTideSynthesis, default=OceanTideSynthesis.iers_hardisp
    Expanded IERS synthesis or the explicit legacy eleven-constituent convention.

Returns
-------
BodyDeformationSettings
    Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived class

)doc" );

    m.def( "ocean_tidal_from_data",
           &tss::oceanTideBodyShapeDeformationFromData,
           py::arg( "data" ),
           py::arg( "synthesis" ) = tudat::basic_astrodynamics::OceanTideSynthesis::iersHardisp,
           R"doc(
Create ocean tidal loading settings from validated, preloaded BLQ records.

Parameters
----------
data : tudatpy.data_input.environment_data.ocean_loading.OceanLoadingData
    Provider coefficients read by ``read_ocean_loading_data`` or ``read_ocean_loading_data_text``.
synthesis : OceanTideSynthesis, default=OceanTideSynthesis.iers_hardisp
    Expanded IERS 342-constituent synthesis or the explicit legacy convention.

Returns
-------
BodyDeformationSettings
    Settings consuming the supplied data without additional BLQ file access.
)doc" );
}

}  // namespace shape_deformation
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
