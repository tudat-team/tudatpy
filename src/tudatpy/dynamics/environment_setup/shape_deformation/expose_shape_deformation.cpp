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
    py::class_< tss::StationLoadingBodyDeformationSettings,
                std::shared_ptr< tss::StationLoadingBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m, "StationLoadingBodyDeformationSettings", R"doc(

         Class for defining model settings for station loading deformation.

         `BodyDeformationSettings` derived class for station loading deformation. Instances of this class are created by
         the :func:`station_loading` function.

      )doc" );
    m.def( "station_loading",
           &tss::stationLoadingBodyShapeDeformation,
           py::arg( "data" ),
           py::arg( "station_bindings" ) = std::map< std::string, std::string >( ),
           py::arg( "maximum_gap_seconds" ) = 86400.0,
           R"doc(

 Function for creating station loading deformation from a precomputed series.

 Function for creating a deformation model from a precomputed station loading series (for instance non-tidal
 atmospheric, oceanic or hydrological loading from the International Mass Loading Service). Each ground station is
 associated with a site of the series: the site given for its name in ``station_bindings`` or, if no bindings are given,
 the nearest site within whose validity radius it lies. The displacements are interpolated linearly in TAI, and
 converted to body-fixed Cartesian coordinates; their rate is the derivative of the interpolant. Epochs outside the
 series, or in a gap longer than the maximum gap, raise an exception.

 The series must be given in the body-fixed frame of the body, and with respect to the reference origin of the station
 coordinates. Avoiding signals that are modelled twice (for instance tidal signals also modelled by a tidal loading
 model, or seasonal loading when the station coordinates contain fitted seasonal terms) is the responsibility of the
 user.


 Parameters
 ----------
 data : ~tudatpy.data_input.environment_data.station_loading.StationLoadingData
     Station loading series.
 station_bindings : dict[str, str], default = {}
     Site identifier of the series for each ground station name. The station must lie within the validity radius of
     the site.
 maximum_gap_seconds : float, default = 86400
     Maximum interval between two samples over which the series is interpolated, in seconds.

 Returns
 -------
 StationLoadingBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.StationLoadingBodyDeformationSettings` class.


 Examples
 --------
 .. code-block:: python

     from tudatpy.data_input.environment_data import station_loading
     from tudatpy.dynamics.environment_setup import shape_deformation

     data = station_loading.read_ephedisp_station_loading_data( "provider.eph" )
     body_settings.get( "Earth" ).shape_deformation_settings.append(
         shape_deformation.station_loading( data, station_bindings={ "MyStation": "weuc" } ) )

 )doc" );
    py::class_< ba::AtmosphericTideEvaluationSettings >( m, "AtmosphericTideEvaluationSettings", R"doc(

Settings for the evaluation of an atmospheric S1/S2 tidal loading model (see :func:`atmospheric_tide_from_grid`).

Conventions that are not set are taken from the metadata of the data. For a grid identified as the Ray and Ponte (2003)
model with Farrell Green functions, conventions that are set in neither default to geocentric. For other data, they
must be set in one of the two.

      )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "lookup_latitude_convention", &ba::AtmosphericTideEvaluationSettings::lookupLatitudeConvention, R"doc(

         Latitude convention with which stations are located in the grid.

         :type: ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideLatitudeConvention
      )doc" )
            .def_readwrite( "displacement_basis", &ba::AtmosphericTideEvaluationSettings::displacementBasis, R"doc(

         Local basis of the displacement components of the coefficients.

         :type: ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideDisplacementBasis
      )doc" )
            .def_readwrite( "convert_to_centre_of_mass", &ba::AtmosphericTideEvaluationSettings::convertToCentreOfMass, R"doc(

         Whether displacements given with respect to the centre of the solid Earth are converted to the centre of mass, by adding the translation of Table 7.6 of the IERS Conventions 2010 (default false). Only possible for the identified Ray and Ponte (2003) product.

         :type: bool
      )doc" );
    py::class_< tss::AtmosphericTideBodyDeformationSettings,
                std::shared_ptr< tss::AtmosphericTideBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m, "AtmosphericTideBodyDeformationSettings", R"doc(

         Class for defining model settings for atmospheric S1/S2 tidal loading deformation.

         `BodyDeformationSettings` derived class for atmospheric tidal loading. Instances of this class are created by the
         ``atmospheric_tide_*`` functions of this module.

      )doc" );
    m.def( "atmospheric_tide_from_grid",
           &tss::atmosphericTideBodyShapeDeformationFromGrid,
           py::arg( "grid" ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           R"doc(

 Function for creating atmospheric S1/S2 tidal loading deformation from a grid.

 Function for creating the atmospheric tidal loading deformation model of Section 7.1.3 of the IERS Conventions 2010,
 using coefficients interpolated (bilinearly) from a global grid at the station location. The diurnal (S1) and
 semi-diurnal (S2) terms are evaluated with UT1 (Eqs. 7.19a/b), using the time converter of the Earth rotation model.
 The rate of the displacement is computed analytically. The duplicated 0/360 degree longitude columns and the polar
 rows of the grid are averaged (the differences are below 1 micrometre for the Ray and Ponte grid); a grid whose
 duplicated values differ by more than 10 micrometres is rejected.

 The displacements are given with respect to the reference origin of the grid, unless the conversion to the centre of
 mass is selected in the evaluation settings.


 Parameters
 ----------
 grid : ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideGrid
     Grid of loading coefficients (see :func:`~tudatpy.data_input.environment_data.atmospheric_loading.read_atmospheric_tide_grid`).
 evaluation_settings : AtmosphericTideEvaluationSettings, default = AtmosphericTideEvaluationSettings()
     Settings for the evaluation of the model.

 Returns
 -------
 AtmosphericTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideBodyDeformationSettings` class.

 )doc" );
    m.def( "atmospheric_tide_from_stations",
           &tss::atmosphericTideBodyShapeDeformationFromStations,
           py::arg( "stations" ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           R"doc(

 Function for creating atmospheric S1/S2 tidal loading deformation from station-specific coefficients.

 Function for creating the atmospheric tidal loading deformation model (see :func:`atmospheric_tide_from_grid`), using
 precomputed coefficients for each station. A ground station is associated with the coefficients whose identifier is
 equal to the station name. The displacement basis must be set in the data or in the evaluation settings.


 Parameters
 ----------
 stations : ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideStationData
     Station-specific loading coefficients (see :func:`~tudatpy.data_input.environment_data.atmospheric_loading.read_atmospheric_tide_station_data`).
 evaluation_settings : AtmosphericTideEvaluationSettings, default = AtmosphericTideEvaluationSettings()
     Settings for the evaluation of the model.

 Returns
 -------
 AtmosphericTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideBodyDeformationSettings` class.

 )doc" );
    m.def( "atmospheric_tide_from_grid_file",
           &tss::atmosphericTideBodyShapeDeformationFromGridFile,
           py::arg( "file_path" ),
           py::arg( "metadata" ) = tudat::input_output::AtmosphericTideMetadata( ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           R"doc(

 Function for creating atmospheric S1/S2 tidal loading deformation from a grid file.

 Same as :func:`atmospheric_tide_from_grid`, but reading the grid from a file (when the settings are created).


 Parameters
 ----------
 file_path : str
     Path of the grid file.
 metadata : ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideMetadata, default = AtmosphericTideMetadata()
     Declared properties of the grid (see :func:`~tudatpy.data_input.environment_data.atmospheric_loading.read_atmospheric_tide_grid`).
 evaluation_settings : AtmosphericTideEvaluationSettings, default = AtmosphericTideEvaluationSettings()
     Settings for the evaluation of the model.

 Returns
 -------
 AtmosphericTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideBodyDeformationSettings` class.

 )doc" );
    m.def( "atmospheric_tide_from_station_file",
           &tss::atmosphericTideBodyShapeDeformationFromStationFile,
           py::arg( "file_path" ),
           py::arg( "metadata" ) = tudat::input_output::AtmosphericTideMetadata( ),
           py::arg( "evaluation_settings" ) = ba::AtmosphericTideEvaluationSettings( ),
           R"doc(

 Function for creating atmospheric S1/S2 tidal loading deformation from a station coefficient file.

 Same as :func:`atmospheric_tide_from_stations`, but reading the coefficients from a file (when the settings are
 created).


 Parameters
 ----------
 file_path : str
     Path of the station coefficient file.
 metadata : ~tudatpy.data_input.environment_data.atmospheric_loading.AtmosphericTideMetadata, default = AtmosphericTideMetadata()
     Declared properties of the data (see :func:`~tudatpy.data_input.environment_data.atmospheric_loading.read_atmospheric_tide_station_data`).
 evaluation_settings : AtmosphericTideEvaluationSettings, default = AtmosphericTideEvaluationSettings()
     Settings for the evaluation of the model.

 Returns
 -------
 AtmosphericTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.AtmosphericTideBodyDeformationSettings` class.

 )doc" );
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

 The displacement includes the time-independent (permanent) part of the degree-2 tide (equivalent to steps 1 and 2 of
 the IERS DEHANTTIDEINEL software, without the step 3 correction). The nominal station coordinates must therefore be
 conventional tide-free coordinates.

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

    py::class_< tss::OceanPoleTideBodyDeformationSettings,
                std::shared_ptr< tss::OceanPoleTideBodyDeformationSettings >,
                tss::BodyDeformationSettings >( m, "OceanPoleTideBodyDeformationSettings", R"doc(

         Class for defining model settings for ocean pole tide loading deformation.

         `BodyDeformationSettings` derived class for ocean pole tide loading. Instances of this class are created by the
         :func:`ocean_pole_tidal` and :func:`ocean_pole_tidal_from_data` functions.

      )doc" );
    m.def( "ocean_pole_tidal",
           &tss::oceanPoleTideBodyShapeDeformation,
           py::arg( "grid_file" ),
           R"doc(

 Function for creating ocean pole tide loading deformation.

 Function for creating the ocean pole tide loading deformation model of Section 7.1.5 (Eq. 7.29) of the IERS
 Conventions 2010, using the 2018 secular pole and the IERS 2010 constants. The coefficients are interpolated
 (bilinearly) from the IERS grid at the station location. The official grid covers latitudes from -89.75 to 89.75
 degrees; at latitudes beyond these, the edge values are used. The official grid already includes the correction to
 the centre of mass, so no further correction is applied. The model requires the GCRS-ITRS rotation model of the Earth,
 from which the polar motion is taken.


 Parameters
 ----------
 grid_file : str
     Path of the (decompressed) IERS grid file, which is read when the settings are created. The file is available
     from the `IERS Conventions Centre <https://iers-conventions.obspm.fr/content/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_.

 Returns
 -------
 OceanPoleTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.OceanPoleTideBodyDeformationSettings` class.

 )doc" );
    m.def( "ocean_pole_tidal_from_data",
           &tss::oceanPoleTideBodyShapeDeformationFromData,
           py::arg( "grid" ),
           R"doc(

 Function for creating ocean pole tide loading deformation from a grid that has already been read.

 Same as :func:`ocean_pole_tidal`, but using a grid read with
 :func:`~tudatpy.data_input.environment_data.ocean_loading.read_ocean_pole_tide_grid`.


 Parameters
 ----------
 grid : ~tudatpy.data_input.environment_data.ocean_loading.OceanPoleTideGrid
     Ocean pole tide loading grid.

 Returns
 -------
 OceanPoleTideBodyDeformationSettings
     Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived
     :class:`~tudatpy.dynamics.environment_setup.shape_deformation.OceanPoleTideBodyDeformationSettings` class.

 )doc" );

    py::enum_< tudat::basic_astrodynamics::OceanTideSynthesis >( m, "OceanTideSynthesis", R"doc(

Enumeration of methods to compute the ocean tidal loading displacement from the eleven constituents of a BLQ file.

      )doc" )
            .value( "iers_hardisp",
                    tudat::basic_astrodynamics::OceanTideSynthesis::iersHardisp,
                    "Expansion to 342 constituents by admittance interpolation, as in the IERS HARDISP software (default)." )
            .value( "legacy_eleven_constituents",
                    tudat::basic_astrodynamics::OceanTideSynthesis::legacyElevenConstituents,
                    "Synthesis of the eleven constituents only (Eq. 7.16 of the IERS Conventions 2010), as used by Tudat "
                    "before the HARDISP synthesis was available." );

    m.def( "ocean_tidal",
           &tss::oceanTideBodyShapeDeformation,
           py::arg( "blq_files" ),
           py::arg( "synthesis" ) = tudat::basic_astrodynamics::OceanTideSynthesis::iersHardisp,
           R"doc(

Function for creating IERS 2010 ocean tide deformation model

Function for creating the IERS 2010 ocean tidal loading deformation model, as described in Section 7.1.2 of the IERS
Conventions 2010. This function requires BLQ files for the stations under consideration, which can be queried on the
`Ocean Loading Provider <https://barre.oso.chalmers.se/loading/l.php>`_. The files are read when the settings are
created. By default, the eleven constituents of the BLQ files are expanded to 342 constituents, as in the IERS HARDISP
software; the synthesis of the eleven constituents only (the method used by Tudat previously) can be selected instead.
The reference origin of the coefficients is not converted.

Parameters
----------
blq_files : list[ string ]
    List of blq files defining ocean tide amplitudes and phases for the requested stations
synthesis : OceanTideSynthesis, default = OceanTideSynthesis.iers_hardisp
    Method to compute the displacement from the BLQ constituents.

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

Function for creating IERS 2010 ocean tide deformation model from BLQ data that has already been read.

Same as :func:`ocean_tidal`, but using data read with :func:`~tudatpy.data_input.environment_data.ocean_loading.read_ocean_loading_data`.

Parameters
----------
data : ~tudatpy.data_input.environment_data.ocean_loading.OceanLoadingData
    Ocean tidal loading coefficients.
synthesis : OceanTideSynthesis, default = OceanTideSynthesis.iers_hardisp
    Method to compute the displacement from the BLQ constituents.

Returns
-------
BodyDeformationSettings
    Instance of the :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings` derived class

)doc" );
}

}  // namespace shape_deformation
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
