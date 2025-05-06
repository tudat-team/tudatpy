/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "tudat/interface/spice.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/basic_astro.h>

#include "expose_spice.h"

namespace py = pybind11;
namespace tsi = tudat::spice_interface;
namespace tba = tudat::basic_astrodynamics;

namespace tudat
{

namespace spice_interface
{

void loadStandardDepracatedSpiceKernels(
        const std::vector< std::string > alternativeEphemerisKernels )
{
    std::string kernelPath = paths::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelPath + "/pck00010.tpc" );
    loadSpiceKernelInTudat( kernelPath + "/gm_de431.tpc" );

    if( alternativeEphemerisKernels.size( ) == 0 )
    {
        loadSpiceKernelInTudat( kernelPath + "/tudat_merged_spk_kernel.bsp" );
    }
    else
    {
        for( unsigned int i = 0; i < alternativeEphemerisKernels.size( ); i++ )
        {
            loadSpiceKernelInTudat( alternativeEphemerisKernels.at( i ) );
        }
    }
    loadSpiceKernelInTudat( kernelPath + "/naif0012.tls" );
}

}  // namespace spice_interface

}  // namespace tudat

namespace tudatpy
{
namespace interface
{
namespace spice
{

void expose_spice( py::module &m )
{
    // time related
    m.def( "convert_julian_date_to_ephemeris_time",
           &tudat::spice_interface::convertJulianDateToEphemerisTime,
           py::arg( "julian_date" ),
           R"doc(

 Convert a Julian date to ephemeris time (equivalent to TDB in Spice).

 Function to convert a Julian date to ephemeris time, which is
 equivalent to barycentric dynamical time. A leap second kernel
 must have been loaded to use this function.


 Parameters
 ----------
 julian_date : int
     Julian date that is to be converted to ephemeris time.
 Returns
 -------
 ephemeris_time : float    Julian date calculated from ephemeris time.






     )doc" );

    //  m.def("jd2tdb",
    //  m.attr("convert_julian_date_to_ephemeris_time"));

    m.def( "convert_ephemeris_time_to_julian_date",
           &tudat::spice_interface::convertEphemerisTimeToJulianDate,
           py::arg( "ephemeris_time" ),
           R"doc(

 Convert ephemeris time (equivalent to TDB) to a Julian date.

 Function to convert ephemeris time, which is nearly equal to
 barycentric dynamical time, to the Julian date. A leap second
 kernel must have been loaded to use this function.


 Parameters
 ----------
 ephemeris_time : float
     Ephemeris time that is to be converted to Julian date.
 Returns
 -------
 julian_date : float    Julian date calculated from ephemeris time.






     )doc" );

    //  m.def("tdb2jd",
    //  m.attr("convert_ephemeris_time_to_julian_date"));

    m.def( "convert_date_string_to_ephemeris_time",
           &tudat::spice_interface::convertDateStringToEphemerisTime,
           py::arg( "date_string" ),
           R"doc(

 Converts a date string to ephemeris time.

 Function to convert a date string, for instance
 1988 June 13, 3:29:48 to ephemeris time, wrapper for `str2et_c`_
 spice function.

 .. _`str2et_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/str2et_c.html

 Parameters
 ----------
 date_string : str
     String representing the date. See documentation of spice
     function `str2et_c` for details on supported formats.

 Returns
 -------
 ephemeris_time : str    Ephemeris time corresponding to given date_string.






     )doc" );

    //  m.def("dstr2jd",
    //  m.attr("convert_date_string_to_ephemeris_time"));

    // positional state related
    m.def( "get_body_cartesian_position_at_epoch",
           &tudat::spice_interface::getBodyCartesianPositionAtEpoch,
           py::arg( "target_body_name" ),
           py::arg( "observer_body_name" ),
           py::arg( "reference_frame_name" ),
           py::arg( "aberration_corrections" ),
           py::arg( "ephemeris_time" ),
           R"doc(

 Get Cartesian position of a body, as observed from another body.

 This function returns the position of a body, relative to another
 body, in a frame specified by the user. Corrections for light-time
 correction and stellar aberration can be applied to obtain the
 state of one of the bodies, as observed from the other. Wrapper
 for `spkpos_c`_ spice function.

 .. _`spkpos_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkpos_c.html

 Parameters
 ----------
 target_body_name : str
     Name of the body of which the state is to be obtained. A kernel
     with the ephemeris of this body must have been loaded. The
     string must be a spice-recognized name or ID.

 observer_body_name : str
     Name of the body relative to which the state is to be obtained.
     A kernel with the ephemeris of this body must have been loaded.
     The string must be a spice-recognized name or ID.

 reference_frame_name : str
     The spice-recognized name of the reference frame in which the
     state is to be returned. Spice kernel(s) required to perform
     the necessary conversion from the states of the target and
     observer bodies to this frame need to have been loaded.

 aberration_corrections : str
     Setting for correction for setting corrections. See Spice
     documentation for extended discussion.
     Short summary:

     - NONE: none
     - LT: light time corrected (one iteration for calculation)
     - CN: light time corrected (multiple iterations, max 3) for calculation,
     - S: Stellar aberration corrected.
     - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").

 ephemeris_time : float
     Observation time (or transmission time of observed light, see description
     of aberrationCorrections).






     )doc" );

    m.def( "get_body_cartesian_state_at_epoch",
           &tudat::spice_interface::getBodyCartesianStateAtEpoch,
           py::arg( "target_body_name" ),
           py::arg( "observer_body_name" ),
           py::arg( "reference_frame_name" ),
           py::arg( "aberration_corrections" ),
           py::arg( "ephemeris_time" ),
           R"doc(

 Get Cartesian state of a body, as observed from another body.

 This function returns the state of a body, relative to another
 body, in a frame specified by the user. Corrections for light-time
 correction and stellar aberration can be applied to obtain the
 state of one of the bodies, as observed from the other. Wrapper
 for `spkezr_c`_ spice function.

 .. _`spkezr_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html

 Parameters
 ----------
 target_body_name : str
     Name of the body of which the state is to be obtained. A kernel
     with the ephemeris of this body must have been loaded. The
     string must be a spice-recognized name or ID.

 observer_body_name : str
     Name of the body relative to which the state is to be obtained.
     A kernel with the ephemeris of this body must have been loaded.
     The string must be a spice-recognized name or ID.

 reference_frame_name : str
     The spice-recognized name of the reference frame in which the
     state is to be returned. Spice kernel(s) required to perform
     the necessary conversion from the states of the target and
     observer bodies to this frame need to have been loaded.

 aberration_corrections : str
     Setting for correction for setting corrections. See Spice
     documentation for extended discussion.
     Short summary:

     - NONE: none
     - LT: light time corrected (one iteration for calculation)
     - CN: light time corrected (multiple iterations, max 3) for calculation
     - S: Stellar aberration corrected.
     - XLT and XCN: can be provided to make the ephemeris time input argument the transmission time, instead of reception time. Arguments can be combined (i.e."LT+S" or "XCN+S").

 ephemeris_time : float
     Observation time (or transmission time of observed light, see description
     of aberrationCorrections).

 Returns
 -------
 cartesian_state_vector : numpy.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).






     )doc" );

    m.def( "get_cartesian_state_from_tle_at_epoch",
           &tudat::spice_interface::getCartesianStateFromTleAtEpoch,
           py::arg( "epoch" ),
           py::arg( "tle" ),
           R"doc(

 Get Cartesian state of a satellite from its two-line element set at a specified epoch.

 This function retrieves the state of a satellite at a certain epoch
 by propagating the SGP or SDP models (near-Earth resp. deep space)
 with the given two-line elements (TLE). This function serves as a
 wrapper for the `ev2lin` function in CSpice.


 Parameters
 ----------
 epoch : float
     Time in seconds since J2000 at which the state is to be retrieved.
 tle : :class:`~tudatpy.kernel.astro.ephemerides.Tle`
     Shared pointer to a Tle object containing the SGP/SDP model parameters as derived from the element set.
 Returns
 -------
 cartesian_state_vector : numpy.ndarray[6,]    Cartesian state vector (x,y,z, position+velocity).






     )doc" );

    // rotational state related
    m.def( "compute_rotation_matrix_between_frames",
           &tudat::spice_interface::computeRotationMatrixBetweenFrames,
           py::arg( "original_frame" ),
           py::arg( "new_frame" ),
           py::arg( "ephemeris_time" ),
           R"doc(

 Computes rotation matrix between two frames.

 This function computes the rotation matrix
 between two frames at a given time instant. Kernels defining the
 two frames, as well as any required intermediate frames, at the
 requested time must have been loaded. Wrapper for `pxform_c`_ spice function.

 .. _`pxform_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/pxform_c.html

 Parameters
 ----------
 original_frame
     Reference frame from which the rotation is made.
 new_frame
     Reference frame to which the rotation is made.
 ephemeris_time
     Value of ephemeris time at which rotation is to be determined.
 Returns
 -------
 Rotation matrix from original to new frame at given time.






     )doc" );

    //   m.def("compute_rotation_quaternion_between_frames",
    //         &tudat::spice_interface::computeRotationQuaternionBetweenFrames,
    //         py::arg("original_frame"),
    //         py::arg("new_frame"),
    //         py::arg("ephemeris_time"),
    //         get_docstring("compute_rotation_quaternion_between_frames").c_str());

    m.def( "compute_rotation_matrix_derivative_between_frames",
           &tudat::spice_interface::computeRotationMatrixDerivativeBetweenFrames,
           py::arg( "original_frame" ),
           py::arg( "new_frame" ),
           py::arg( "ephemeris_time" ),
           R"doc(

 Computes time derivative of rotation matrix between two frames.

 This function computes the derivative of the rotation matrix
 between two frames at a given time instant. kernels defining the
 two frames, as well as any required intermediate frames, at the
 requested time must have been loaded. Wrapper for (part of) `sxform_c`_ spice function.

 .. _`sxform_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sxform_c.html

 Parameters
 ----------
 original_frame
     Reference frame from which the rotation is made.
 new_frame
     Reference frame to which the rotation is made.
 ephemeris_time
     Value of ephemeris time at which rotation is to be determined.
 Returns
 -------
 Time derivative of rotation matrix from original to new frame at given time.






     )doc" );

    m.def( "get_angular_velocity_vector_of_frame_in_original_frame",
           &tudat::spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame,
           py::arg( "original_frame" ),
           py::arg( "new_frame" ),
           py::arg( "ephemeris_time" ),
           R"doc(

 Computes the angular velocity of one frame w.r.t. to another frame.

 Computes the angular velocity of one frame w.r.t. to another frame.
 at a given time instant. kernels defining the two frames, as well
 as any required intermediate frames, at the requested time must
 have been loaded. Wrapper for `xf2rav_c`_ spice function (utilizing `sxform_c`_).

 .. _`xf2rav_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/xf2rav_c.html
 .. _`sxform_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sxform_c.html


 Parameters
 ----------
 original_frame
     Reference frame from which the rotation is made.
 new_frame
     Reference frame to which the rotation is made.
 ephemeris_time
     Value of ephemeris time at which rotation is to be determined.
 Returns
 -------
 Angular velocity of newFrame w.r.t. originalFrame, expressed in originalFrame.






     )doc" );

    m.def( "compute_rotation_quaternion_and_rotation_matrix_"
           "derivative_between_frames",
           &tudat::spice_interface::
                   computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames,
           py::arg( "original_frame" ),
           py::arg( "new_frame" ),
           py::arg( "ephemeris_time" ),
           R"doc(No documentation found.)doc" );

    m.def( "get_body_properties",
           &tudat::spice_interface::getBodyProperties,
           py::arg( "body_name" ),
           py::arg( "property" ),
           py::arg( "max_n_val" ),
           R"doc(

 Get property of a body from Spice.

 Function to retrieve a property of a body from Spice, wraps the `bodvrd_c`_ Spice function.

 .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html

 Parameters
 ----------
 body_name
     Name of the body of which the property is to be retrieved.
 property
     Name of the property that is to be retrieved. Naming conventions can be found
     in the `bodvrd_c`_ function documentation.

     .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html

 maximum_number_of_values : int
     Number of values by which the property is expressed (i.e. 1 for
     gravitational parameter, 3 for tri-axial ellipsoid principal axes).

 Returns
 -------
 Property value(s) expressed in an STL vector of doubles.



 Notes
 -----
 Function returns values with distance unit km, not m!




     )doc" );

    m.def( "get_body_gravitational_parameter",
           &tudat::spice_interface::getBodyGravitationalParameter,
           py::arg( "body_name" ),
           R"doc(

 Get gravitational parameter of a body.

 This function retrieves the gravitational parameter of a body.
 Wraps the `bodvrd_c`_ spice function with "GM" as property type.

 .. note::

     The gravitational parameter returned is directly retrieved from the loaded SPICE kernels.
     This value does *not* necessarily correspond to the gravitational parameter used in Tudat's default models, as retrieved from e.g. :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`.
     For more information on the default gravity models, see the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models.html#gravity-field>`_.

 .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html


 Parameters
 ----------
 body
     Name of the body of which the parameter is to be retrieved.
 Returns
 -------
 Gravitational parameter of requested body.






     )doc" );

    m.def( "get_average_radius",
           &tudat::spice_interface::getAverageRadius,
           py::arg( "body_name" ),
           R"doc(

 Get the (arithmetic) mean of the three principal axes of the tri-axial ellipsoid shape.

 Returns the (arithmetic) mean of the three principal axes of the
 tri-axial ellipsoid shape of the requested body. Uses the `bodvrd_c`_ spice function with "RADII" as property type.

 .. _`bodvrd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodvrd_c.html


 Parameters
 ----------
 body
     Name of the body of which the average radius is to be retrieved.
 Returns
 -------
 Arithmetic mean of principal axes of tri-axial ellipsoid shape model of body.






     )doc" );

    m.def( "convert_body_name_to_naif_id",
           &tudat::spice_interface::convertBodyNameToNaifId,
           py::arg( "body_name" ),
           R"doc(

 Convert a body name to its NAIF identification number.

 This function converts a body name to its NAIF identification
 number. The NAIF id number is required for a number of spice
 functions, whereas the name is easily interpretable by the user.
 Wrapper for the `bods2c_c`_ function.

 .. _`bods2c_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bods2c_c.html

 Parameters
 ----------
 body_name
     Name of the body for which NAIF id is to be retrieved.
 Returns
 -------
 NAIF id number for the body with bodyName.






     )doc" );

    m.def( "convert_naif_id_to_body_name",
           &tudat::spice_interface::convertNaifIdToBodyName,
           py::arg( "naif_id" ),
           R"doc(No documentation found.)doc" );

    //        // kernel pool related
    //        m.def("get_standard_kernels",
    //              &tudat::spice_interface::getStandardSpiceKernels,
    //              get_docstring("get_standard_kernels").c_str());

    m.def( "load_standard_kernels",
           &tudat::spice_interface::loadStandardSpiceKernels,
           py::arg( "alternative_kernels" ) = std::vector< std::string >( ),  // <pybind11/stl.h>
           R"doc(

 Loads the default spice kernels shopped with tudat.

 Loads the default spice kernels shopped with tudat. The kernels that are loaded are (in order):

 - pck00010.tpc - Orientation and size/shape data for natural bodies, based mainly on  IAU Working Group on Cartographic Coordinates and Rotational Elements, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/>`__
 - inpop19a_TDB_m100_p100_spice.tpc - Masses of solar system planets and large asteroids, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`__
 - NOE-4-2020.tpc - Mars and Martian moon masses; Mars rotation model, as determined/used in NOE Martian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`__
 - NOE-5-2021.tpc - Jupiter and selected Jovian moon (Io, Europa, Ganymede, Callisto, Amalthea) masses; Jupiter rotation model, as determined/used in NOE Jovian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`__
 - NOE-6-2018-MAIN-v2.tpc - Saturn and selected Saturnian moon (Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus) masses; Saturn rotation model, as determined/used in NOE Saturnian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`__
 - codes_300ast_20100725.bsp - Ephemerides of 300 of of the largest asteroids, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/>`__
 - inpop19a_TDB_m100_p100_spice.bsp - Ephemerides of Solar system planetary system barycenters, Sun, Moon, Earth and Pluto, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`__
 - NOE-4-2020.bsp - Mars, Phobos and Deimos ephemerides (w.r.t. Martian system barycenter), as determined/used in NOE Martian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`__
 - NOE-5-2021.bsp - Jupiter, Io, Europa, Ganymede, Callisto, Amalthea ephemerides (w.r.t. Jovian system barycenter), as determined/used in NOE Jovian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`__
 - NOE-6-2018-MAIN-v2.bsp - Saturn, Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus ephemerides (w.r.t. Saturnian system barycenter), as determined/used in NOE Saturnian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`__
 - naif0012.tls - Leap second kernel, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/>`__


 Parameters
 ----------
 kernel_paths : list[str], default = None
     Optional alternative ephemeris kernels to be loaded, instead of the default ephemeris kernels. Note that using this input automatically prevents all of the above .bsp (ephemeris) kernels from loading.





     )doc" );

    m.def( "load_standard_deprecated_kernels",
           &tudat::spice_interface::loadStandardDepracatedSpiceKernels,
           py::arg( "alternative_kernels" ) = std::vector< std::string >( ),  // <pybind11/stl.h>
           R"doc(

 Loads the default spice kernels shopped with tudat.

 Loads the default spice kernels shopped with tudat. The kernels that are loaded are (in order):

 - pck00010.tpc - Orientation and size/shape data for natural bodies, based mainly on  IAU Working Group on Cartographic Coordinates and Rotational Elements, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/>`_
 - inpop19a_TDB_m100_p100_spice.tpc - Masses of solar system planets and large asteroids, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`_
 - NOE-4-2020.tpc - Mars and Martian moon masses; Mars rotation model, as determined/used in NOE Martian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`_
 - NOE-5-2021.tpc - Jupiter and selected Jovian moon (Io, Europa, Ganymede, Callisto, Amalthea) masses; Jupiter rotation model, as determined/used in NOE Jovian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`_
 - NOE-6-2018-MAIN-v2.tpc - Saturn and selected Saturnian moon (Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus) masses; Saturn rotation model, as determined/used in NOE Saturnian satellite ephemerides, obtained from `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`_
 - codes_300ast_20100725.bsp - Ephemerides of 300 of of the largest asteroids, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/>`_
 - inpop19a_TDB_m100_p100_spice.bsp - Ephemerides of Solar system planetary system barycenters, Sun, Moon, Earth and Pluto, as determined in INPOP19a ephemerides, obtained from `here <https://www.imcce.fr/recherche/equipes/asd/inpop/download19a>`_
 - NOE-4-2020.bsp - Mars, Phobos and Deimos ephemerides (w.r.t. Martian system barycenter), as determined/used in NOE Martian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/MARS/2020/>`_
 - NOE-5-2021.bsp - Jupiter, Io, Europa, Ganymede, Callisto, Amalthea ephemerides (w.r.t. Jovian system barycenter), as determined/used in NOE Jovian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/JUPITER/2021/>`_
 - NOE-6-2018-MAIN-v2.bsp - Saturn, Mimas, Enceladus, Tethys, Dione, Rhea, Titan, Hyperion, Iapetus ephemerides (w.r.t. Saturnian system barycenter), as determined/used in NOE Saturnian satellite ephemerides, obtained from  `here <ftp://ftp.imcce.fr/pub/ephem/satel/NOE/SATURNE/2018/>`_
 - naif0012.tls - Leap second kernel, obtained from `here <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/>`_


 Parameters
 ----------
 kernel_paths : list[str], default = None
     Optional alternative ephemeris kernels to be loaded, instead of the default ephemeris kernels. Note that using this input automatically prevents all of the above .bsp (ephemeris) kernels from loading.





     )doc" );

    m.def( "get_total_count_of_kernels_loaded",
           &tudat::spice_interface::getTotalCountOfKernelsLoaded,
           R"doc(

 Get the number of spice kernels currently loaded.

 This function returns the amount of Spice kernels that are loaded
 into the kernel pool. The same kernel can be loaded multiple times.
 Wrapper for the `ktotal_c`_ function.

 .. _`ktotal_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/ktotal_c.html

 Returns
 -------
 n_kernels : int    Number of spice kernels currently loaded.






     )doc" );

    m.def( "check_body_property_in_kernel_pool",
           &tudat::spice_interface::checkBodyPropertyInKernelPool,
           py::arg( "body_name" ),
           py::arg( "body_property" ),
           R"doc(

 Check if a certain property of a body is in the kernel pool.

 This function checks if a certain property of a body is in the
 kernel pool. These properties are defined in PCK kernels. Their
 names are given in the kernel file, typical names can be found in
 the Spice documentation. Wrapper for the `bodfnd_c`_ function.

 .. _`bodfnd_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/bodfnd_c.html


 Parameters
 ----------
 body_name
     Name of the body of which the property is to be checked.
 body_property
     Name of the property of which the presence is to be checked, not case-sensitive.
 Returns
 -------
 bool
     True if property is in pool, false if not.






     )doc" );

    m.def( "load_kernel",
           &tudat::spice_interface::loadSpiceKernelInTudat,
           py::arg( "kernel_file" ),
           R"doc(

 Loads a Spice kernel into the pool.

 This function loads a Spice kernel into the kernel pool, from which
 it can be used by the various internal spice routines. Matters
 regarding the manner in which Spice handles different kernels
 containing the same information can be found in the spice required
 reading documentation, kernel section. Wrapper for the `furnsh_c`_
 function.

 .. _`furnsh_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html


 Parameters
 ----------
 file_path : str
     Path to the spice kernel to be loaded.





     )doc" );

    m.def( "clear_kernels",
           &tudat::spice_interface::clearSpiceKernels,
           R"doc(

 Clear all loaded spice kernels.

 This function removes all Spice kernels from the kernel pool.
 Wrapper for the `kclear_c`_ function.

 .. _`kclear_c`: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/kclear_c.html

 Returns
 -------
 None
     None






     )doc" );

    m.def( "continue_after_errors",
           &tudat::spice_interface::toggleErrorReturn,
           R"doc(No documentation found.)doc" );

    m.def( "suppress_error_output",
           &tudat::spice_interface::suppressErrorOutput,
           R"doc(No documentation found.)doc" );

    //      py::class_<tudat::ephemerides::SpiceEphemeris,
    //            std::shared_ptr<tudat::ephemerides::SpiceEphemeris>>(m,
    //            "SpiceEphemeris",
    //                                    get_docstring("SpiceEphemeris").c_str())
    //            .def(py::init<
    //                  const std::string &,
    //                  const std::string &,
    //                  const bool,
    //                  const bool,
    //                  const bool,
    //                  const std::string &,
    //                  const double>(),
    //                  py::arg("target_body_name"),
    //                  py::arg("observer_body_name"),
    //                  py::arg("correct_for_stellar_aberration") =
    //                  false,
    //                  py::arg("correct_for_light_time_aberration")
    //                  = true,
    //                  py::arg("converge_light_time_aberration") =
    //                  false, py::arg("reference_frame_name") =
    //                  "ECLIPJ2000",
    //                  py::arg("reference_julian_day") =
    //                  tba::JULIAN_DAY_ON_J2000,
    //                  get_docstring("SpiceEphemeris.ctor").c_str())
    //            .def("get_cartesian_state",
    //                  &tudat::ephemerides::SpiceEphemeris::getCartesianState,
    //                  py::arg("seconds_since_epoch"),
    //                  get_docstring("SpiceEphemeris.get_cartesian_state").c_str());
};

}  // namespace spice
}  // namespace interface
}  // namespace tudatpy
