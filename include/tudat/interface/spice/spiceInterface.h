/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      In order to use the Spice interface, the C-Spice toolkit must be installed on your machine.
 *      It can be downloaded from http://naif.jpl.nasa.gov/naif/toolkit_C.html for all common
 *      operating systems. By placing the cspice folder in the external directory, or 1, 2 or 3
 *      folder levels above the project source directory, it will be automatically located by the
 *      cmake list. In order to use Spice with Tudat, please run the makeall file provided with
 *      Spice to compile the static library.
 *      IMPORTANT: Before being able to use it, the cspice.a file in the cspice/lib folder needs to
 *      be renamed to libcspice.a.
 *
 *      In addition, the USE_CSPICE variable needs to be set to 1 in the top-level CMakeLists.txt.
 *
 */

#ifndef TUDAT_SPICE_INTERFACE_H
#define TUDAT_SPICE_INTERFACE_H

#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/io/basicInputOutput.h"

extern "C" {
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
}

#include "tudat/astro/ephemerides/tleEphemeris.h"

namespace tudat {

namespace spice_interface {

//! @get_docstring(convert_julian_date_to_ephemeris_time)
double convertJulianDateToEphemerisTime(const double julianDate);

//! @get_docstring(convert_ephemeris_time_to_julian_date)
double convertEphemerisTimeToJulianDate(const double ephemerisTime);

//! @get_docstring(convert_date_string_to_ephemeris_time)
double convertDateStringToEphemerisTime(const std::string &dateString);

//! @get_docstring(get_body_cartesian_state_at_epoch)
Eigen::Vector6d getBodyCartesianStateAtEpoch(
    const std::string &targetBodyName, const std::string &observerBodyName,
    const std::string &referenceFrameName, const std::string &aberrationCorrections,
    const double ephemerisTime);

//! @get_docstring(get_body_cartesian_position_at_epoch)
Eigen::Vector3d getBodyCartesianPositionAtEpoch(const std::string &targetBodyName,
                                                const std::string &observerBodyName,
                                                const std::string &referenceFrameName,
                                                const std::string &aberrationCorrections,
                                                const double ephemerisTime);

//! @get_docstring(get_cartesian_state_from_tle_at_epoch)
Eigen::Vector6d getCartesianStateFromTleAtEpoch(double epoch, std::shared_ptr<ephemerides::Tle> tle);

//! @get_docstring(compute_rotation_quaternion_between_frames)
Eigen::Quaterniond computeRotationQuaternionBetweenFrames(const std::string &originalFrame,
                                                          const std::string &newFrame,
                                                          const double ephemerisTime);

Eigen::Matrix3d computeRotationMatrixBetweenFrames(const std::string &originalFrame,
                                                   const std::string &newFrame,
                                                   const double ephemerisTime);

Eigen::Matrix6d computeStateRotationMatrixBetweenFrames(const std::string &originalFrame,
                                                   const std::string &newFrame,
                                                   const double ephemerisTime);

//! @get_docstring(compute_rotation_matrix_derivative_between_frames)
Eigen::Matrix3d computeRotationMatrixDerivativeBetweenFrames(const std::string &originalFrame,
                                                             const std::string &newFrame,
                                                             const double ephemerisTime);

//! @get_docstring(get_angular_velocity_vector_of_frame_in_original_frame)
Eigen::Vector3d getAngularVelocityVectorOfFrameInOriginalFrame(const std::string &originalFrame,
                                                               const std::string &newFrame,
                                                               const double ephemerisTime);

std::pair<Eigen::Quaterniond, Eigen::Matrix3d> computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames(
    const std::string &originalFrame, const std::string &newFrame, const double ephemerisTime);

//! @get_docstring(get_body_properties)
std::vector<double> getBodyProperties(const std::string &body,
                                      const std::string &property,
                                      const int maximumNumberOfValues = 1);

//! @get_docstring(get_body_gravitational_parameter)
double getBodyGravitationalParameter(const std::string &body);

//! @get_docstring(get_average_radius)
double getAverageRadius(const std::string &body);

double getAverageEquatorialRadius( const std::string& body );

double getPolarRadius( const std::string& body );

//! @get_docstring(convert_body_name_to_naif_id)
int convertBodyNameToNaifId(const std::string &bodyName);

//! Convert a NAIF identification number to its body name.
std::string convertNaifIdToBodyName( int bodyNaifId );

//! @get_docstring(check_body_property_in_kernel_pool)
bool checkBodyPropertyInKernelPool(const std::string &bodyName, const std::string &bodyProperty);

//! @get_docstring(load_kernel)
void loadSpiceKernelInTudat(const std::string &fileName);

//! @get_docstring(get_total_count_of_kernels_loaded)
int getTotalCountOfKernelsLoaded();

//! @get_docstring(clear_kernels)
void clearSpiceKernels();

//! @get_docstring(get_standard_kernels)
std::vector<std::string> getStandardSpiceKernels(const std::vector<std::string> alternativeEphemerisKernels =
                                                     std::vector<std::string>());

//! @get_docstring(load_standard_kernels)
void loadStandardSpiceKernels(const std::vector<std::string> alternativeEphemerisKernels =
                                  std::vector<std::string>());

}// namespace spice_interface
}// namespace tudat

#endif// TUDAT_SPICE_INTERFACE_H
