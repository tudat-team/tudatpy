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


Eigen::Matrix3d getRotationFromJ2000ToEclipJ2000( );

Eigen::Matrix3d getRotationFromEclipJ2000ToJ2000( );

/**
 * @brief Sets the default error handling action of the CSPICE toolkit to "RETURN".
 *
 * This function uses the CSPICE function `erract_c` to set the error handling mode
 * to "RETURN". In this mode, CSPICE routines do not abort when an error is detected.
 * Instead, they output error messages (if enabled) and immediately return upon entry.
 *
 * This mode allows for safe handling of errors in programs where continued execution
 * without program termination is desired. It is especially useful when multiple
 * CSPICE routines are called sequentially, as it prevents unexpected crashes and
 * ensures error messages are retained for later retrieval.
 *
 * @note The CSPICE toolkit provides several error handling modes. By using this
 * function, the error handling is explicitly set to "RETURN". Ensure this behavior
 * aligns with the error handling strategy of your program.
 *
 * Example Usage:
 * @code
 * toggleErrorReturn();
 * // Call CSPICE routines here. Errors will cause immediate returns instead of aborting.
 * @endcode
 *
 * @see erract_c
 */
void toggleErrorReturn( );

/**
 * @brief Sets the CSPICE error output device to "ABORT".
 *
 * This function uses the CSPICE function `errdev_c` to configure the error output device
 * such that errors cause the program to terminate execution immediately. When set to 
 * "ABORT", CSPICE routines will output error messages (if enabled) and then stop execution.
 *
 * This setting is useful for non-interactive programs or applications where errors should
 * result in immediate termination to avoid undefined behavior or further erroneous operations.
 *
 * @note The CSPICE toolkit provides multiple error output options. By using this function,
 * the error output behavior is explicitly set to "ABORT".
 *
 * Example Usage:
 * @code
 * toggleErrorAbort();
 * // CSPICE routines will now terminate the program if errors occur.
 * @endcode
 *
 * @see errdev_c
 */
void toggleErrorAbort( );

/**
 * @brief Suppresses all CSPICE error output.
 *
 * This function uses the CSPICE function `errdev_c` to configure the error output device
 * to "NULL". When set to "NULL", no error messages will be output, effectively silencing
 * all CSPICE error reporting.
 *
 * This setting is useful in scenarios where error messages are not desired, such as when
 * testing or debugging, or when errors are handled programmatically without needing
 * console or file outputs.
 *
 * @warning Suppressing error output can make debugging difficult, as errors will not be
 * displayed. Use this setting with caution.
 *
 * Example Usage:
 * @code
 * suppressErrorOutput();
 * // CSPICE routines will no longer output error messages.
 * @endcode
 *
 * @see errdev_c
 */
void suppressErrorOutput( );

/**
 * @brief Retrieves the current CSPICE long error message if an error condition exists.
 *
 * This function checks for an active error condition using the CSPICE function `failed_c`.
 * If an error is detected (i.e., `failed_c` returns `SPICETRUE`), it retrieves the detailed
 * long error message using `getmsg_c` and returns it as a standard C++ string. If no error
 * condition exists (`failed_c` returns `SPICEFALSE`), the function returns an empty string.
 *
 * The long error message provides a detailed explanation of the error, including potentially
 * specific data about the error's occurrence. This function is useful for debugging and
 * understanding error causes in CSPICE-based programs.
 *
 * @return A string containing the CSPICE long error message if an error exists, or an empty
 * string if no error condition is present.
 *
 * Example Usage:
 * @code
 * std::string errorMessage = getErrorMessage();
 * if (!errorMessage.empty()) {
 *     std::cerr << "Error: " << errorMessage << std::endl;
 * }
 * @endcode
 *
 * @note This function does not alter the CSPICE error state. If error handling needs to
 * continue, ensure that the CSPICE environment is appropriately managed after retrieving
 * the error message.
 *
 * @see failed_c
 * @see getmsg_c
 */
std::string getErrorMessage( );

/**
 * @brief Checks for a CSPICE error condition and resets the error status if one exists.
 *
 * This function uses `failed_c` to check whether a CSPICE error condition has been signaled
 * (i.e., an error was detected). If an error condition exists (`failed_c` returns `SPICETRUE`),
 * the function resets the error status using `reset_c` and returns `true`. If no error condition
 * exists (`failed_c` returns `SPICEFALSE`), the function returns `false`.
 *
 * Resetting the error status clears the error messages and re-enables the setting of long
 * error messages, allowing CSPICE routines to continue normal execution. This function is
 * particularly useful in applications that use the "RETURN" error action mode, where CSPICE
 * routines return immediately upon entry if an error condition exists.
 *
 * @return `true` if an error condition exists and has been reset, otherwise `false`.
 *
 * Example Usage:
 * @code
 * if (checkFailure()) {
 *     std::cerr << "An error occurred and was reset." << std::endl;
 * }
 * @endcode
 *
 * @note Use this function to handle and reset CSPICE error conditions in cases where
 *       continued program execution is desired. Ensure that appropriate diagnostic
 *       actions are taken before calling this function if required.
 *
 * @see failed_c
 * @see reset_c
 */
bool checkFailure( );



}// namespace spice_interface
}// namespace tudat

#endif// TUDAT_SPICE_INTERFACE_H
