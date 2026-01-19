/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"

#include <tudat/interface/spice.h>

namespace tsi = tudat::spice_interface;

WASM_MODULE_PATH("interface_spice")

EMSCRIPTEN_BINDINGS(tudatpy_interface_spice) {
    using namespace emscripten;

    // Time conversion functions
    function("interface_spice_convert_julian_date_to_ephemeris_time",
        &tsi::convertJulianDateToEphemerisTime);

    function("interface_spice_convert_ephemeris_time_to_julian_date",
        &tsi::convertEphemerisTimeToJulianDate);

    function("interface_spice_convert_date_string_to_ephemeris_time",
        &tsi::convertDateStringToEphemerisTime);

    function("interface_spice_get_approximate_utc_from_tdb",
        &tsi::getApproximateUtcFromTdb);

    // Position/state functions
    function("interface_spice_get_body_cartesian_position_at_epoch",
        &tsi::getBodyCartesianPositionAtEpoch);

    function("interface_spice_get_body_cartesian_state_at_epoch",
        &tsi::getBodyCartesianStateAtEpoch);

    function("interface_spice_get_cartesian_state_from_tle_at_epoch",
        &tsi::getCartesianStateFromTleAtEpoch);

    // Body properties
    function("interface_spice_get_body_gravitational_parameter",
        &tsi::getBodyGravitationalParameter);

    function("interface_spice_get_average_radius",
        &tsi::getAverageRadius);

    // Kernel management
    function("interface_spice_load_kernel",
        &tsi::loadSpiceKernelInTudat);

    function("interface_spice_clear_kernels",
        &tsi::clearSpiceKernels);

    function("interface_spice_get_total_count_of_kernels_loaded",
        &tsi::getTotalCountOfKernelsLoaded);

    // Frame conversion
    function("interface_spice_compute_rotation_quaternion_between_frames",
        &tsi::computeRotationQuaternionBetweenFrames);

    function("interface_spice_compute_rotation_matrix_derivative_between_frames",
        &tsi::computeRotationMatrixDerivativeBetweenFrames);

    function("interface_spice_get_angular_velocity_vector_of_frame_in_original_frame",
        &tsi::getAngularVelocityVectorOfFrameInOriginalFrame);

    // Body identification
    function("interface_spice_get_body_properties",
        &tsi::getBodyProperties);

    function("interface_spice_check_body_property_in_kernel_pool",
        &tsi::checkBodyPropertyInKernelPool);

    function("interface_spice_convert_body_name_to_naif_id",
        &tsi::convertBodyNameToNaifId);

    function("interface_spice_convert_naif_id_to_body_name",
        &tsi::convertNaifIdToBodyName);

    // Frame rotation
    function("interface_spice_compute_rotation_matrix_between_frames",
        &tsi::computeRotationMatrixBetweenFrames);

    function("interface_spice_compute_rotation_quaternion_and_rotation_matrix_derivative_between_frames",
        &tsi::computeRotationQuaternionAndRotationMatrixDerivativeBetweenFrames);

    // Error handling
    function("interface_spice_continue_after_errors",
        &tsi::toggleErrorReturn);

    function("interface_spice_suppress_error_output",
        &tsi::suppressErrorOutput);
}

#endif
