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

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

namespace trf = tudat::reference_frames;

WASM_MODULE_PATH("astro_frame_conversion")

namespace {
using namespace tudatpy_wasm;

Matrix3dWrapper inertialToRswRotationMatrix(const Vector6dWrapper& state) {
    return Matrix3dWrapper(trf::getInertialToRswSatelliteCenteredFrameRotationMatrix(state.eigen()));
}

Matrix3dWrapper rswToInertialRotationMatrix(const Vector6dWrapper& state) {
    return Matrix3dWrapper(trf::getRswSatelliteCenteredToInertialFrameRotationMatrix(state.eigen()));
}

Matrix3dWrapper inertialToTnwRotationMatrix(const Vector6dWrapper& state, bool nAxisPointsAway) {
    return Matrix3dWrapper(trf::getInertialToTnwRotation(state.eigen(), nAxisPointsAway));
}

Matrix3dWrapper tnwToInertialRotationMatrix(const Vector6dWrapper& state, bool nAxisPointsAway) {
    return Matrix3dWrapper(trf::getTnwToInertialRotation(state.eigen(), nAxisPointsAway));
}

Matrix3dWrapper inertialToBodyFixedRotationMatrix(double poleDeclination, double poleRightAscension, double primeMeridianLongitude) {
    return Matrix3dWrapper(trf::getInertialToPlanetocentricFrameTransformationMatrix(
        poleDeclination, poleRightAscension, primeMeridianLongitude));
}

Matrix3dWrapper bodyFixedToInertialRotationMatrix(double poleDeclination, double poleRightAscension, double primeMeridianLongitude) {
    return Matrix3dWrapper(trf::getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        poleDeclination, poleRightAscension, primeMeridianLongitude));
}

}

EMSCRIPTEN_BINDINGS(tudatpy_astro_frame_conversion) {
    using namespace emscripten;

    function("astro_frame_conversion_inertial_to_rsw_rotation_matrix", &inertialToRswRotationMatrix);
    function("astro_frame_conversion_rsw_to_inertial_rotation_matrix", &rswToInertialRotationMatrix);
    function("astro_frame_conversion_inertial_to_tnw_rotation_matrix", &inertialToTnwRotationMatrix);
    function("astro_frame_conversion_tnw_to_inertial_rotation_matrix", &tnwToInertialRotationMatrix);
    function("astro_frame_conversion_inertial_to_body_fixed_rotation_matrix", &inertialToBodyFixedRotationMatrix);
    function("astro_frame_conversion_body_fixed_to_inertial_rotation_matrix", &bodyFixedToInertialRotationMatrix);
}

#endif
