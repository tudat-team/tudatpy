/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INERTIATENSORPARTIALS_H
#define TUDAT_INERTIATENSORPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicSineCoefficients.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 =
        ( Eigen::Matrix3d( ) << 1.0 / 3.0, 0.0, 0.0,
          0.0, 1.0 / 3.0, 0.0,
          0.0, 0.0, -2.0 / 3.0 ).finished( );


const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 =
        ( Eigen::Matrix3d( ) << 0.0, 0.0, -1.0,
          0.0, 0.0, 0.0,
          -1.0, 0.0, 0.0 ).finished( );

const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 =
        ( Eigen::Matrix3d( ) << -2.0, 0.0, 0.0,
          0.0, 2.0, 0.0,
          0.0, 0.0, 0.0 ).finished( );

const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 =
        ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
          0.0, 0.0, -1.0,
          0.0, -1.0, 0.0 ).finished( );

const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 =
        ( Eigen::Matrix3d( ) << 0.0, -2.0, 0.0,
          -2.0, 0.0, 0.0,
          0.0, 0.0, 0.0 ).finished( );

const static Eigen::Matrix3d UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_MEAN_MOMENT =
        Eigen::Matrix3d::Identity( );

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_INERTIATENSORPARTIALS_H
