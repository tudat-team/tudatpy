/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate the partial of position of a point on or in a body wrt position of that body.
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyPosition( )
{
    return Eigen::Matrix3d::Identity( );
}

//! Function to calculate the partial of position of a point on a body wrt its body-fixed position
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyFixedPointPosition(
        const Eigen::Matrix3d& rotationMatrixToInertialFrame )
{
    return rotationMatrixToInertialFrame;
}

}

}


