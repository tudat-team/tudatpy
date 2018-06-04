/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Function to calculate the altitude of a point over a central body
//! from positions of both the point and the body (in any frame)
double getAltitudeFromNonBodyFixedPosition(
        const std::shared_ptr< BodyShapeModel > bodyShapeModel, const Eigen::Vector3d& position,
        const Eigen::Vector3d& bodyPosition, const Eigen::Quaterniond& toBodyFixedFrame )
{
    return bodyShapeModel->getAltitude( toBodyFixedFrame * ( position - bodyPosition ) );
}

//! Function to calculate the altitude of a point over a central body
//! from positions of both the point and the body (in any frame)
double getAltitudeFromNonBodyFixedPositionFunctions(
        const std::shared_ptr< BodyShapeModel > bodyShapeModel, const Eigen::Vector3d& position,
        const std::function< Eigen::Vector3d( ) > bodyPositionFunction,
        const std::function< Eigen::Quaterniond( ) > toBodyFixedFrameFunction )
{
    return getAltitudeFromNonBodyFixedPosition( bodyShapeModel, position,  bodyPositionFunction( ),
                                                toBodyFixedFrameFunction( ) );
}

} // namespace basic_astrodynamics
} // namespace tudat
