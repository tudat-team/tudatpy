/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/positionAngleAndSeparationPartial.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace observation_partials
{

std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > calculatePositionAngleAndSeparationPartialWrtLinesOfSight(
        const Eigen::Vector3d& firstLineOfSight,
        const Eigen::Vector3d& secondLineOfSight,
        const Eigen::Vector3d& northPoleDirection,
        const bool calculatePositionAngle )
{
    Eigen::Matrix< double, 2, 3 > firstUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    Eigen::Matrix< double, 2, 3 > secondUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    if( calculatePositionAngle )
    {
        const Eigen::Vector3d unnormalisedEastDirection = northPoleDirection.cross( firstLineOfSight );
        const double eastDirectionNorm = unnormalisedEastDirection.norm( );
        const Eigen::Vector3d eastDirection = unnormalisedEastDirection / eastDirectionNorm;
        const Eigen::Vector3d northDirection = firstLineOfSight.cross( eastDirection );

        const double eastProjection = secondLineOfSight.dot( eastDirection );
        const double northProjection = secondLineOfSight.dot( northDirection );
        const double positionAngleDenominator = eastProjection * eastProjection + northProjection * northProjection;

        const Eigen::Matrix3d eastDirectionPartial = ( Eigen::Matrix3d::Identity( ) - eastDirection * eastDirection.transpose( ) ) *
                linear_algebra::getCrossProductMatrix( northPoleDirection ) / eastDirectionNorm;
        const Eigen::Matrix3d northDirectionPartial = -linear_algebra::getCrossProductMatrix( eastDirection ) +
                linear_algebra::getCrossProductMatrix( firstLineOfSight ) * eastDirectionPartial;

        firstUnitVectorPartial.row( 0 ) = ( northProjection * secondLineOfSight.transpose( ) * eastDirectionPartial -
                                            eastProjection * secondLineOfSight.transpose( ) * northDirectionPartial ) /
                positionAngleDenominator;
        secondUnitVectorPartial.row( 0 ) =
                ( northProjection * eastDirection.transpose( ) - eastProjection * northDirection.transpose( ) ) / positionAngleDenominator;
    }

    const double cosineOfSeparation = firstLineOfSight.dot( secondLineOfSight );
    const double sineOfSeparation = firstLineOfSight.cross( secondLineOfSight ).norm( );
    firstUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * firstLineOfSight - secondLineOfSight ).transpose( ) / sineOfSeparation;
    secondUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * secondLineOfSight - firstLineOfSight ).transpose( ) / sineOfSeparation;

    return std::make_pair( firstUnitVectorPartial, secondUnitVectorPartial );
}

}  // namespace observation_partials

}  // namespace tudat
