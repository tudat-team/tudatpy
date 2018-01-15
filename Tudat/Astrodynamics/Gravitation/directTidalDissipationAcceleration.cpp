/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace gravitation
{

//! Function to compute the acceleration acting on a satellite due to tidal deformation caused by this satellite on host planet.
Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnPlanet(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag,
        const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );
    double distanceSquared = distance * distance;
    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                 radialComponentMultiplier * relativePosition + timeLag * (
                    2.0 * ( relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) +
                    ( relativePosition.cross( planetAngularVelocityVector ) + relativeVelocity ) ) );

}

//! Function to compute the acceleration acting on a satellite due to tidal deformation caused in this satellite by host planet.
Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier,
        const double timeLag,const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double distance = relativePosition.norm( );
    double distanceSquared = distance * distance;
    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                 2.0 * radialComponentMultiplier * relativePosition + timeLag * (
                   7.0 * relativePosition.dot( relativeVelocity ) * relativePosition / distanceSquared ) );

}

//! Function to retrieve all DirectTidalDissipationAcceleration from an AccelerationMap, for specific deformed/deforming bodies
std::vector< boost::shared_ptr< DirectTidalDissipationAcceleration > > getTidalDissipationAccelerationModels(
        const basic_astrodynamics::AccelerationMap accelerationModelList, const std::string bodyBeingDeformed,
        const std::vector< std::string >& bodiesCausingDeformation )
{
    // Iterate over all bodies undergoing acceleration
    std::vector< boost::shared_ptr< DirectTidalDissipationAcceleration > > selectedDissipationModels;
    for( basic_astrodynamics::AccelerationMap::const_iterator modelIterator1 = accelerationModelList.begin( );
         modelIterator1 != accelerationModelList.end( ); modelIterator1++ )
    {
        std::string bodyUndergoingAcceleration = modelIterator1->first;

        // Iterate over all bodies exerting acceleration
        basic_astrodynamics::SingleBodyAccelerationMap singleBodyAccelerationList = modelIterator1->second;
        for( basic_astrodynamics::SingleBodyAccelerationMap::const_iterator modelIterator2 = singleBodyAccelerationList.begin( );
             modelIterator2 != singleBodyAccelerationList.end( ); modelIterator2++ )
        {
            std::string bodyExertingAcceleration = modelIterator2->first;

            // Iterate over all accelerations being exerted by bodyExertingAcceleration on bodyUndergoingAcceleration
            for( unsigned int i = 0; i < modelIterator2->second.size( ); i++ )
            {
                // Check if acceleration model is due to tidal dissipations
                boost::shared_ptr< DirectTidalDissipationAcceleration > currentDissipationAcceleration =
                    boost::dynamic_pointer_cast< DirectTidalDissipationAcceleration >( modelIterator2->second.at( i ) );
                if( currentDissipationAcceleration != NULL )
                {
                    // Check whether model correspionds to input requirements
                    if( currentDissipationAcceleration->getModelTideOnPlanet( ) &&
                            ( bodyExertingAcceleration == bodyBeingDeformed ) )
                    {
                        if( ( std::find( bodiesCausingDeformation.begin( ), bodiesCausingDeformation.end( ),
                                       bodyUndergoingAcceleration ) != bodiesCausingDeformation.end( ) ) ||
                                bodiesCausingDeformation.size( ) == 0 )
                        {
                            selectedDissipationModels.push_back( currentDissipationAcceleration );
                        }
                    }
                    else if( !currentDissipationAcceleration->getModelTideOnPlanet( ) &&
                             ( bodyUndergoingAcceleration == bodyBeingDeformed ) )
                    {
                        if( ( std::find( bodiesCausingDeformation.begin( ), bodiesCausingDeformation.end( ),
                                       bodyExertingAcceleration ) != bodiesCausingDeformation.end( ) ) ||
                                bodiesCausingDeformation.size( ) == 0 )
                        {
                            selectedDissipationModels.push_back( currentDissipationAcceleration );
                        }
                    }
                }
            }
        }
    }
    return selectedDissipationModels;
}


} // namespace gravitation

} // namespace tudat
