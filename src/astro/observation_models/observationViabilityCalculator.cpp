/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/observationViabilityCalculator.h"

namespace tudat
{

namespace observation_models
{


double getEvaluationEpochOfViabilityBody(
    const std::vector< Eigen::Vector6d >& linkEndStates,
    const std::vector< double >& linkEndTimes,
    const std::pair< int, int > linkEndIndexPair,
    const std::function< Eigen::Vector6d( const double ) > viabilityBodyStateFunction )
{

    double firstEndEpoch = linkEndTimes.at( linkEndIndexPair.first );
    double secondEndEpoch = linkEndTimes.at( linkEndIndexPair.second );

    Eigen::Vector3d positionOfFirstEnd = linkEndStates.at( linkEndIndexPair.first  ).segment( 0, 3 );
    Eigen::Vector3d positionOfSecondEnd = linkEndStates.at( linkEndIndexPair.second ).segment( 0, 3 );

    // Get position of occulting body
    Eigen::Vector3d positionOfOccultingBodyAtFirstInstant = viabilityBodyStateFunction(
        linkEndTimes.at( linkEndIndexPair.first )).segment( 0, 3 );
    Eigen::Vector3d positionOfOccultingBodyAtSecondInstant = viabilityBodyStateFunction(
        linkEndTimes.at( linkEndIndexPair.second )).segment( 0, 3 );


    double distanceToFirstEnd = (positionOfOccultingBodyAtFirstInstant - positionOfFirstEnd).norm();
    double distanceToSecondEnd = (positionOfOccultingBodyAtSecondInstant - positionOfSecondEnd).norm();

    return firstEndEpoch * (distanceToSecondEnd / (distanceToFirstEnd + distanceToSecondEnd)) +
           secondEndEpoch * (distanceToFirstEnd / (distanceToFirstEnd + distanceToSecondEnd));
}

//! Function to check whether an observation is viable
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times, const LinkEnds& linkEnds,
        const std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > >& viabilityCalculators )
{
    bool isObservationFeasible = 1;

    if( viabilityCalculators.count( linkEnds ) > 0 )
    {
        isObservationFeasible = isObservationViable( states, times, viabilityCalculators.at( linkEnds ) );
    }

    return isObservationFeasible;
}

//! Function to check whether an observation is viable
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times,
        const std::vector< std::shared_ptr< ObservationViabilityCalculator > >& viabilityCalculators )
{
    bool isObservationFeasible = 1;

    for( unsigned int i = 0; i < viabilityCalculators.size( ); i++ )
    {
        if( viabilityCalculators.at( i )->isObservationViable( states, times ) == 0 )
        {
            isObservationFeasible = 0;
            break;
        }
    }

    return isObservationFeasible;
}

//! Function for determining whether the elevation angle at station is sufficient to allow observation
bool MinimumElevationAngleCalculator::isObservationViable(
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;

    // Iterate over all sets of entries of input vector for which elvation angle is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        double elevationAngle = ground_stations::calculateGroundStationElevationAngle(
                    pointingAngleCalculator_, linkEndStates, linkEndTimes, linkEndIndices_.at( i ) );
        // Check if elevation angle criteria is met for current link.
        if( elevationAngle < minimumElevationAngle_ )
        {
            isObservationPossible = false;
        }
    }

    return isObservationPossible;
}

double computeMinimumLinkDistanceToPoint( const Eigen::Vector3d& observingBody,
                                          const Eigen::Vector3d& transmittingBody,
                                          const Eigen::Vector3d& relativePoint )
{
    Eigen::Vector3d observerToPoint = relativePoint - observingBody;
    Eigen::Vector3d transmitterToPoint = relativePoint - transmittingBody;
    Eigen::Vector3d transmitterToObserver = observingBody - transmittingBody;
    Eigen::Vector3d observerToTransmitter = transmittingBody - observingBody;

    double transmitterAngle = transmitterToObserver.dot( transmitterToPoint );
    double observerAngle = transmitterToObserver.dot( observerToPoint );

    double minimumDistance = TUDAT_NAN;

    // Minimum distance is at one of the extremes
    if( transmitterAngle * observerAngle > 0.0 )
    {
        if( std::fabs( transmitterAngle ) < std::fabs( observerAngle ) )
        {
            minimumDistance = transmitterToPoint.norm( );
        }
        else
        {
            minimumDistance = observerToPoint.norm( );
        }
    }

    // Minimum distance is impact parameter
    else
    {
        // Compute as average value with both points as reference, to minimuze numerical errors
        minimumDistance =
            ( ( transmitterToObserver.cross( transmitterToPoint ) ).norm( ) / transmitterToObserver.norm( ) +
             ( observerToTransmitter.cross( observerToPoint ) ).norm( ) / observerToTransmitter.norm( ) ) / 2.0;
    }

    return minimumDistance;

}

double computeCosineBodyAvoidanceAngle( const Eigen::Vector3d& observingBody,
                                        const Eigen::Vector3d& transmittingBody,
                                        const Eigen::Vector3d& bodyToAvoid )
{
    return linear_algebra::computeCosineOfAngleBetweenVectors(
                bodyToAvoid - observingBody,
                transmittingBody - observingBody );
}

double computeCosineBodyAvoidanceAngle( const std::vector< Eigen::Vector6d >& linkEndStates,
                                        const std::pair< int, int > observingAndTransmittingIndex,
                                        const Eigen::Vector3d& bodyToAvoid )
{
    return computeCosineBodyAvoidanceAngle(
                linkEndStates.at( observingAndTransmittingIndex.first ).segment< 3 >( 0 ),
                linkEndStates.at( observingAndTransmittingIndex.second ).segment< 3 >( 0 ),
                bodyToAvoid );
}

//! Function for determining whether the avoidance angle to a given body at station is sufficient to allow observation.
bool BodyAvoidanceAngleCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                        const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;
    Eigen::Vector3d positionOfBodyToAvoid;
    double currentCosineOfAngle;

    // Iterate over all sets of entries of input vector for which avoidance angle is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {
        double evaluationEpoch = getEvaluationEpochOfViabilityBody( linkEndStates, linkEndTimes, linkEndIndices_.at( i ), stateFunctionOfBodyToAvoid_ );

        // Compute cosine of avoidance angles
        positionOfBodyToAvoid = stateFunctionOfBodyToAvoid_( evaluationEpoch  )
                .segment( 0, 3 );
        currentCosineOfAngle = computeCosineBodyAvoidanceAngle(
                    linkEndStates, linkEndIndices_.at( i ), positionOfBodyToAvoid );

        // Check if avoidance angle is sufficiently large
        if( currentCosineOfAngle > std::cos( bodyAvoidanceAngle_ ) )
        {
            isObservationPossible = 0;
            break;
        }
    }

    return isObservationPossible;
}

bool computeOccultation(
    const Eigen::Vector3d observer1Position,
    const Eigen::Vector3d observer2Position,
    const Eigen::Vector3d occulterPosition,
    const double radius )
{

    double observerRelativeDistance = (observer2Position - observer1Position ).norm( );
    double observer1OcculterDistance = (occulterPosition - observer1Position ).norm( );
    double observer2OcculterDistance = (occulterPosition - observer2Position ).norm( );

    double cosineBody1Angle =
       - ( observer1OcculterDistance * observer1OcculterDistance - observer2OcculterDistance * observer2OcculterDistance - observerRelativeDistance * observerRelativeDistance ) /
            ( 2.0 * observer2OcculterDistance * observerRelativeDistance );
    double cosineBody2Angle =
        - ( observer2OcculterDistance * observer2OcculterDistance - observer1OcculterDistance * observer1OcculterDistance - observerRelativeDistance * observerRelativeDistance ) /
        ( 2.0 * observer1OcculterDistance * observerRelativeDistance );

    if( cosineBody1Angle < 0.0 || cosineBody2Angle < 0.0 )
    {
        return false;
    }
    else
    {
        double distanceToTest = observer2OcculterDistance * std::sqrt( 1.0 - cosineBody1Angle * cosineBody1Angle );
        if ( distanceToTest < radius )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

//! Function for determining whether the link is occulted during the observataion.
bool OccultationCalculator::isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                 const std::vector< double >& linkEndTimes )
{
    bool isObservationPossible = 1;


    // Iterate over all sets of entries of input vector for which occultation is to be checked.
    for( unsigned int i = 0; i < linkEndIndices_.size( ); i++ )
    {

        double evaluationEpoch = getEvaluationEpochOfViabilityBody( linkEndStates, linkEndTimes, linkEndIndices_.at( i ), stateFunctionOfOccultingBody_ );

//        std::cout<<"Epoch "<<evaluationEpoch<<" "<<linkEndTimes.at( 0 )<<" "<<linkEndTimes.at( 1 )<<std::endl;
//        std::cout<<"Epoch diff "<<" "<<linkEndTimes.at( 0 ) - evaluationEpoch<<" "<<linkEndTimes.at( 1 )- evaluationEpoch<<std::endl;
//        std::cout<<"Position diff "<<( linkEndStates.at( 1 ) - stateFunctionOfOccultingBody_( linkEndTimes.at( 1 ) ) ).segment( 0, 3 ).norm( )<<std::endl;
//        // Check if observing link end is occulted by body.
//        if( mission_geometry::computeShadowFunction(
//                    linkEndStates.at( linkEndIndices_.at( i ).first ).segment( 0, 3 ), 0.0,
//                    stateFunctionOfOccultingBody_( evaluationEpoch ).segment( 0, 3 ),
//                    radiusOfOccultingBody_,
//                    linkEndStates.at( linkEndIndices_.at( i ).second ).segment( 0, 3 ) ) < 1.0E-10 )
        if( computeOccultation(
            linkEndStates.at( linkEndIndices_.at( i ).first ).segment( 0, 3 ),
            linkEndStates.at( linkEndIndices_.at( i ).second ).segment( 0, 3 ),
            stateFunctionOfOccultingBody_( evaluationEpoch ).segment( 0, 3 ),
            radiusOfOccultingBody_ ) )
        {
            isObservationPossible = 0;
            break;
        }
    }

    return isObservationPossible;
}



}

}
