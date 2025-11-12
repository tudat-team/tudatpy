/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/differencedObservationPartial.h"
namespace tudat
{

namespace observation_partials
{


std::pair< observation_models::LinkEndType, observation_models::LinkEndType > getDefaultDifferencedReferenceLinkEndTypes(
        const observation_models::LinkEndType& undifferencedReferenceLinkEndType )
{
    return std::make_pair( undifferencedReferenceLinkEndType, undifferencedReferenceLinkEndType );
}

std::pair< observation_models::LinkEndType, observation_models::LinkEndType > getDifferencedTimeOfArrivalDifferencedReferenceLinkEndTypes(
        const observation_models::LinkEndType& undifferencedReferenceLinkEndType )
{
    if( undifferencedReferenceLinkEndType != observation_models::receiver )
    {
        throw std::runtime_error( "Error when getting differenced reference linke ends for differenced time of arrival, input is not supported" );
    }
    return std::make_pair( observation_models::receiver, observation_models::transmitter );
}


void DifferencedObservablePartialScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                  const std::vector< double >& times,
                                                  const observation_models::LinkEndType fixedLinkEnd,
                                                  const Eigen::VectorXd currentObservation )
{
    try
    {
        if( customCheckFunction_ != nullptr )
        {
            customCheckFunction_( fixedLinkEnd );
        }

        std::pair< observation_models::LinkEndType, observation_models::LinkEndType > differencedReferenceLinkEndTypes =
                undifferencedToDifferencedReferenceLinkEndType_( fixedLinkEnd );

        firstPartialScaling_->update( utilities::getVectorEntries( linkEndStates, firstIndices_ ),
                                      utilities::getVectorEntries( times, firstIndices_ ),
                                      differencedReferenceLinkEndTypes.first,
                                      Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );

        secondPartialScaling_->update( utilities::getVectorEntries( linkEndStates, secondIndices_ ),
                                       utilities::getVectorEntries( times, secondIndices_ ),
                                       differencedReferenceLinkEndTypes.second,
                                       Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );
    }
    catch( const std::exception& caughtException )
    {
        std::string exceptionText = std::string( caughtException.what( ) );
        throw std::runtime_error( "Error when computing differenced observation partial scaling, error: " + exceptionText );
    }
}

}  // namespace observation_partials

}  // namespace tudat
