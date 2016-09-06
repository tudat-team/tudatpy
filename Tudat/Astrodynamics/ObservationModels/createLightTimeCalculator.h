/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATELIGHTTIMECALCULATOR_H
#define TUDAT_CREATELIGHTTIMECALCULATOR_H

#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{
namespace observation_models
{

//! Function to create a light-time calculation object
/*!
 *  Function to create a light-time calculation object from light time correction settings environment and link end
 *  state functions.
 *  \param transmitterCompleteEphemeris Function returning the transmitter Cartesian state as a function of time.
 *  \param receiverCompleteEphemeris Function returning the receiver Cartesian state as a function of time.
 *  \param bodyMap List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 */
template< typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > >
createLightTimeCalculator(
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >& transmitterCompleteEphemeris,
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >& receiverCompleteEphemeris,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd )
{
    std::vector< LightTimeCorrectionFunction > lightTimeCorrectionFunctions;

    // Create lighttime correction functions from lightTimeCorrections
    for( unsigned int i = 0; i < lightTimeCorrections.size( ); i++ )
    {
        LightTimeCorrectionFunction correctionFunction =
                boost::bind( &LightTimeCorrection::calculateLightTimeCorrection, createLightTimeCorrections(
                                        lightTimeCorrections[ i ], bodyMap, transmittingLinkEnd, receivingLinkEnd ),
                             _1, _2, _3, _4 );
        lightTimeCorrectionFunctions.push_back(  correctionFunction );
    }

    // Create light time calculator.
    return boost::make_shared< LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > >
            ( transmitterCompleteEphemeris, receiverCompleteEphemeris, lightTimeCorrectionFunctions );
}

//! Function to create a light-time calculation object
/*!
 *  Function to create a light-time calculation object from light time correction settings environment and link end
 *  identifiers.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 *  \param bodyMap List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 */
template< typename ObservationScalarType = double, typename TimeType = double,
          typename StateScalarType = ObservationScalarType >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > >
createLightTimeCalculator(
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections )
{
    // Check input consistency
   if( transmittingLinkEnd.second != "" || receivingLinkEnd.second != "" )
   {
       throw std::runtime_error( "Error, body reference points not yet supported when creating light time calculator" );
   }

   // Get link end state functions and create light time calculator.
   return createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
               boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                            bodyMap.at( transmittingLinkEnd.first ), _1 ),
               boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                            bodyMap.at( receivingLinkEnd.first ), _1 ),
               bodyMap, lightTimeCorrections, transmittingLinkEnd, receivingLinkEnd );
}

}

}
#endif // TUDAT_CREATELIGHTTIMECALCULATOR_H
