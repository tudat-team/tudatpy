/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/ephemerides/compositeEphemeris.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

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
 *  \param bodies List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& transmitterCompleteEphemeris,
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& receiverCompleteEphemeris,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd )
{
    std::vector< std::shared_ptr< LightTimeCorrection > > lightTimeCorrectionFunctions;

    // Create lighttime correction functions from lightTimeCorrections
    for( unsigned int i = 0; i < lightTimeCorrections.size( ); i++ )
    {

        lightTimeCorrectionFunctions.push_back(
                    createLightTimeCorrections(
                        lightTimeCorrections[ i ], bodies, transmittingLinkEnd, receivingLinkEnd ) );
    }

    // Create light time calculator.
    return std::make_shared< LightTimeCalculator< ObservationScalarType, TimeType > >
            ( transmitterCompleteEphemeris, receiverCompleteEphemeris, lightTimeCorrectionFunctions );
}

//! Function to create a light-time calculation object
/*!
 *  Function to create a light-time calculation object from light time correction settings environment and link end
 *  identifiers.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 *  \param bodies List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{

    // Get link end state functions and create light time calculator.
    return createLightTimeCalculator< ObservationScalarType, TimeType >(
                simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, ObservationScalarType >(
                    transmittingLinkEnd, bodies ),
                simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, ObservationScalarType >(
                    receivingLinkEnd, bodies ),
                bodies, lightTimeCorrections, transmittingLinkEnd, receivingLinkEnd );
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATELIGHTTIMECALCULATOR_H
