/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrection.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{
namespace observation_models
{

//! Function to create an ephemeris for a reference point on a body
/*!
 *  Function to create an ephemeris for a reference point on a body, taking into account the time-variable rotation
 *  of the body and its global ephemeris
 *  \param bodyWithReferencePoint Body on which reference point is located
 *  \param bodyRotationModel Rotation model that is to be used for going from body-fixed to inertial frame
 *  \param referencePointStateFunction Function returning the state of the reference point on the body (in a body-fixed
 *  frame).
 *  \return Reference point ephemeris in global coordinates.
 */
template< typename TimeType = double, typename StateScalarType = double >
boost::shared_ptr< ephemerides::Ephemeris > createReferencePointEphemeris(
        const boost::shared_ptr< simulation_setup::Body > bodyWithReferencePoint,
        const boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel,
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > referencePointStateFunction )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    // Create list of state/rotation functions that are to be used
    std::map< int, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > > stationEphemerisVector;
    stationEphemerisVector[ 2 ] = boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris
                                               < StateScalarType, TimeType >, bodyWithReferencePoint, _1 );
    stationEphemerisVector[ 0 ] = referencePointStateFunction;

    std::map< int, boost::function< StateType( const TimeType, const StateType& ) > > stationRotationVector;
    stationRotationVector[ 1 ] =  boost::bind( &ephemerides::transformStateToGlobalFrame
                                               < StateScalarType, TimeType >, _2, _1, bodyRotationModel );

    // Create and return ephemeris
    return boost::make_shared< ephemerides::CompositeEphemeris< TimeType, StateScalarType > >(
                stationEphemerisVector, stationRotationVector, "SSB", "ECLIPJ2000" );
}

//! Function to retrieve a state function for a link end (either a body center of mass or ground station).
/*!
 *  Function to retrieve a state function for a link end (either a body center of mass or ground station).
 *  \param bodyWithLinkEnd Body on/in which link end is situated.
 *  \param linkEndId Id of link end for which state function is to be created. First: name of body, second: name of
 *  reference point (empty if center of mass is to be used
 *  \return Requested state function
 */
template< typename TimeType = double, typename StateScalarType = double >
boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > getLinkEndCompleteEphemerisFunction(
        const boost::shared_ptr< simulation_setup::Body > bodyWithLinkEnd,
        const std::pair< std::string, std::string >& linkEndId )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    boost::function< StateType( const TimeType& ) > linkEndCompleteEphemerisFunction;

    // Checking transmitter if a reference point is to be used
    if( linkEndId.second != "" )
    {
        if( bodyWithLinkEnd->getGroundStationMap( ).count( linkEndId.second ) == 0 )
        {
            std::string errorMessage = "Error when making ephemeris function for " + linkEndId.first + ", " +
                    linkEndId.second + ", station not found.";
            throw std::runtime_error( errorMessage );
        }

        // Retrieve function to calculate state of transmitter S/C
        linkEndCompleteEphemerisFunction =
                boost::bind( &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType,TimeType >,
                             createReferencePointEphemeris< TimeType, StateScalarType >(
                                 bodyWithLinkEnd, bodyWithLinkEnd->getRotationalEphemeris( ),
                                 boost::bind( &ground_stations::GroundStation::getStateInPlanetFixedFrame
                                              < StateScalarType, TimeType >,
                                              bodyWithLinkEnd->getGroundStation( linkEndId.second ), _1 ) ), _1 );

    }
    // Else, create state function for center of mass
    else
    {
        // Create function to calculate state of transmitting ground station.
        linkEndCompleteEphemerisFunction =
                boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                                        bodyWithLinkEnd, _1 );
    }
    return linkEndCompleteEphemerisFunction;
}

//! Function to create a state function of a link end, expressed in base frame.
/*!
 *  Function to create a state function of a link end, expressed in base frame.
 *  \param linkEndId Name of the reference point for which state function is to be created.
 *  \param bodyMap List of body objects that comprises the environment
 *  \return Requested state function.
 */
template< typename TimeType = double, typename StateScalarType = double >
boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > getLinkEndCompleteEphemerisFunction(
        const std::pair< std::string, std::string > linkEndId, const simulation_setup::NamedBodyMap& bodyMap )
{
    if( bodyMap.count( linkEndId.first ) == 0  )
    {
        std::string errorMessage = "Error when making ephemeris function for " + linkEndId.first + ", " +
                linkEndId.second + ", body not found.";
        throw std::runtime_error( errorMessage );
    }
    return getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >( bodyMap.at( linkEndId.first ), linkEndId );
}

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
template< typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const boost::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& transmitterCompleteEphemeris,
        const boost::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& receiverCompleteEphemeris,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd )
{
    std::vector< boost::shared_ptr< LightTimeCorrection > > lightTimeCorrectionFunctions;

    // Create lighttime correction functions from lightTimeCorrections
    for( unsigned int i = 0; i < lightTimeCorrections.size( ); i++ )
    {

        lightTimeCorrectionFunctions.push_back(
                    createLightTimeCorrections(
                        lightTimeCorrections[ i ], bodyMap, transmittingLinkEnd, receivingLinkEnd ) );
    }

    // Create light time calculator.
    return boost::make_shared< LightTimeCalculator< ObservationScalarType, TimeType > >
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
template< typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
        std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ) )
{

    // Get link end state functions and create light time calculator.
    return createLightTimeCalculator< ObservationScalarType, TimeType >(
                getLinkEndCompleteEphemerisFunction< TimeType, ObservationScalarType >(
                    transmittingLinkEnd, bodyMap ),
                getLinkEndCompleteEphemerisFunction< TimeType, ObservationScalarType >(
                    receivingLinkEnd, bodyMap ),
                bodyMap, lightTimeCorrections, transmittingLinkEnd, receivingLinkEnd );
}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATELIGHTTIMECALCULATOR_H
