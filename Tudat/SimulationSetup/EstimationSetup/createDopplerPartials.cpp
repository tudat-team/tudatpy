
/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EstimationSetup/createDopplerPartials.h"

namespace tudat
{

namespace observation_partials
{

using namespace ephemerides;

//! Function to generate one-way doppler partial wrt an initial position of a body.
boost::shared_ptr< OneWayDopplerPartial > createOneWayDopplerPartialWrtBodyState(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< PositionPartialScaling > oneWayDopplerScaler,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects  )
{
    if( boost::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ) == NULL )
    {
        throw std::runtime_error( "Error, expected one-way doppler scaling when making one-way doppler partial" );
    }

    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( oneWayDopplerLinkEnds, bodyMap, bodyToEstimate );

    // Create one-doppler partials if any position partials are created (i.e. if any dependency exists).
    boost::shared_ptr< OneWayDopplerPartial > oneWayDopplerPartial;
    if( positionPartials.size( ) > 0 )
    {
        oneWayDopplerPartial = boost::make_shared< OneWayDopplerPartial >(
                    boost::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ),
                    positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    // Return doppler partial object (NULL if no dependency exists).
    return oneWayDopplerPartial;
}

//! Function to create an object that computes the scaling of the state partials to obtain proper time rate partials
boost::shared_ptr< OneWayDopplerProperTimeComponentScaling > createDopplerProperTimePartials(
        const boost::shared_ptr< observation_models::DopplerProperTimeRateInterface > dopplerProperTimeInterface,
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const observation_models::LinkEndType linkEndAtWhichPartialIsComputed )
{
    boost::shared_ptr< OneWayDopplerProperTimeComponentScaling >  properTimeRateDopplerPartial = NULL;
    if( dopplerProperTimeInterface == NULL )
    {
        properTimeRateDopplerPartial = NULL;
    }
    else if( boost::dynamic_pointer_cast< observation_models::DirectFirstOrderDopplerProperTimeRateInterface >(
                 dopplerProperTimeInterface ) != NULL )
    {
        bool computeStatePartials = ( oneWayDopplerLinkEnds.at( linkEndAtWhichPartialIsComputed ).first !=
                boost::dynamic_pointer_cast< observation_models::DirectFirstOrderDopplerProperTimeRateInterface >(
                    dopplerProperTimeInterface )->getCentralBody( ) );
        properTimeRateDopplerPartial = boost::make_shared< OneWayDopplerDirectFirstOrderProperTimeComponentScaling >(
                    boost::dynamic_pointer_cast< observation_models::DirectFirstOrderDopplerProperTimeRateInterface >(
                        dopplerProperTimeInterface ), linkEndAtWhichPartialIsComputed, computeStatePartials );
    }
    else
    {
        std::cerr << "Warning, proper time contribution to Doppler observable not incorporated into Doppler partial " << std::endl;
        properTimeRateDopplerPartial = NULL;
    }
    return properTimeRateDopplerPartial;

}

}

}
