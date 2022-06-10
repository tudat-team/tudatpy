/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/estimation_setup/createRelativeAngularPositionPartials.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"

namespace tudat
{


namespace observation_partials
{

//! Function to generate relative angular position partial wrt a position of a body.
std::shared_ptr< RelativeAngularPositionPartial > createRelativeAngularPositionPartialWrtBodyPosition(
        const observation_models::LinkEnds relativeAngularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( relativeAngularPositionLinkEnds, bodies, bodyToEstimate );

    // Create relative angular position partials if any position partials are created (i.e. if any dependency exists).
    std::shared_ptr< RelativeAngularPositionPartial > relativeAngularPositionPartial;
    if( positionPartials.size( ) > 0 )
    {
        relativeAngularPositionPartial = std::make_shared< RelativeAngularPositionPartial >(
                    relativeAngularPositionScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    return relativeAngularPositionPartial;
}

}

}

