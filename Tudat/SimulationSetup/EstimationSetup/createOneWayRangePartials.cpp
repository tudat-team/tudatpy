#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"

namespace tudat
{

namespace observation_partials
{

//! Function to generate one-way range partial wrt an initial position of a body.
boost::shared_ptr< OneWayRangePartial > createOneWayRangePartialWrtBodyPosition(
        const observation_models::LinkEnds oneWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >& lightTimeCorrectionPartialObjects  )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartials =
            createPositionPartialsWrtBodyPosition( oneWayRangeLinkEnds, bodyMap, bodyToEstimate );

    // Create one-range partials if any position partials are created (i.e. if any dependency exists).
    boost::shared_ptr< OneWayRangePartial > oneWayRangePartial;
    if( positionPartials.size( ) > 0 )
    {
        oneWayRangePartial = boost::make_shared< OneWayRangePartial >(
                    oneWayRangeScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    // Return range partial object (NULL if no dependency exists).
    return oneWayRangePartial;
}


}

}


