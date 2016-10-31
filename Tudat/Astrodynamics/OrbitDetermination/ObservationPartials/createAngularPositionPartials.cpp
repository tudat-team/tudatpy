#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createAngularPositionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/createPositionPartials.h"

namespace tudat
{


namespace observation_partials
{

boost::shared_ptr< AngularPositionPartial > createAngularPositionPartialWrtInitialPosition(
        const observation_models::LinkEnds angularPositionLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< AngularPositionScaling > angularPositionScaler,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >& lightTimeCorrectionPartialObjects )
{
    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartials =
            createPositionPartialsWrtBodyPosition( angularPositionLinkEnds, bodyMap, bodyToEstimate );
    boost::shared_ptr< AngularPositionPartial > angularPositionPartial;

    if( positionPartials.size( ) > 0 )
    {
        angularPositionPartial = boost::make_shared< AngularPositionPartial >(
                    angularPositionScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    return angularPositionPartial;
}

}

}

