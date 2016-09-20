#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
        const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool > getLightTimeParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const boost::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial )
{
    std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
                const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool > partialFunction;
    partialFunction.second = 0;

    switch( lightTimeCorrectionPartial->getCorrectionType( ) )
    {
    case observation_models::first_order_relativistic:
    {
        boost::shared_ptr< FirstOrderRelativisticLightTimeCorrectionPartial > currentLightTimeCorrectorPartial =
                boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionPartial >( lightTimeCorrectionPartial );
        if( currentLightTimeCorrectorPartial == NULL )
        {
            std::cerr<<"Error when getting light time correction partial function, type "<<lightTimeCorrectionPartial->getCorrectionType( )<<
                       "is inconsistent"<<std::endl;
        }
        if( parameterId.first == estimatable_parameters::gravitational_parameter )
        {
            std::vector< std::string > perturbingBodies = currentLightTimeCorrectorPartial->getPerturbingBodies( );
            std::vector< std::string >::iterator findIterator = std::find(
                        perturbingBodies.begin( ), perturbingBodies.end( ), parameterId.second.first );
            if( findIterator != perturbingBodies.end( ) )
            {
                int bodyIndex = std::distance( perturbingBodies.begin( ),  findIterator );
                partialFunction = std::make_pair(
                            boost::bind( &FirstOrderRelativisticLightTimeCorrectionPartial::wrtBodyGravitationalParameter,
                                         currentLightTimeCorrectorPartial, _1, _2, bodyIndex ), 1 );
            }
        }
        break;
    }
    }

    return partialFunction;
}

}

}
