#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/firstOrderRelativisticLightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to get the function returning the light-time correction partial for given correction partial and parameter.

std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
        const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool >
getLightTimeParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const boost::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial )
{
    // Declare return type, set second part to 0 (no dependency found).
    std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
                const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool > partialFunction;
    partialFunction.second = 0;

    // Check type of light-time correction
    switch( lightTimeCorrectionPartial->getCorrectionType( ) )
    {
    // Correction type of 1st-order relativistic
    case observation_models::first_order_relativistic:
    {
        // Check consistency of input.
        boost::shared_ptr< FirstOrderRelativisticLightTimeCorrectionPartial > currentLightTimeCorrectorPartial =
                boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionPartial >( lightTimeCorrectionPartial );
        if( currentLightTimeCorrectorPartial == NULL )
        {
            std::string errorMessage = "Error when getting light time correction partial function, type " +
                    boost::lexical_cast< std::string >( lightTimeCorrectionPartial->getCorrectionType( ) ) +
                       "is inconsistent.";
            throw std::runtime_error( errorMessage );
        }

        // Set partial of gravitational parameter
        if( parameterId.first == estimatable_parameters::gravitational_parameter )
        {
            // Retrieve function from FirstOrderRelativisticLightTimeCorrectionPartial if correction depends on
            // body associated with parameter.
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
