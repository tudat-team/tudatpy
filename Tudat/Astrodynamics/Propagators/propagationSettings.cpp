#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

//! Get size of state for single propagated state of given type.
int getSingleIntegrationSize( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 6;
        break;
    case body_mass_state:
        singleStateSize = 1;
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + boost::lexical_cast< std::string >( stateType ) + "when getting size";
       throw std::runtime_error( errorMessage );
    }
    return singleStateSize;
}

//! Get order of differential equation for governing equations of dynamics of given type.
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 2;
        break;
    case body_mass_state:
        singleStateSize = 1;
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + boost::lexical_cast< std::string >( stateType ) + "when getting order";
       throw std::runtime_error( errorMessage );
    }
    return singleStateSize;
}


}

}
