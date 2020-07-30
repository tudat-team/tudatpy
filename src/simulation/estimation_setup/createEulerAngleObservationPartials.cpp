#include "tudat/simulation/estimation_setup/createEulerAngleObservationPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/eulerAngleObservablePartials.h"

namespace tudat
{

namespace observation_partials
{


//! Function to generate one-way range partial wrt an initial position of a body.
std::shared_ptr< ObservationPartial< 3 > >  createEulerAngleObservablePartialWrtCurrentRotationalState(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier )
{
    return std::make_shared< EulerAngleObervationPartialWrtCurrentRotationalState >( parameterIdentifier );

}

}

}

