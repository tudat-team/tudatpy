#include "Tudat/SimulationSetup/EstimationSetup/createEulerAngleObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/eulerAngleObservablePartials.h"

namespace tudat
{

namespace observation_partials
{


//! Function to generate one-way range partial wrt an initial position of a body.
std::shared_ptr< ObservationPartial< 3 > >  createEulerAngleObservablePartialWrtCurrentOrientation(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier )
{
    std::cout<<"Creating "<<parameterIdentifier.first<<" "<<parameterIdentifier.second.first<<" "<<
               parameterIdentifier.second.second<<std::endl;
    return std::make_shared< EulerAngleObervationPartialWrtCurrentRotationalState >( parameterIdentifier );

}

}

}

