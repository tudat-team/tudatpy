#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace simulation_setup
{

#if TUDAT_BUILD_ALL_TESTS
template std::shared_ptr< propagators::VariationalEquationsSolver< double, double > > createVariationalEquationsSolver< double, double >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
        const bool integrateEquationsOnCreation );
#endif

}  // namespace simulation_setup

}  // namespace tudat
