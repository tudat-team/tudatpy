#include "tudat/simulation/environment_setup/createBodiesFactory.h"

namespace tudat
{

namespace simulation_setup
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
template SystemOfBodies createSystemOfBodies< double, double >( const BodyListSettings& bodySettings );
#endif

}  // namespace simulation_setup

}  // namespace tudat
