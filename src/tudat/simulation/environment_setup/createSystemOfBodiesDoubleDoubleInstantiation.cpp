#include "tudat/simulation/environment_setup/createBodiesFactory.h"

namespace tudat
{

namespace simulation_setup
{

template SystemOfBodies createSystemOfBodies< double, double >( const BodyListSettings& bodySettings );

}  // namespace simulation_setup

}  // namespace tudat
