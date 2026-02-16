/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef ORBITDETERMINATIONTESTCASEUTILITIES_H
#define ORBITDETERMINATIONTESTCASEUTILITIES_H

#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/compareEstimationAndCovarianceResultsTestCase.h"

namespace tudat
{
namespace unit_tests
{

// Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::physical_constants;

Eigen::VectorXd getDefaultInitialParameterPerturbation( );

}  // namespace unit_tests

}  // namespace tudat

#endif  // ORBITDETERMINATIONTESTCASEUTILITIES_H
