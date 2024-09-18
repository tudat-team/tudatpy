/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include <Eigen/Core>

#include "tudat/simulation/simulation.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/io/applicationOutput.h"

int main()
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::celestial_body_constants;
    using namespace tudat::basic_astrodynamics;

    double area = 0.28483;
    double coefficient = 1.13;
    double lageosMass = 406.9;

    const std::vector<std::string> bodiesToPropagate{"LAGEOS"};
    const std::vector<std::string> centralBodies{"Earth"};
    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Get initial state from two-line elements
    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::VectorXd initialState = tleEphemeris->getCartesianState(tle->getEpoch());

    // Set propagation period to 10 revolutions
    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    Eigen::Vector6d(initialState), EARTH_GRAVITATIONAL_PARAMETER )
            [ orbital_element_conversions::semiMajorAxisIndex ],
            EARTH_GRAVITATIONAL_PARAMETER, lageosMass );

    const auto startTime = tle->getEpoch();
    const auto endTime = startTime + 10 * orbitalPeriod;
    const auto printInterval = 2 * orbitalPeriod;

    // Get settings for celestial bodies
    // Earth's default radiation source model is from Knocke (1988) and includes albedo + thermal
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth"}, globalFrameOrigin, globalFrameOrientation);

    // Get settings for LAGEOS
    bodySettings.addSettings( "LAGEOS" );
    bodySettings.at("LAGEOS")->constantMass = lageosMass;
    bodySettings.at("LAGEOS")->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);

    // Create bodies
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    // Create accelerations
    SelectedAccelerationMap accelerationMap{
            {"LAGEOS", {
                    {"Earth", {
                            pointMassGravityAcceleration(),
                            radiationPressureAcceleration()
                    }}
            }}
    };
    auto accelerations = createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);

    // Set up simulation
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LAGEOS", "Earth"),
                    relativeVelocityDependentVariable("LAGEOS", "Earth"),
                    singleAccelerationDependentVariable(radiation_pressure, "LAGEOS", "Earth"),
            };

    auto integratorSettings = rungeKuttaFixedStepSettings(
            orbitalPeriod / 100, rungeKuttaFehlberg78);

    auto outputProcessingSettings = std::make_shared<SingleArcPropagatorProcessingSettings>(
            false, false, 1, TUDAT_NAN,
            std::make_shared<PropagationPrintSettings>(
                    true, false, printInterval, 0, false, false, false, true));

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies, accelerations, bodiesToPropagate, initialState, startTime, integratorSettings,
            propagationTimeTerminationSettings(endTime), cowell, dependentVariablesList, outputProcessingSettings);

    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);
    auto propagationResults =
            std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(dynamicsSimulator->getPropagationResults());

    // Store results
    std::string outputSubFolder = "output_LageosRadiationPressureAcceleration/";
    input_output::writeDataMapToTextFile(propagationResults->getDependentVariableHistory(),
                                         "dependent_variable_history.csv",
                                         outputSubFolder,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");
    input_output::writeIdMapToTextFile(propagationResults->getDependentVariableId(),
                                       "dependent_variable_names.csv",
                                       outputSubFolder,
                                       ",");
    std::cout << "Results stored to " << "output_LageosRadiationPressureAcceleration/" << std::endl;
}
