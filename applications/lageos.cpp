#include <iostream>

#include <Eigen/Core>

#include "tudat/simulation/simulation.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"


void checkTLE()
{
    using namespace tudat::ephemerides;

    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::Vector6d state1 = tleEphemeris->getCartesianState(tle->getEpoch());
    Eigen::Vector6d state2 = tleEphemeris->getCartesianState(tle->getEpoch() + 113. * 60 / 4);
    std::cout << state1.segment(3, 3) << "\n\n " << state2.segment(3, 3) << std::endl;
    std::cout << state1.segment(0, 3).norm() << "\n\n " << state2.segment(0, 3).norm() << std::endl;
}

// Simulation that is replicated on Orekit for validation
void orekitLike()
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
    double bodyMass = 406.9;

    const std::vector<std::string> bodiesToPropagate{"Vehicle"};
    const std::vector<std::string> centralBodies{"Earth"};
    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::VectorXd initialState = tleEphemeris->getCartesianState(tle->getEpoch());

    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    Eigen::Vector6d(initialState), EARTH_GRAVITATIONAL_PARAMETER )
            [ orbital_element_conversions::semiMajorAxisIndex ],
            EARTH_GRAVITATIONAL_PARAMETER, bodyMass );

    const auto startTime = spice_interface::convertDateStringToEphemerisTime("1970-01-18");
    const auto endTime = startTime + 10 * orbitalPeriod;

    // Get settings for celestial bodies
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth"}, globalFrameOrigin, globalFrameOrientation);

    // Get settings for vehicle
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at("Vehicle")->constantMass = bodyMass;
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);

    // Create bodies
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"Vehicle", {
                    {"Earth", {
                            pointMassGravityAcceleration(),
                            radiationPressureAcceleration()
                    }}
            }}
    };
    auto accelerations = createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("Vehicle", "Earth"),
                    relativeVelocityDependentVariable("Vehicle", "Earth"),
                    singleAccelerationDependentVariable(radiation_pressure, "Vehicle", "Earth"),
            };

    auto integratorSettings = rungeKuttaFixedStepSettings(
            orbitalPeriod / 100, rungeKuttaFehlberg78);

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies, accelerations, bodiesToPropagate, initialState, startTime, integratorSettings,
            propagationTimeTerminationSettings(endTime), cowell, dependentVariablesList);

    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);
    auto propagationResults =
            std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(dynamicsSimulator->getPropagationResults());
    auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
    auto dependentVariableNames = propagationResults->getDependentVariableId();

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         "results/vv/knocke_tudat",
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       "results/vv/knocke_tudat",
                                       ";");
}

void dynamicPanelingConvergence(int nRings)
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
    double bodyMass = 406.9;

    const std::vector<std::string> bodiesToPropagate{"Vehicle"};
    const std::vector<std::string> centralBodies{"Earth"};
    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::VectorXd initialState = tleEphemeris->getCartesianState(tle->getEpoch());

    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    Eigen::Vector6d(initialState), EARTH_GRAVITATIONAL_PARAMETER )
            [ orbital_element_conversions::semiMajorAxisIndex ],
            EARTH_GRAVITATIONAL_PARAMETER, bodyMass );

    const auto startTime = tle->getEpoch();
    const auto endTime = startTime + 10 * orbitalPeriod;

    // Get settings for celestial bodies
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth"}, globalFrameOrigin, globalFrameOrientation);

    std::vector<int> numberOfPanelsPerRing;
    for (int i = 1; i <= nRings; ++i)
    {
        numberOfPanelsPerRing.push_back(i * 6);
    }

    bodySettings.at("Earth")->radiationSourceModelSettings =
            extendedRadiationSourceModelSettings("Sun", {
                    albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke),
                    delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke)
            }, numberOfPanelsPerRing);

    // Get settings for vehicle
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at("Vehicle")->constantMass = bodyMass;
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);

    // Create bodies
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"Vehicle", {
                    {"Earth", {
                            pointMassGravityAcceleration(),
                            radiationPressureAcceleration()
                    }}
            }}
    };
    auto accelerations = createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("Vehicle", "Earth"),
                    relativeVelocityDependentVariable("Vehicle", "Earth"),
                    singleAccelerationDependentVariable(radiation_pressure, "Vehicle", "Earth"),
            };

    auto integratorSettings = rungeKuttaFixedStepSettings(
            orbitalPeriod / 100, rungeKuttaFehlberg78);

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies, accelerations, bodiesToPropagate, initialState, startTime, integratorSettings,
            propagationTimeTerminationSettings(endTime), cowell, dependentVariablesList);

    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);
    auto propagationResults =
            std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(dynamicsSimulator->getPropagationResults());
    auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
    auto dependentVariableNames = propagationResults->getDependentVariableId();

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         "results/vv/convergence/" + std::to_string(nRings),
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       "results/vv/convergence/" + std::to_string(nRings),
                                       ";");
}

void orekitLikeFindDate(const std::string& startDate)
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
    double bodyMass = 406.9;

    const std::vector<std::string> bodiesToPropagate{"Vehicle"};
    const std::vector<std::string> centralBodies{"Earth"};
    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::VectorXd initialState = tleEphemeris->getCartesianState(tle->getEpoch());

    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    Eigen::Vector6d(initialState), EARTH_GRAVITATIONAL_PARAMETER )
            [ orbital_element_conversions::semiMajorAxisIndex ],
            EARTH_GRAVITATIONAL_PARAMETER, bodyMass );

    const auto startTime = spice_interface::convertDateStringToEphemerisTime(startDate);
    const auto endTime = startTime + 3 * orbitalPeriod;

    // Get settings for celestial bodies
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth"}, globalFrameOrigin, globalFrameOrientation);

    bodySettings.at("Earth")->radiationSourceModelSettings =
            extendedRadiationSourceModelSettings("Sun", {
                    albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke),
                    delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke)
            }, {6, 12});

    // Get settings for vehicle
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at("Vehicle")->constantMass = bodyMass;
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);

    // Create bodies
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"Vehicle", {
                    {"Earth", {
                            pointMassGravityAcceleration(),
                            radiationPressureAcceleration()
                    }}
            }}
    };
    auto accelerations = createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("Vehicle", "Earth"),
                    relativeVelocityDependentVariable("Vehicle", "Earth"),
                    singleAccelerationDependentVariable(radiation_pressure, "Vehicle", "Earth"),
            };

    auto integratorSettings = rungeKuttaFixedStepSettings(
            orbitalPeriod / 100, rungeKuttaFehlberg78);

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies, accelerations, bodiesToPropagate, initialState, startTime, integratorSettings,
            propagationTimeTerminationSettings(endTime), cowell, dependentVariablesList);

    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);
    auto propagationResults =
            std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(dynamicsSimulator->getPropagationResults());
    auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
    auto dependentVariableNames = propagationResults->getDependentVariableId();

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         "results/vv/finddate/" + startDate + "_aprime",
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       "results/vv/finddate/" + startDate + "_aprime",
                                       ";");
}

// Compare static and dynamic paneling
void lageosSimulation(std::string name)
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
    double bodyMass = 406.9;

    const std::vector<std::string> bodiesToPropagate{"Vehicle"};
    const std::vector<std::string> centralBodies{"Earth"};
    const auto globalFrameOrigin = "Earth";
    const auto globalFrameOrientation = "ECLIPJ2000";

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    const auto tle = std::make_shared<Tle>(
            "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
            "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402");
    const auto tleEphemeris = std::make_shared<TleEphemeris>(globalFrameOrigin, globalFrameOrientation, tle);

    Eigen::VectorXd initialState = tleEphemeris->getCartesianState(tle->getEpoch());

    if (name == "sh_static_20000_low" || name == "sh_dynamic_hires_low" || name == "sh_dynamic_double_low")
    {
        Eigen::Vector3d initialPos = initialState.segment(0, 3);
        Eigen::Vector3d initialVel = initialState.segment(3, 3);

        initialPos /= initialPos.norm();
        initialPos *= (1736 + 50.) / 1736. * 6378e3;
        initialVel /= initialVel.norm();
        initialVel *= sqrt(EARTH_GRAVITATIONAL_PARAMETER / initialPos.norm());

        initialState.segment(0, 3) = initialPos;
        initialState.segment(3, 3) = initialVel;
    }

    double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    Eigen::Vector6d(initialState), EARTH_GRAVITATIONAL_PARAMETER )
            [ orbital_element_conversions::semiMajorAxisIndex ],
            EARTH_GRAVITATIONAL_PARAMETER, bodyMass );

    const auto startTime = tle->getEpoch();
    const auto endTime = startTime + 10 * orbitalPeriod;

    // Get settings for celestial bodies
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth"}, globalFrameOrigin, globalFrameOrientation);

    if (name == "const_dynamic_double")
    {
        bodySettings.at("Earth")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings("Sun", {
                        constantPanelRadiosityModelSettings(1000)
                }, {6, 12});
    }
    else if (name == "const_dynamic_hires")
    {
        bodySettings.at("Earth")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings("Sun", {
                        constantPanelRadiosityModelSettings(1000)
                }, {12, 24, 36});
    }
    else if (name == "const_dynamic_superres")
    {
        bodySettings.at("Earth")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings("Sun", {
                        constantPanelRadiosityModelSettings(1000)
                }, {24, 36, 48, 60, 72, 84});
    }
    else if (name == "sh_dynamic_double" || name == "sh_dynamic_double_low")
    {
         bodySettings.at("Earth")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings("Sun", {
                    albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke),
                    delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke)
                }, {6, 12});
    }
    else if (name == "sh_dynamic_hires" || name == "sh_dynamic_hires_low")
    {
        bodySettings.at("Earth")->radiationSourceModelSettings =
                extendedRadiationSourceModelSettings("Sun", {
                    albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke),
                    delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke)
                }, {12, 24, 36});
    }
    else
    {
        throw std::invalid_argument("Invalid name");
    }

    // Get settings for vehicle
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at("Vehicle")->constantMass = bodyMass;
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);

    // Create bodies
    auto bodies = createSystemOfBodies( bodySettings );
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"Vehicle", {
                    {"Earth", {
                            pointMassGravityAcceleration(),
                            radiationPressureAcceleration()
                    }}
            }}
    };
    auto accelerations = createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("Vehicle", "Earth"),
                    relativeVelocityDependentVariable("Vehicle", "Earth"),
                    singleAccelerationDependentVariable(radiation_pressure, "Vehicle", "Earth"),
                    visibleSourceAreaDependentVariable("Vehicle", "Earth")
            };

    auto integratorSettings = rungeKuttaFixedStepSettings(
            orbitalPeriod / 100, rungeKuttaFehlberg78);

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies, accelerations, bodiesToPropagate, initialState, startTime, integratorSettings,
            propagationTimeTerminationSettings(endTime), cowell, dependentVariablesList);

    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);
    auto propagationResults =
            std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(dynamicsSimulator->getPropagationResults());
    auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
    auto dependentVariableNames = propagationResults->getDependentVariableId();

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         "results/static_vs_dynamic/" + name,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       "results/static_vs_dynamic/" + name,
                                       ";");

//    for (const auto& row : dependentVariableHistory)
//    {
//        Eigen::Vector3d position = row.second.segment(0, 3);
//        Eigen::Vector3d velocity = row.second.segment(3, 3);
//        Eigen::Vector3d acceleration = row.second.segment(6, 3);
//
//        Eigen::Vector3d radialUnit = position.normalized();
//        Eigen::Vector3d alongTrackUnit = (velocity - radialUnit * velocity.dot(radialUnit)).normalized();
//        Eigen::Vector3d crossTrackUnit = radialUnit.cross(alongTrackUnit);
//
//        auto radialAcceleration = acceleration.dot(radialUnit);
//        auto alongTrackAcceleration = acceleration.dot(alongTrackUnit);
//        auto crossTrackAcceleration = acceleration.dot(crossTrackUnit);
//    }
}

int main()
{
//    checkTLE();

//    for (std::string name : {
////        "const_dynamic_double", "const_dynamic_hires",
//        "const_dynamic_superres",
////        "sh_dynamic_hires",
//        "sh_dynamic_hires",
////        "sh_dynamic_hires_low",
////        "sh_dynamic_double", "sh_dynamic_double_low"
//    })
//    {
//        std::cout << name << std::endl;
//        lageosSimulation(name);
//    }

//    for (std::string startDate : {
//        "1977-04-22",
//        "1977-02-22",
//        "1977-03-03",
//        "1977-08-25",
//        "1977-08-31"
//    })
//    {
//        std::cout << startDate << std::endl;
//        orekitLikeFindDate(startDate);
//    }

//    for (int nRings = 1; nRings <= 10; ++nRings)
//    {
//        std::cout << nRings << std::endl;
//        dynamicPanelingConvergence(nRings);
//    }

    orekitLike();
}
