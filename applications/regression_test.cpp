#include <iostream>

#include <Eigen/Core>

#include "tudat/simulation/simulation.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "liblro.h"


using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::celestial_body_constants;
using namespace tudat::basic_astrodynamics;



Eigen::VectorXd createSimulationInitialState();
std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation
    (const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&,
     const std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>>&);



namespace simulation_constants
{
    const auto simulationDuration = 5 * 113 * 60;  // 565 min, about 5 orbital revolutions
    const auto simulationStart = "2010 SEP 26 06:00:00";
    double simulationStartEpoch;
    double simulationEndEpoch;
    const auto stepSize = 10.0;

    const auto globalFrameOrigin = "Moon";
    const auto globalFrameOrientation = "ECLIPJ2000";
    const auto moonFrame = "MOON_PA";

    const auto lroMass = 1208.0;
    const auto lroArea = 15.38;
    const auto lroCp = 1.41;

    const std::vector<std::string> bodiesToPropagate{"LRO"};
    const std::vector<std::string> centralBodies{"Moon"};
}
using namespace simulation_constants;

void cannonballTargetNew(int i)
{
    const auto resultsFolder = "results/regression/cannonball_new/" + std::to_string(i);

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Moon"}, globalFrameOrigin, globalFrameOrientation);
    bodySettings.at("Sun")->radiationSourceModelSettings =
            // Sun luminosity used by old models is fixed to this value
            isotropicPointRadiationSourceModelSettings(constantLuminosityModelSettings(3.839E26));

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = lroMass;
    bodySettings.at("LRO")->radiationPressureTargetModelSettings =
            cannonballRadiationPressureTargetModelSettings(lroArea, lroCp, {"Moon"});

    // Create bodies
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        pointMassGravityAcceleration(),
                }},
                {"Sun", {
                        radiationPressureAcceleration()
                }},
            }}
    };

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Sun"),
                    relativePositionDependentVariable("Sun", "Moon"),
            };

    auto propagationResults = createAndRunSimulation(
            bodies,
            createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies),
            createSimulationInitialState(),
            dependentVariablesList);
    saveSimulationResults(propagationResults, resultsFolder);
}

void paneledTargetNew(int i)
{
    const auto resultsFolder = "results/regression/paneled_new/" + std::to_string(i);

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Moon"}, globalFrameOrigin, globalFrameOrientation);
    bodySettings.at("Sun")->radiationSourceModelSettings =
            isotropicPointRadiationSourceModelSettings(constantLuminosityModelSettings(3.839E26));

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = lroMass;
    bodySettings.at("LRO")->rotationModelSettings = spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS", "LRO_SC_BUS");
    bodySettings.at("LRO")->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettings({
            TargetPanelSettings(2.82, 0.29, 0.22, false, Eigen::Vector3d::UnitX()),
            TargetPanelSettings(2.82, 0.39, 0.19, false, -Eigen::Vector3d::UnitX()),
            TargetPanelSettings(3.69, 0.32, 0.23, false, Eigen::Vector3d::UnitY()),
            TargetPanelSettings(3.69, 0.32, 0.18, false, -Eigen::Vector3d::UnitY()),
            TargetPanelSettings(5.14, 0.32, 0.18, false, Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(5.14, 0.54, 0.15, false, -Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(11.0, 0.05, 0.05, false, Eigen::Vector3d::UnitX()),
            TargetPanelSettings(11.0, 0.05, 0.05, false, -Eigen::Vector3d::UnitX()),
            TargetPanelSettings(1.0, 0.18, 0.28, false, Eigen::Vector3d::UnitY()),
            TargetPanelSettings(1.0, 0.019, 0.0495, false, -Eigen::Vector3d::UnitY()),
    }, {"Moon"});

    // Create bodies
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        pointMassGravityAcceleration(),
                }},
                {"Sun", {
                        radiationPressureAcceleration()
                }},
            }}
    };

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Sun"),
                    relativePositionDependentVariable("Sun", "Moon"),
            };

    auto propagationResults = createAndRunSimulation(
            bodies,
            createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies),
            createSimulationInitialState(),
            dependentVariablesList);
    saveSimulationResults(propagationResults, resultsFolder);
}

void cannonballTargetOld(int i)
{
    const auto resultsFolder = "results/regression/cannonball_old/" + std::to_string(i);

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Moon"}, globalFrameOrigin, globalFrameOrientation);

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = lroMass;
    bodySettings.at("LRO")->radiationPressureSettings["Sun"] =
            cannonBallRadiationPressureSettings("Sun", lroArea, lroCp, {"Moon"});

    // Create bodies
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        pointMassGravityAcceleration(),
                }},
                {"Sun", {
                        cannonBallRadiationPressureAcceleration()
                }},
            }}
    };

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    singleAccelerationDependentVariable(cannon_ball_radiation_pressure, "LRO", "Sun"),
                    relativePositionDependentVariable("Sun", "Moon"),
            };

    auto propagationResults = createAndRunSimulation(
            bodies,
            createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies),
            createSimulationInitialState(),
            dependentVariablesList);
    saveSimulationResults(propagationResults, resultsFolder);
}

void paneledTargetOld(int i)
{
    const auto resultsFolder = "results/regression/paneled_old/" + std::to_string(i);

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Moon"}, globalFrameOrigin, globalFrameOrientation);

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = lroMass;
    bodySettings.at("LRO")->rotationModelSettings = spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS", "LRO_SC_BUS");
    bodySettings.at("LRO")->radiationPressureSettings["Sun"] =
            panelledRadiationPressureInterfaceSettings("Sun",
               // Emissivities
               {0.29, 0.39, 0.32, 0.32, 0.32, 0.54, 0.05, 0.05, 0.18, 0.019},
               // Areas
               {2.82, 2.82, 3.69, 3.69, 5.14, 5.14, 11.0, 11.0, 1.0, 1.0},
               // Diffuse reflectivity
               {0.22, 0.19, 0.23, 0.18, 0.18, 0.15, 0.05, 0.05, 0.28, 0.0495},
               // Surface normal
               {
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitY(),
                -Eigen::Vector3d::UnitY(),
                Eigen::Vector3d::UnitZ(),
                -Eigen::Vector3d::UnitZ(),
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitY(),
                -Eigen::Vector3d::UnitY()
                }, {"Moon"});

    // Create bodies
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        pointMassGravityAcceleration(),
                }},
                {"Sun", {
                        std::make_shared<AccelerationSettings>(panelled_radiation_pressure_acceleration)
                }},
            }}
    };

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    singleAccelerationDependentVariable(panelled_radiation_pressure_acceleration, "LRO", "Sun"),
                    relativePositionDependentVariable("Sun", "Moon"),
            };

    auto propagationResults = createAndRunSimulation(
            bodies,
            createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies),
            createSimulationInitialState(),
            dependentVariablesList);
    saveSimulationResults(propagationResults, resultsFolder);
}

Eigen::VectorXd createSimulationInitialState()
{
    auto ephemerisLRO = std::make_shared<ephemerides::SpiceEphemeris>(
            "LRO", globalFrameOrigin, false, false, false, globalFrameOrientation);
    auto initialState = ephemerisLRO->getCartesianState(simulationStartEpoch);
    return initialState;
}

std::shared_ptr<propagators::SingleArcSimulationResults<>> createAndRunSimulation(
        const SystemOfBodies& bodies, const AccelerationMap& accelerations, const Eigen::VectorXd& initialState,
        const std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>>& dependentVariablesList)
{
    auto integratorSettings = rungeKuttaFixedStepSettings(
            stepSize, rungeKuttaFehlberg78);

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies,
            accelerations,
            bodiesToPropagate,
            initialState,
            simulationStartEpoch,
            integratorSettings,
            propagationTimeTerminationSettings(simulationEndEpoch),
            cowell,
            dependentVariablesList);

    // Create and run simulation
    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);

    auto propagationResults = dynamicsSimulator->getPropagationResults();
    return std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(propagationResults);
}


int main()
{
    loadLROSpiceKernels();

    int nRepetitions = 10;

    std::cout << "cannonballTargetNew" << std::endl;
    for (int i = 0; i < nRepetitions; ++i)
    {
        std::cout << "  " << i << std::endl;
        cannonballTargetNew(i);
    }

    std::cout << "paneledTargetNew" << std::endl;
    for (int i = 0; i < nRepetitions; ++i)
    {
        std::cout << "  " << i << std::endl;
        paneledTargetNew(i);
    }

    std::cout << "cannonballTargetOld" << std::endl;
    for (int i = 0; i < nRepetitions; ++i)
    {
        std::cout << "  " << i << std::endl;
        cannonballTargetOld(i);
    }

    std::cout << "paneledTargetOld" << std::endl;
    for (int i = 0; i < nRepetitions; ++i)
    {
        std::cout << "  " << i << std::endl;
        paneledTargetOld(i);
    }
}
