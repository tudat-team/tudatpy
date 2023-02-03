#include <iostream>

#include <Eigen/Core>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptors.hpp>

#include <tudat/simulation/simulation.h>
#include "tudat/interface/spice/spiceEphemeris.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::physical_constants;
using namespace tudat::gravitation;
using namespace tudat::numerical_integrators;


void loadLROSpiceKernels();
SystemOfBodies createSimulationBodies();
AccelerationMap createSimulationAccelerations(const SystemOfBodies&);
Eigen::VectorXd createSimulationInitialState();
std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation(const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&);
void saveSimulationResults(const std::shared_ptr<propagators::SingleArcSimulationResults<>>& propagationResults);


namespace simulation_constants
{
    const auto resultsFolder = "results/baseline";
    const auto simulationDuration = 5 * 113 * 60;  // 565 min, about 5 orbital revolutions
//    const auto simulationDuration = 500;
    const auto simulationStart = "2010 JUN 26 06:00:00";
    double simulationStartEpoch;
    double simulationEndEpoch;
    const auto printInterval = simulationDuration / 10;
    const auto stepSize = 10.0;

    const auto globalFrameOrigin = "Moon";
    const auto globalFrameOrientation = "ECLIPJ2000";
    const auto moonFrame = "MOON_PA";

    const std::vector<std::string> bodiesToPropagate{"LRO"};
    const std::vector<std::string> centralBodies{"Moon"};
}
using namespace simulation_constants;


int main()
{
    loadLROSpiceKernels();

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    auto bodies = createSimulationBodies();
    auto accelerations = createSimulationAccelerations(bodies);
    auto initialState = createSimulationInitialState();
    auto dynamicsSimulator = createAndRunSimulation(bodies, accelerations, initialState);
    saveSimulationResults(dynamicsSimulator);

    return EXIT_SUCCESS;
}

// Loads SPICE kernels for LRO
void loadLROSpiceKernels()
{
    using namespace tudat::spice_interface;

    std::string path = "spice/lro/data";

    // Leap seconds
    loadSpiceKernelInTudat(path + "/lsk/naif0012.tls");
    // Planetary orientation shapes
    loadSpiceKernelInTudat(path + "/pck/pck00010.tpc");
    // Planetary gravitational parameters
    loadSpiceKernelInTudat(path + "/pck/gm_de431.tpc");
    // Lunar frame
    loadSpiceKernelInTudat(path + "/fk/moon_080317.tf");
    loadSpiceKernelInTudat(path + "/pck/moon_pa_de421_1900_2050.bpc");

    // LRO spacecraft bus and instrument frames
    loadSpiceKernelInTudat(path + "/fk/lro_frames_2012255_v02.tf");
    // LRO spacecraft clock
    loadSpiceKernelInTudat(path + "/sclk/lro_clkcor_2022075_v00.tsc");

    // LRO ephemeris
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path + "/spk"), {})) {
        if (entry.path().extension() == ".bsp")
        {
            loadSpiceKernelInTudat(entry.path().string());
        }
    }

    // LRO orientation
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path + "/ck"), {})) {
        if (entry.path().extension() == ".bc")
        {
            loadSpiceKernelInTudat(entry.path().string());
        }
    }
}

SystemOfBodies createSimulationBodies()
{
    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);

//    bodySettings.at("Sun")->radiationSourceModelSettings =
//            isotropicPointRadiationSourceModelSettings(
//                    constantLuminosityModelSettings(1e33));

    bodySettings.at("Moon")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, moonFrame);
    std::dynamic_pointer_cast<SphericalHarmonicsGravityFieldSettings>(
            bodySettings.at("Moon")->gravityFieldSettings)->resetAssociatedReferenceFrame(moonFrame);
    bodySettings.at("Moon")->radiationSourceModelSettings =
            staticallyPaneledRadiationSourceModelSettings("Sun", {
                albedoPanelRadiosityModelSettings(albedo_dlam1),
                angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95)
            }, 2000, {"Earth"});

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = 1208.0;
    bodySettings.at("LRO")->rotationModelSettings = spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS");
//    bodySettings.at("LRO")->radiationPressureTargetModelSettings =
//            cannonballRadiationPressureTargetModelSettings(15.38, 1.41);
    bodySettings.at("LRO")->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettingsWithOccultationMap({
            TargetPanelSettings(2.82, 0.29, 0.22, Eigen::Vector3d::UnitX()),
            TargetPanelSettings(2.82, 0.39, 0.19, -Eigen::Vector3d::UnitX()),
            TargetPanelSettings(3.69, 0.32, 0.23, Eigen::Vector3d::UnitY()),
            TargetPanelSettings(3.69, 0.32, 0.18, -Eigen::Vector3d::UnitY()),
            TargetPanelSettings(5.14, 0.32, 0.18, Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(5.14, 0.54, 0.15, -Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(11.0, 0.05, 0.05, "Sun"),
            TargetPanelSettings(11.0, 0.05, 0.05, "Sun", false),  // not officially given
            TargetPanelSettings(1.0, 0.18, 0.28, "Earth"),
            TargetPanelSettings(1.0, 0.019, 0.0495, "Earth", false),
    }, {
        // Moon is never occulted as seen from LRO
        {"Sun", {"Earth", "Moon"}}
    });

    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    return bodies;
}

AccelerationMap createSimulationAccelerations(const SystemOfBodies& bodies)
{
    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        sphericalHarmonicAcceleration(100, 100),
                        radiationPressureAcceleration()
                }},
                {"Earth", {
                        sphericalHarmonicAcceleration(50, 50)
                }},
                {"Sun", {
                        pointMassGravityAcceleration(),
                        radiationPressureAcceleration()
                }},
            }}
    };
    return createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);
}

Eigen::VectorXd createSimulationInitialState()
{
    auto ephemerisLRO = std::make_shared<ephemerides::SpiceEphemeris>(
            "LRO", globalFrameOrigin, false, false, false, globalFrameOrientation);
    auto initialState = ephemerisLRO->getCartesianState(simulationStartEpoch);
    return initialState;
}

std::shared_ptr<propagators::SingleArcSimulationResults<>> createAndRunSimulation(
        const SystemOfBodies& bodies, const AccelerationMap& accelerations, const Eigen::VectorXd& initialState)
{
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    keplerianStateDependentVariable("LRO", "Moon"),
                    relativePositionDependentVariable("Sun", "Moon"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "LRO", "Moon"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "LRO", "Earth"),
                    singleAccelerationDependentVariable(point_mass_gravity, "LRO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Moon"),
                    receivedIrradianceDependentVariable("LRO", "Sun"),
                    receivedIrradianceDependentVariable("LRO", "Moon"),
                    receivedFractionDependentVariable("LRO", "Sun"),
                    visibleSourcePanelDependentVariable("LRO", "Moon"),
                    illuminatedSourcePanelDependentVariable("LRO", "Moon"),
                    visibleAndIlluminatedSourcePanelDependentVariable("LRO", "Moon"),
            };

    auto integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances(
            stepSize,
            rungeKuttaFehlberg78,
            stepSize,
            stepSize,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());

    auto outputProcessingSettings = std::make_shared<SingleArcPropagatorProcessingSettings>(
            false, false, 1, TUDAT_NAN,
            std::make_shared<PropagationPrintSettings>(
                    true, false, printInterval, 0, false, false, false, true));

    auto propagatorSettings = translationalStatePropagatorSettings(
            centralBodies,
            accelerations,
            bodiesToPropagate,
            initialState,
            simulationStartEpoch,
            integratorSettings,
            propagationTimeTerminationSettings(simulationEndEpoch),
            cowell,
            dependentVariablesList,
            outputProcessingSettings);

    // Create and run simulation
    auto dynamicsSimulator = createDynamicsSimulator<double, double>(
            bodies, propagatorSettings);

    auto propagationResults = dynamicsSimulator->getPropagationResults();
    return std::dynamic_pointer_cast<SingleArcSimulationResults<double, double>>(propagationResults);
}

void saveSimulationResults(const std::shared_ptr<SingleArcSimulationResults<>>& propagationResults)
{
    auto stateHistory = propagationResults->getEquationsOfMotionNumericalSolution();
    auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
    auto cpuTimeHistory = propagationResults->getCumulativeComputationTimeHistory();
    auto dependentVariableNames = propagationResults->getDependentVariableId();
    auto stateNames = propagationResults->getProcessedStateIds();

    input_output::writeDataMapToTextFile(stateHistory,
                                         "state_history.csv",
                                         resultsFolder,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         resultsFolder,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeDataMapToTextFile(cpuTimeHistory,
                                         "cpu_time.csv",
                                         resultsFolder,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       resultsFolder,
                                       ";");

    input_output::writeIdMapToTextFile(stateNames,
                                       "state_names.csv",
                                       resultsFolder,
                                       ";");
}
