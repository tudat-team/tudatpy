#include "lro_json.h"

#include <iostream>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <Eigen/Core>

#include <tudat/simulation/simulation.h>
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceInterface.h"


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

using namespace simulation_constants;



static Settings settings;


int main(int argc, char* argv[])
{
    loadLROSpiceKernels();
    settings = loadSettings(argv[1]);

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

    bodySettings.at("Moon")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, moonFrame);
    std::dynamic_pointer_cast<SphericalHarmonicsGravityFieldSettings>(
            bodySettings.at("Moon")->gravityFieldSettings)->resetAssociatedReferenceFrame(moonFrame);

    if (settings.useMoonRadiation) {
        std::vector<std::shared_ptr<PanelRadiosityModelSettings>> panelRadiosityModels;
        panelRadiosityModels.push_back(albedoPanelRadiosityModelSettings(0.12, settings.useInstantaneousReradiation));

        if (settings.thermalType == "Delayed")
        {
            panelRadiosityModels.push_back(delayedThermalPanelRadiosityModelSettings(0.95));
        }
        else if (settings.thermalType == "AngleBased")
        {
            panelRadiosityModels.push_back(angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95));
        }
        else
        {
            throw std::runtime_error("Invalid thermal_type");
        }

        std::vector<std::string> occultingBodiesForMoon{};
        if (settings.useOccultation)
        {
            occultingBodiesForMoon = {"Earth"};
        }

        bodySettings.at("Moon")->radiationSourceModelSettings =
                std::make_shared<StaticallyPaneledRadiationSourceModelSettings>(
                        "Sun", panelRadiosityModels, settings.numberOfPanelsMoon, occultingBodiesForMoon);
    }

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = 1208.0;
    bodySettings.at("LRO")->rotationModelSettings = spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS");
    if (settings.targetType == "Cannonball") {
        bodySettings.at("LRO")->radiationPressureTargetModelSettings =
                cannonballRadiationPressureTargetModelSettings(15.38, 1.41);
    }
    else if (settings.targetType == "Paneled")
    {
        std::map<std::string, std::vector<std::string>> occultingBodiesForLRO{};
        if (settings.useOccultation)
        {
            // TODO once multiple occultation is supported, add occultation by Earth as well
            occultingBodiesForLRO = {{"Sun", {"Moon"}}};
        }

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
        }, occultingBodiesForLRO);
    }
    else
    {
        throw std::runtime_error("Invalid target_type");
    }

    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    return bodies;
}

AccelerationMap createSimulationAccelerations(const SystemOfBodies& bodies)
{
    SelectedAccelerationMap accelerationMap{
            {"LRO", {
                {"Moon", {
                        sphericalHarmonicAcceleration(100, 100)
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

    if (settings.useMoonRadiation)
    {
        accelerationMap["LRO"]["Moon"].push_back(radiationPressureAcceleration());
    }

    return createAccelerationModelsMap(bodies, accelerationMap, bodiesToPropagate, centralBodies);
}

Eigen::VectorXd createSimulationInitialState()
{
    auto ephemerisLRO = std::make_shared<ephemerides::SpiceEphemeris>(
            "LRO", globalFrameOrigin, false, false, false, globalFrameOrientation);
    auto initialState = ephemerisLRO->getCartesianState(settings.simulationStartEpoch);
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
                    receivedIrradianceDependentVariable("LRO", "Sun"),
                    receivedFractionDependentVariable("LRO", "Sun"),
            };
    if (settings.useMoonRadiation)
    {
        dependentVariablesList.insert(dependentVariablesList.end(), {
                singleAccelerationDependentVariable(radiation_pressure, "LRO", "Moon"),
                receivedIrradianceDependentVariable("LRO", "Moon"),
                visibleSourcePanelDependentVariable("LRO", "Moon"),
                illuminatedSourcePanelDependentVariable("LRO", "Moon"),
                visibleAndIlluminatedSourcePanelDependentVariable("LRO", "Moon"),
        });
    }

    auto integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances(
            settings.stepSize,
            rungeKuttaFehlberg78,
            settings.stepSize,
            settings.stepSize,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());

    auto outputProcessingSettings = std::make_shared<SingleArcPropagatorProcessingSettings>(
            false,
            false,
            std::make_shared<PropagationPrintSettings>(
                    true, false, false, settings.printInterval));

    auto propagatorSettings =  translationalStatePropagatorSettings(
            centralBodies,
            accelerations,
            bodiesToPropagate,
            initialState,
            settings.simulationStartEpoch,
            integratorSettings,
            propagationTimeTerminationSettings(settings.simulationEndEpoch),
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
    auto stateNames = propagationResults->getStateIds();

    input_output::writeDataMapToTextFile(stateHistory,
                                         "state_history.csv",
                                         settings.saveDir,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeDataMapToTextFile(dependentVariableHistory,
                                         "dependent_variable_history.csv",
                                         settings.saveDir,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeDataMapToTextFile(cpuTimeHistory,
                                         "cpu_time.csv",
                                         settings.saveDir,
                                         "",
                                         std::numeric_limits< double >::digits10,
                                         std::numeric_limits< double >::digits10,
                                         ",");

    input_output::writeIdMapToTextFile(dependentVariableNames,
                                       "dependent_variable_names.csv",
                                       settings.saveDir,
                                       ",");

    input_output::writeIdMapToTextFile(stateNames,
                                       "state_names.csv",
                                       settings.saveDir,
                                       ",");
}
