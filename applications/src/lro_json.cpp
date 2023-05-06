#include "lro_json.h"

#include <iostream>
#include <memory>

#include <Eigen/Core>

#include <tudat/simulation/simulation.h>
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "liblro.h"


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
    saveSimulationResults(dynamicsSimulator, settings.saveDir);

    return EXIT_SUCCESS;
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

        if (settings.albedoDistributionMoon == "Constant")
        {
            panelRadiosityModels.push_back(albedoPanelRadiosityModelSettings(0.12));
        }
        else if (settings.albedoDistributionMoon == "DLAM1")
        {
            panelRadiosityModels.push_back(
                    albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1));
        }
        else
        {
            throw std::runtime_error("Invalid albedo_distribution_moon");
        }

        if (settings.thermalTypeMoon == "Delayed")
        {
            panelRadiosityModels.push_back(delayedThermalPanelRadiosityModelSettings(0.95));
        }
        else if (settings.thermalTypeMoon == "AngleBased")
        {
            panelRadiosityModels.push_back(angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95));
        }
        else
        {
            throw std::runtime_error("Invalid thermal_type_moon");
        }

        std::vector<std::string> occultingBodiesForMoon{};
        if (settings.useOccultation)
        {
            occultingBodiesForMoon = {"Earth"};
        }

        if (settings.panelingMoon == "Static")
        {
            bodySettings.at("Moon")->radiationSourceModelSettings =
                    std::make_shared<StaticallyPaneledRadiationSourceModelSettings>(
                            "Sun", panelRadiosityModels, settings.numberOfPanelsMoon, occultingBodiesForMoon);
        }
        else if (settings.panelingMoon == "Dynamic")
        {
            bodySettings.at("Moon")->radiationSourceModelSettings =
                    std::make_shared<DynamicallyPaneledRadiationSourceModelSettings>(
                            "Sun", panelRadiosityModels, settings.numberOfPanelsPerRingMoon, occultingBodiesForMoon);
        }
        else
        {
            throw std::runtime_error("Invalid paneling_moon");
        }
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
            // Moon is never occulted as seen from LRO
            occultingBodiesForLRO = {{"Sun", {"Earth", "Moon"}}};
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
            false, false, 1, TUDAT_NAN,
            std::make_shared<PropagationPrintSettings>(
                    true, false, settings.printInterval, 0, false, false, false, true));

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
