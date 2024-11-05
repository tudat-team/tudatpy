#include <iostream>
#include <memory>
#include <string>

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


SystemOfBodies createSimulationBodies();
AccelerationMap createSimulationAccelerations(const SystemOfBodies&);
Eigen::VectorXd createSimulationInitialState();
std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation(const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&);


namespace simulation_constants
{
    const auto resultsFolder = "results/baseline";
    const auto simulationDuration = 2 * 113 * 60;  // 565 min, about 5 orbital revolutions
//    const auto simulationDuration = 500;
    const std::string simulationStart = "2010 JUN 29 07:47:00";
    double simulationStartEpoch;
    double simulationEndEpoch;
    const auto printInterval = simulationDuration / 10;
    const auto stepSize = 10.0;

    const auto globalFrameOrigin = "Moon";
    const auto globalFrameOrientation = "ECLIPJ2000";
    const auto moonFrame = "MOON_ME";

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
    saveSimulationResults(dynamicsSimulator, resultsFolder, true);

    return EXIT_SUCCESS;
}

SystemOfBodies createSimulationBodies()
{
    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);

//    bodySettings.at("Sun")->radiationSourceModelSettings =
//            isotropicPointRadiationSourceModelSettings(
//                    constantLuminosityModelSettings(1e33));

    bodySettings.at("Moon")->shapeModelSettings = sphericalBodyShapeSettings(1737.4e3);
    bodySettings.at("Moon")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, moonFrame, moonFrame);
    std::dynamic_pointer_cast<SphericalHarmonicsGravityFieldSettings>(
            bodySettings.at("Moon")->gravityFieldSettings)->resetAssociatedReferenceFrame(moonFrame);
//    bodySettings.at("Moon")->radiationSourceModelSettings =
//            staticallyPaneledRadiationSourceModelSettings("Sun", {
//                albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1),
//                angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95)
//            }, 2000, {"Earth"});
    bodySettings.at("Moon")->radiationSourceModelSettings =
            extendedRadiationSourceModelSettings({
                albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun"),
                angleBasedThermalPanelRadiosityModelSettings(95, 385, 0.95, "Sun")
            }, {6, 12, 18, 24, 30}, {"Earth"});

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = 1208.0;
    bodySettings.at("LRO")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS", "LRO_SC_BUS");
//    bodySettings.at("LRO")->radiationPressureTargetModelSettings =
//            cannonballRadiationPressureTargetModelSettingsWithOccultationMap(11.52, 1.04, {{"Sun", {"Earth", "Moon"}}});

    // Sun is only tracked for beta < 30 deg (Mazarico 2018)
    bool trackSun = false;
    for (const auto& month: {"MAR", "APR", "SEP", "OCT"})
    {
        if (simulationStart.find(month) != std::string::npos) {
            trackSun = true;
        }
    }

    if (trackSun)
    {
        bodySettings.at("LRO")->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettingsWithOccultationMap({
                TargetPanelSettings(2.82, 0.29, 0.22, true, Eigen::Vector3d::UnitX()),
                TargetPanelSettings(2.82, 0.39, 0.19, true, -Eigen::Vector3d::UnitX()),
                TargetPanelSettings(3.69, 0.32, 0.23, true, Eigen::Vector3d::UnitY()),
                TargetPanelSettings(3.69, 0.32, 0.18, true, -Eigen::Vector3d::UnitY()),
                TargetPanelSettings(5.14, 0.32, 0.18, true, Eigen::Vector3d::UnitZ()),
                TargetPanelSettings(5.14, 0.54, 0.15, true, -Eigen::Vector3d::UnitZ()),
                TargetPanelSettings(11.0, 0.05, 0.05, true, "Sun"),
                TargetPanelSettings(11.0, 0.30, 0.20, true, "Sun", false),  // not officially given
                TargetPanelSettings(1.0, 0.18, 0.28, true, "Earth"),
                TargetPanelSettings(1.0, 0.019, 0.0495, true, "Earth", false),
        }, {
            // Moon is never occulted as seen from LRO
            {"Sun", {"Earth", "Moon"}}
        });
    }
    else
    {
        bodySettings.at("LRO")->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettingsWithOccultationMap({
                TargetPanelSettings(2.82, 0.29, 0.22, true, Eigen::Vector3d::UnitX()),
                TargetPanelSettings(2.82, 0.39, 0.19, true, -Eigen::Vector3d::UnitX()),
                TargetPanelSettings(3.69, 0.32, 0.23, true, Eigen::Vector3d::UnitY()),
                TargetPanelSettings(3.69, 0.32, 0.18, true, -Eigen::Vector3d::UnitY()),
                TargetPanelSettings(5.14, 0.32, 0.18, true, Eigen::Vector3d::UnitZ()),
                TargetPanelSettings(5.14, 0.54, 0.15, true, -Eigen::Vector3d::UnitZ()),
                TargetPanelSettings(11.0, 0.05, 0.05, true, Eigen::Vector3d(-1, -1, 0).normalized()),
                TargetPanelSettings(11.0, 0.30, 0.20, true, -Eigen::Vector3d(-1, -1, 0).normalized()),
                TargetPanelSettings(1.0, 0.18, 0.28, true, "Earth"),
                TargetPanelSettings(1.0, 0.019, 0.0495, true, "Earth", false),
        }, {
            // Moon is never occulted as seen from LRO
            {"Sun", {"Earth", "Moon"}}
        });
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
                        sphericalHarmonicAcceleration(100, 100),
                        radiationPressureAcceleration()
                }},
                {"Earth", {
                        pointMassGravityAcceleration(),
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
                    altitudeDependentVariable("LRO", "Moon"),
                    relativePositionDependentVariable("Sun", "Moon"),
                    relativePositionDependentVariable("Earth", "Moon"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "LRO", "Moon"),
                    singleAccelerationDependentVariable(point_mass_gravity, "LRO", "Earth"),
                    singleAccelerationDependentVariable(point_mass_gravity, "LRO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "LRO", "Moon"),
                    latitudeDependentVariable("LRO", "Moon"),
                    longitudeDependentVariable("LRO", "Moon"),
                    receivedIrradianceDependentVariable("LRO", "Sun"),
                    receivedIrradianceDependentVariable("LRO", "Moon"),
                    receivedFractionDependentVariable("LRO", "Sun"),
                    visibleAndEmittingSourcePanelCountDependentVariable("LRO", "Moon"),
                    visibleSourceAreaDependentVariable("LRO", "Moon")
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
