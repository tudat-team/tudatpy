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


SystemOfBodies createSimulationBodies();
AccelerationMap createSimulationAccelerations(const SystemOfBodies&);
Eigen::VectorXd createSimulationInitialState();
std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation(const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&);


namespace simulation_constants
{
    const auto simulationDuration = 1 * 113 * 60;  // 565 min, about 5 orbital revolutions
//    const auto simulationDuration = 500;
    const auto simulationStart = "2010 SEP 29 07:47:00";
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

    return EXIT_SUCCESS;
}

SystemOfBodies createSimulationBodies()
{
    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);

//    bodySettings.at("Sun")->radiationSourceModelSettings =
//            isotropicPointRadiationSourceModelSettings(
//                    constantLuminosityModelSettings(1e33));

    bodySettings.at("Moon")->shapeModelSettings = oblateSphericalBodyShapeSettings(1738.1e3, 0.0012);
    bodySettings.at("Moon")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, moonFrame, moonFrame);
//    bodySettings.at("Moon")->radiationSourceModelSettings =
//            staticallyPaneledRadiationSourceModelSettings("Sun", {
//                albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1),
//                angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95)
//            }, 2000, {"Earth"});
    bodySettings.at("Moon")->radiationSourceModelSettings =
            extendedRadiationSourceModelSettings("Sun", {
                albedoPanelRadiosityModelSettings(SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1),
                angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95)
            }, {6, 12, 18, 24, 30});

    // Create LRO
    bodySettings.addSettings("LRO");
    bodySettings.at("LRO")->constantMass = 1208.0;
    bodySettings.at("LRO")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, "LRO_SC_BUS", "LRO_SC_BUS");
//    bodySettings.at("LRO")->radiationPressureTargetModelSettings =
//            cannonballRadiationPressureTargetModelSettings(11.52, 1.04);
    bodySettings.at("LRO")->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettings({
            TargetPanelSettings(2.82, 0.29, 0.22, false, Eigen::Vector3d::UnitX()),
            TargetPanelSettings(2.82, 0.39, 0.19, false, -Eigen::Vector3d::UnitX()),
            TargetPanelSettings(3.69, 0.32, 0.23, false, Eigen::Vector3d::UnitY()),
            TargetPanelSettings(3.69, 0.32, 0.18, false, -Eigen::Vector3d::UnitY()),
            TargetPanelSettings(5.14, 0.32, 0.18, false, Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(5.14, 0.54, 0.15, false, -Eigen::Vector3d::UnitZ()),
            TargetPanelSettings(11.0, 0.05, 0.05, false, "Sun"),
            TargetPanelSettings(11.0, 0.05, 0.05, false, "Sun", false),  // not officially given
            TargetPanelSettings(1.0, 0.18, 0.28, false, "Earth"),
            TargetPanelSettings(1.0, 0.019, 0.0495, false, "Earth", false),
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
                        pointMassGravityAcceleration(),
                        radiationPressureAcceleration()
                }},
                {"Sun", {
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
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList {};

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
