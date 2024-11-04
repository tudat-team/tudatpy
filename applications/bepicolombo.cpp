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
using namespace tudat::unit_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::physical_constants;
using namespace tudat::gravitation;
using namespace tudat::numerical_integrators;


void loadMPOSpiceKernels();
SystemOfBodies createSimulationBodies();
AccelerationMap createSimulationAccelerations(const SystemOfBodies&);
Eigen::VectorXd createSimulationInitialState();
std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation(const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&);


namespace simulation_constants
{
    const auto resultsFolder = "results/mpo_1year";
//    const auto simulationDuration = 2 * 8355.0; // 2 orbits
//    const auto simulationDuration = 13 * 365 * 24 * 3600; // 13 years
    const auto simulationDuration = 1 * 365 * 24 * 3600; // 1 year
    const auto simulationStart = "2026 MAR 15";  // one day after orbit insertion
    double simulationStartEpoch;
    double simulationEndEpoch;
    const auto printInterval = simulationDuration / 10;
    const auto stepSize = 8355.0 / 360;  // ~23 s, 1° of orbit

    const auto globalFrameOrigin = "Mercury";
    const auto globalFrameOrientation = "ECLIPJ2000";
    const auto mercuryFrame = "IAU_Mercury";

    const auto mercuryGravitationalParameter = 3.3e23 * GRAVITATIONAL_CONSTANT;

    const std::vector<std::string> bodiesToPropagate{"MPO"};
    const std::vector<std::string> centralBodies{"Mercury"};
}
using namespace simulation_constants;


int main()
{
    loadMPOSpiceKernels();

    simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);
    simulationEndEpoch = simulationStartEpoch + simulationDuration;

    auto bodies = createSimulationBodies();
    auto accelerations = createSimulationAccelerations(bodies);
    auto initialState = createSimulationInitialState();
    auto dynamicsSimulator = createAndRunSimulation(bodies, accelerations, initialState);
    saveSimulationResults(dynamicsSimulator, resultsFolder);

    return EXIT_SUCCESS;
}

void loadMPOSpiceKernels()
{
    using namespace tudat::spice_interface;

    std::string path = "spice/bepicolombo/kernels";

    // Leap seconds
    loadSpiceKernelInTudat(path + "/lsk/naif0012.tls");
    // Planetary orientation shapes
    loadSpiceKernelInTudat(path + "/pck/pck00010.tpc");
    // Planetary gravitational parameters
    loadSpiceKernelInTudat(path + "/pck/gm_de431.tpc");

    // MPO spacecraft bus and instrument frames
    loadSpiceKernelInTudat(path + "/fk/bc_mpo_v33.tf");
    // MPO science frames
    loadSpiceKernelInTudat(path + "/fk/bc_sci_v12.tf");

    // MPO ephemeris
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path + "/spk"), {})) {
        if (entry.path().extension() == ".bsp")
        {
            loadSpiceKernelInTudat(entry.path().string());
        }
    }
}

SystemOfBodies createSimulationBodies()
{
    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Mercury"}, globalFrameOrigin, globalFrameOrientation);

    bodySettings.at("Mercury")->rotationModelSettings =
            spiceRotationModelSettings(globalFrameOrientation, mercuryFrame, mercuryFrame);
    bodySettings.at("Mercury")->radiationSourceModelSettings =
            extendedRadiationSourceModelSettings("Sun", {
                albedoPanelRadiosityModelSettings(0.12),
            }, {6, 12, 18, 24, 30});


    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -2.7e-5;
    cosineCoefficients( 2, 2 ) = 1.6e-5;

    bodySettings.at("Mercury")->gravityFieldSettings = sphericalHarmonicsGravitySettings(
                mercuryGravitationalParameter, 2439e3, cosineCoefficients, sineCoefficients, mercuryFrame );

    // Create MPO
    bodySettings.addSettings("MPO");
    bodySettings.at("MPO")->constantMass = 357;
//    bodySettings.at("MPO")->radiationPressureTargetModelSettings =
//            cannonballRadiationPressureTargetModelSettings(1.9e-2 * 357, 1.0);
    bodySettings.at("MPO")->radiationPressureTargetModelSettings =
            cannonballRadiationPressureTargetModelSettingsWithOccultationMap(1.9e-2 * 357, 1.0, {{"Sun", {"Mercury"}}});

    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    return bodies;
}

AccelerationMap createSimulationAccelerations(const SystemOfBodies& bodies)
{
    SelectedAccelerationMap accelerationMap{
            {"MPO", {
                {"Mercury", {
                        sphericalHarmonicAcceleration(2, 2),
                        radiationPressureAcceleration()
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
//    auto ephemerisMPO = std::make_shared<ephemerides::SpiceEphemeris>(
//            "MPO", globalFrameOrigin, false, false, false, globalFrameOrientation);
//    auto initialStateInGlobalCartesianElements = ephemerisMPO->getCartesianState(simulationStartEpoch);

    Eigen::Vector6d initialStateInMercuryKeplerianElements(6);
    initialStateInMercuryKeplerianElements(semiMajorAxisIndex ) = 3389e3;
    initialStateInMercuryKeplerianElements(eccentricityIndex ) = 0.162;
    initialStateInMercuryKeplerianElements(inclinationIndex ) = convertDegreesToRadians(90.);
    initialStateInMercuryKeplerianElements(argumentOfPeriapsisIndex )
            = convertDegreesToRadians(0.7);
    initialStateInMercuryKeplerianElements(longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians(180);
    initialStateInMercuryKeplerianElements(trueAnomalyIndex ) = convertDegreesToRadians(0);

    const Eigen::Vector6d initialStateInMercuryCartesianElements = convertKeplerianToCartesianElements(
            initialStateInMercuryKeplerianElements,
            mercuryGravitationalParameter );

    const Eigen::Matrix6d rotationFromMercuryToGlobalFrame = spice_interface::computeStateRotationMatrixBetweenFrames(
            mercuryFrame, globalFrameOrientation, simulationStartEpoch);

    const Eigen::Vector6d initialStateInGlobalCartesianElements =
            rotationFromMercuryToGlobalFrame * initialStateInMercuryCartesianElements;

    return initialStateInGlobalCartesianElements;
}

std::shared_ptr<propagators::SingleArcSimulationResults<>> createAndRunSimulation(
        const SystemOfBodies& bodies, const AccelerationMap& accelerations, const Eigen::VectorXd& initialState)
{
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings>> dependentVariablesList
            {
                    relativePositionDependentVariable("MPO", "Mercury"),
                    relativeVelocityDependentVariable("MPO", "Mercury"),
                    keplerianStateDependentVariable("MPO", "Mercury"),
                    altitudeDependentVariable("MPO", "Mercury"),
                    relativePositionDependentVariable("Sun", "Mercury"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "MPO", "Mercury"),
                    singleAccelerationDependentVariable(point_mass_gravity, "MPO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "MPO", "Sun"),
                    singleAccelerationDependentVariable(radiation_pressure, "MPO", "Mercury"),
                    receivedIrradianceDependentVariable("MPO", "Sun"),
                    receivedIrradianceDependentVariable("MPO", "Mercury"),
                    receivedFractionDependentVariable("MPO", "Sun"),
                    visibleAndEmittingSourcePanelCountDependentVariable("MPO", "Mercury"),
                    visibleSourceAreaDependentVariable("MPO", "Mercury")
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
