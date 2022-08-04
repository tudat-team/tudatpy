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
SingleArcDynamicsSimulator<> createSimulator(const SystemOfBodies&, const AccelerationMap&, const Eigen::VectorXd&);
void runSimulationAndSaveResults(SingleArcDynamicsSimulator<>& dynamicsSimulator);


const auto outputFolder = "/home/dominik/dev/tudat-bundle/output/lro/baseline";
const auto simulationDuration = 226 * 60;  // 226 min, about 2 orbital revolutions
const auto simulationStartEpoch =
        (convertCalendarDateToJulianDay(2010, 1, 1, 0, 0, 0) - JULIAN_DAY_ON_J2000) * JULIAN_DAY;
const auto simulationEndEpoch = simulationStartEpoch + simulationDuration;
const auto printInterval = simulationDuration / 10;
const auto minStepSize = 10.0;

const auto globalFrameOrigin = "Moon";
const auto globalFrameOrientation = "ECLIPJ2000";

const std::vector< std::string > bodiesToPropagate{"LRO"};
const std::vector< std::string > centralBodies{"Moon"};


int main()
{
    loadLROSpiceKernels();

    auto bodies = createSimulationBodies();
    auto accelerations = createSimulationAccelerations(bodies);
    auto initialState = createSimulationInitialState();
    auto dynamicsSimulator = createSimulator(bodies, accelerations, initialState);
    runSimulationAndSaveResults(dynamicsSimulator);

    return EXIT_SUCCESS;
}

// Loads reconstructed LRO and DE421 ephemerides
void loadLROSpiceKernels()
{
    std::string path = "/home/dominik/dev/tudat-bundle/spice/lro/data/spk";

    std::vector<std::string> ephemerisKernels;

    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
        if (entry.path().extension() == ".bsp")
        {
            ephemerisKernels.push_back(entry.path().string());
        }
    }

    spice_interface::loadStandardSpiceKernels(ephemerisKernels);
}

SystemOfBodies createSimulationBodies()
{
    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"});
    for(auto const& item : bodySettings.getMap()){
        auto singleBodySettings = item.second;
        singleBodySettings->ephemerisSettings->resetFrameOrientation(globalFrameOrientation);
        singleBodySettings->rotationModelSettings->resetOriginalFrame(globalFrameOrientation);
    }
    auto bodies = createSystemOfBodies(bodySettings);

    // Create LRO
    bodies.createEmptyBody("LRO");
    bodies.getBody("LRO")->setConstantBodyMass(1916.0);
    bodies.getBody("LRO")->setRadiationPressureTargetModel(
            createRadiationPressureTargetModel(
                    std::make_shared<CannonballRadiationPressureTargetModelSettings>(15.38, 1.41),
            "LRO"));

    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    return bodies;
}

AccelerationMap createSimulationAccelerations(const SystemOfBodies& bodies)
{
    SelectedAccelerationMap accelerationMap {
            {"LRO", {
                {"Moon", {
                        sphericalHarmonicAcceleration(2, 2)
                }},
                {"Earth", {
                        sphericalHarmonicAcceleration(2, 2)
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

SingleArcDynamicsSimulator<> createSimulator(
        const SystemOfBodies& bodies, const AccelerationMap& accelerations, const Eigen::VectorXd& initialState)
{
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList
            {
                    relativePositionDependentVariable("LRO", "Moon"),
                    relativeVelocityDependentVariable("LRO", "Moon"),
                    keplerianStateDependentVariable("LRO", "Moon"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "LRO", "Moon"),
                    singleAccelerationDependentVariable(spherical_harmonic_gravity, "LRO", "Earth"),
                    singleAccelerationDependentVariable(point_mass_gravity, "LRO", "Sun")
            };

    auto propagatorSettings =  translationalStatePropagatorSettings (
            centralBodies, accelerations, bodiesToPropagate, initialState, simulationEndEpoch, cowell,
            createDependentVariableSaveSettings(dependentVariablesList), printInterval);

    auto integratorSettings = adamsBashforthMoultonSettings(
                    simulationStartEpoch, minStepSize, minStepSize, 1.0e5);

    SingleArcDynamicsSimulator< > dynamicsSimulator(
            bodies, integratorSettings, propagatorSettings, true, false, false );
    return dynamicsSimulator;
}

void runSimulationAndSaveResults(SingleArcDynamicsSimulator<>& dynamicsSimulator)
{
    auto integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution();
    auto dependentVariableResult = dynamicsSimulator.getDependentVariableHistory();

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "lro_propagation_history.dat",
                                          outputFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          ",");

    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "lro_dependent_variable_history.dat",
                                          outputFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          ",");
}
