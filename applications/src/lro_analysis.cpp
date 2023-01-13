#include <iostream>

#include <Eigen/Core>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptors.hpp>

#include <tudat/simulation/simulation.h>
#include "tudat/astro/electromagnetism/radiationSourceModel.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::basic_mathematics;
using namespace tudat::electromagnetism;
using namespace tudat::physical_constants;


void loadLROSpiceKernels();
void analyzeSolarIrradiance();
void analyzeIrradianceAtLRO();
void analyzeVariationOfIrradianceWithSubsolarAngle();
void analyzeVariationOfIrradianceWithNumberOfPanels();


const auto simulationStart = "2010 JUN 26 06:00:00";
const auto globalFrameOrigin = "Moon";
const auto globalFrameOrientation = "ECLIPJ2000";


int main()
{
//    analyzeSolarIrradiance();
    analyzeIrradianceAtLRO();
//    analyzeVariationOfIrradianceWithSubsolarAngle();
//    analyzeVariationOfIrradianceWithNumberOfPanels();
}

// Analyze solar irradiance at Earth and Moon
// Result: at simulation time, Earth is at apoapsis -> minimum solar irradiance
void analyzeSolarIrradiance()
{
    loadLROSpiceKernels();

    auto simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);
    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    bodies.at("Sun")->setStateFromEphemeris(simulationStartEpoch);
    bodies.at("Sun")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
    bodies.at("Moon")->setStateFromEphemeris(simulationStartEpoch);
    bodies.at("Earth")->setStateFromEphemeris(simulationStartEpoch);

    auto sunRadiationSourceModel =
            std::dynamic_pointer_cast<IsotropicPointRadiationSourceModel>(bodies.at("Sun")->getRadiationSourceModel());

    auto sunIrradianceAtEarth = sunRadiationSourceModel->evaluateIrradianceAtPosition(
            bodies.at("Earth")->getPosition() - bodies.at("Sun")->getPosition());
    auto sunIrradianceAtMoon = sunRadiationSourceModel->evaluateIrradianceAtPosition(
            bodies.at("Moon")->getPosition() - bodies.at("Sun")->getPosition());

    std::cout << "Sun at Earth: " << sunIrradianceAtEarth << std::endl;
    std::cout << "Sun at Moon: " << sunIrradianceAtMoon << std::endl;
}

// Analyze lunar irradiance at LRO
// Result: angle-based thermal radiation dominates by far
void analyzeIrradianceAtLRO()
{
    loadLROSpiceKernels();

    auto simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);

    for (auto& radiosityModel : std::vector<std::shared_ptr<PanelRadiosityModelSettings>>{
            albedoPanelRadiosityModelSettings(0.12),
            angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95),
            delayedThermalPanelRadiosityModelSettings(0.95)
    })
    {
        // Create planets
        auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);
        bodySettings.at("Moon")->radiationSourceModelSettings =
                staticallyPaneledRadiationSourceModelSettings("Sun", {radiosityModel}, 25000);

        auto bodies = createSystemOfBodies(bodySettings);
        setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

        bodies.at("Sun")->setStateFromEphemeris(simulationStartEpoch);
        bodies.at("Sun")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
        bodies.at("Moon")->setStateFromEphemeris(simulationStartEpoch);
        bodies.at("Moon")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
        bodies.at("Earth")->setStateFromEphemeris(simulationStartEpoch);

        auto moonRadius = bodies.at("Moon")->getShapeModel()->getAverageRadius();
        auto earthRadius = bodies.at("Earth")->getShapeModel()->getAverageRadius();
        double moonEarthDistance = (bodies.at("Earth")->getPosition() - bodies.at("Moon")->getPosition()).norm();

        auto moonRadiationSourceModel = bodies.at("Moon")->getRadiationSourceModel();
        auto sunRadiationSourceModel =
                std::dynamic_pointer_cast<IsotropicPointRadiationSourceModel>(
                        bodies.at("Sun")->getRadiationSourceModel());

        auto sunIrradianceAtMoon = sunRadiationSourceModel->evaluateIrradianceAtPosition(
                bodies.at("Moon")->getPosition() - bodies.at("Sun")->getPosition());

        double theta = 20 / 180. * mathematical_constants::PI;
        Eigen::Vector3d lroPosition = Eigen::Vector3d(cos(theta), sin(theta), 0) * (50e3 + moonRadius);
        auto moonIrradianceAtLRO = moonRadiationSourceModel->evaluateTotalIrradianceAtPosition(
                lroPosition,
                sunIrradianceAtMoon,
                -Eigen::Vector3d::UnitX());

        std::cout << "Moon at LRO: " << moonIrradianceAtLRO << std::endl;
    }
}

// Analyze how the subsolar angle (angle between Sun-Moon and LRO-Moon vectors) changes the received irradiance of LRO
// Result: Falls of to 50% of maximum at 60°, and constant at minimum beyond 90° -> as expected
void analyzeVariationOfIrradianceWithSubsolarAngle()
{
    loadLROSpiceKernels();

    auto simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);

    // Create planets
    auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);
    bodySettings.at("Moon")->radiationSourceModelSettings =
            staticallyPaneledRadiationSourceModelSettings("Sun", {
                albedoPanelRadiosityModelSettings(0.12),
                angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95),
//                delayedThermalPanelRadiosityModelSettings(0.95)
                }, 25000);

    auto bodies = createSystemOfBodies(bodySettings);
    setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

    bodies.at("Sun")->setStateFromEphemeris(simulationStartEpoch);
    bodies.at("Sun")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
    bodies.at("Moon")->setStateFromEphemeris(simulationStartEpoch);
    bodies.at("Moon")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
    bodies.at("Earth")->setStateFromEphemeris(simulationStartEpoch);

    auto moonRadius = bodies.at("Moon")->getShapeModel()->getAverageRadius();

    auto moonRadiationSourceModel = bodies.at("Moon")->getRadiationSourceModel();
    auto sunRadiationSourceModel =
            std::dynamic_pointer_cast<IsotropicPointRadiationSourceModel>(bodies.at("Sun")->getRadiationSourceModel());

    auto sunIrradianceAtMoon = sunRadiationSourceModel->evaluateIrradianceAtPosition(
            bodies.at("Moon")->getPosition() - bodies.at("Sun")->getPosition());

    // MOve from 0° to 180° in steps of 5°
    for (double theta = 0; theta <= mathematical_constants::PI; theta += mathematical_constants::PI / 36)
    {
        Eigen::Vector3d lroPosition = Eigen::Vector3d(cos(theta), sin(theta), 0) * (50e3 + moonRadius);
        auto moonIrradianceAtLRO = moonRadiationSourceModel->evaluateTotalIrradianceAtPosition(
                lroPosition,
                sunIrradianceAtMoon,
                -Eigen::Vector3d::UnitX());

        std::cout << "Moon at LRO: " << moonIrradianceAtLRO << "  " << theta * 180 / mathematical_constants::PI << std::endl;
    }
}

// Analyze how the number of panels changes the received irradiance of LRO
// Result: no large changes over 25000 panels, but significant deviation even with 5000 panels
void analyzeVariationOfIrradianceWithNumberOfPanels()
{
    loadLROSpiceKernels();

    auto simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(simulationStart);

    for (int nPanels : {50, 100, 500, 5000, 25000, 50000, 100000})
    {
        // Create planets
        auto bodySettings = getDefaultBodySettings({"Sun", "Earth", "Moon"}, globalFrameOrigin, globalFrameOrientation);
        bodySettings.at("Moon")->radiationSourceModelSettings =
                staticallyPaneledRadiationSourceModelSettings("Sun", {
                    albedoPanelRadiosityModelSettings(0.12),
                    angleBasedThermalPanelRadiosityModelSettings(100, 375, 0.95),
    //                delayedThermalPanelRadiosityModelSettings(0.95)
                }, nPanels);

        auto bodies = createSystemOfBodies(bodySettings);
        setGlobalFrameBodyEphemerides(bodies.getMap(), globalFrameOrigin, globalFrameOrientation);

        bodies.at("Sun")->setStateFromEphemeris(simulationStartEpoch);
        bodies.at("Sun")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
        bodies.at("Moon")->setStateFromEphemeris(simulationStartEpoch);
        bodies.at("Moon")->getRadiationSourceModel()->updateMembers(simulationStartEpoch);
        bodies.at("Earth")->setStateFromEphemeris(simulationStartEpoch);

        auto moonRadius = bodies.at("Moon")->getShapeModel()->getAverageRadius();
        auto earthRadius = bodies.at("Earth")->getShapeModel()->getAverageRadius();
        double moonEarthDistance = (bodies.at("Earth")->getPosition() - bodies.at("Moon")->getPosition()).norm();

        auto moonRadiationSourceModel = bodies.at("Moon")->getRadiationSourceModel();
        auto sunRadiationSourceModel =
                std::dynamic_pointer_cast<IsotropicPointRadiationSourceModel>(bodies.at("Sun")->getRadiationSourceModel());

        auto sunIrradianceAtMoon = sunRadiationSourceModel->evaluateIrradianceAtPosition(
                bodies.at("Moon")->getPosition() - bodies.at("Sun")->getPosition());
        auto sunIrradianceAtEarth = sunRadiationSourceModel->evaluateIrradianceAtPosition(
                bodies.at("Earth")->getPosition() - bodies.at("Sun")->getPosition());

        auto moonIrradianceAtLRO = moonRadiationSourceModel->evaluateTotalIrradianceAtPosition(
                Eigen::Vector3d(50e3 + moonRadius, 0, 0),
                sunIrradianceAtMoon,
                -Eigen::Vector3d::UnitX());

//        std::cout << "Sun at Moon: " << sunIrradianceAtMoon << std::endl;
//        std::cout << "Sun at Earth: " << sunIrradianceAtEarth << std::endl;
        std::cout << "Moon at LRO: " << moonIrradianceAtLRO << std::endl;
    }
}

// Loads SPICE kernels for LRO
void loadLROSpiceKernels()
{
    using namespace tudat::spice_interface;

    std::string path = "/home/dominik/dev/tudat-bundle/spice/lro/data";

    // Leap seconds
    loadSpiceKernelInTudat(path + "/lsk/naif0012.tls");
    // Planetary shapes
    loadSpiceKernelInTudat(path + "/pck/pck00010.tpc");
    // Planetary gravitational parameters
    loadSpiceKernelInTudat(path + "/pck/gm_de431.tpc");

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
