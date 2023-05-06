#ifndef TUDATBUNDLE_LRO_JSON_H
#define TUDATBUNDLE_LRO_JSON_H

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <algorithm>
#include <iterator>

#include <nlohmann/json.hpp>

#include <tudat/simulation/simulation.h>

using namespace tudat;



struct Settings;

std::string formatEphemerisTime(double et);

Settings loadSettings(char* path);

std::ostream& operator<<(std::ostream& os, const Settings& settings);


simulation_setup::SystemOfBodies createSimulationBodies();

basic_astrodynamics::AccelerationMap createSimulationAccelerations(const simulation_setup::SystemOfBodies&);

Eigen::VectorXd createSimulationInitialState();

std::shared_ptr< propagators::SingleArcSimulationResults<>> createAndRunSimulation(
        const simulation_setup::SystemOfBodies&, const basic_astrodynamics::AccelerationMap&, const Eigen::VectorXd&);





struct Settings
{
    std::string id{};
    std::string saveDir{};
    std::string simulationStart{};
    double simulationDuration{};
    double simulationStartEpoch{};
    double simulationEndEpoch{};
    double printInterval{};
    std::string targetType{};
    bool useOccultation{};
    bool useMoonRadiation{};
    std::string panelingMoon{};
    std::string albedoDistributionMoon{};
    unsigned int numberOfPanelsMoon{};
    std::vector<int> numberOfPanelsPerRingMoon{};
    std::string thermalTypeMoon{};
    double stepSize{};
};

Settings loadSettings(char* path)
{
    Settings settings;

    std::ifstream f(path);
    nlohmann::json settings_ = nlohmann::json::parse(f);

    settings.id = settings_["id"];
    settings.saveDir = settings_["save_dir"];

    settings.simulationStart = settings_["simulation_start"];
    settings.simulationDuration = settings_["simulation_duration"];
    settings.simulationStartEpoch = spice_interface::convertDateStringToEphemerisTime(settings.simulationStart);
    settings.simulationEndEpoch = settings.simulationStartEpoch + settings.simulationDuration;
    settings.printInterval = settings.simulationDuration / 10;

    settings.targetType = settings_["target_type"];
    settings.useOccultation = settings_["use_occultation"];
    settings.useMoonRadiation = settings_["use_moon_radiation"];
    settings.panelingMoon = settings_["paneling_moon"];
    settings.albedoDistributionMoon = settings_["albedo_distribution_moon"];
    settings.numberOfPanelsMoon = settings_["number_of_panels_moon"];
    settings.numberOfPanelsPerRingMoon = settings_["number_of_panels_per_ring_moon"].get<std::vector<int>>();
    settings.thermalTypeMoon = settings_["thermal_type_moon"];

    settings.stepSize = settings_["step_size"];

    std::cout << settings << std::endl;

    return settings;
}

std::ostream& operator<<(std::ostream& os, const Settings& settings)
{
    os << std::boolalpha;
    os
        << "Settings:" << std::endl
        << " - ID: " << settings.id << std::endl
        << " - Save directory: "<<settings.saveDir<<std::endl
        << " - Simulation start: "<<formatEphemerisTime(settings.simulationStartEpoch)<<std::endl
        << " - Simulation duration: "<<(settings.simulationDuration / 60)<<" min"<<std::endl
        << " - Simulation end: "<<formatEphemerisTime(settings.simulationEndEpoch)<<std::endl
        << " - Target type: "<<settings.targetType<<std::endl
        << " - Use occultation: "<<settings.useOccultation<<std::endl
        << " - Use moon radiation: "<<settings.useMoonRadiation<<std::endl
        << " - Paneling type (Moon): "<<settings.panelingMoon<<std::endl
        << " - Albedo distribution (Moon): "<<settings.albedoDistributionMoon<<std::endl
        << " - Number of panels (Moon): "<<settings.numberOfPanelsMoon<<std::endl
        << " - Number of panels per ring (Moon): ";
    std::copy(settings.numberOfPanelsPerRingMoon.begin(), settings.numberOfPanelsPerRingMoon.end(),
              std::ostream_iterator<int>(os, ", "));
    os
        << std::endl
        << " - Thermal type (moon): "<<settings.thermalTypeMoon<<std::endl
        << " - Step size: " << settings.stepSize << " s" << std::endl;
    os << std::noboolalpha;
    return os;
}

std::string formatEphemerisTime(double et)
{
    SpiceChar out[65];
    timout_c(et, "YYYY-MM-DD HR:MN:SC UTC ::UTC", 65, out);
    return out;
}

namespace simulation_constants
{
    const auto globalFrameOrigin = "Moon";
    const auto globalFrameOrientation = "ECLIPJ2000";
    const auto moonFrame = "MOON_PA";

    const std::vector<std::string> bodiesToPropagate{"LRO"};
    const std::vector<std::string> centralBodies{"Moon"};
}


#endif //TUDATBUNDLE_LRO_JSON_H
