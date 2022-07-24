#include <Eigen/Core>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <tudat/simulation/simulation.h>
#include "tudat/interface/spice/spiceEphemeris.h"

using namespace tudat;

void loadLROSpiceKernels();
void saveLROEphemeris(double, std::shared_ptr<ephemerides::SpiceEphemeris>);

int main()
{
    loadLROSpiceKernels();

    auto coverageSpiceLROStart = basic_astrodynamics::convertCalendarDateToJulianDay(
            2009, 06, 18, 22, 17, 06.184);

    auto ephemerisLRO = std::make_shared<ephemerides::SpiceEphemeris>(
            "LRO", "Moon", true, true, false, "ECLIPJ2000", coverageSpiceLROStart
            );
    saveLROEphemeris(coverageSpiceLROStart, ephemerisLRO);

//    std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormalFunctions{
//            []( const double ) { return Eigen::Vector3d(-1, 0, 0); }
//    };
//    std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions{};
//
//    auto radiationPressureInterface = std::make_shared<PanelledRadiationPressureInterface >(
//                    []() { return 1; },
//                    []() { return Eigen::Vector3d(0, 0, 0); },
//                    []() { return Eigen::Vector3d(1, 0, 0); },
//                    localFrameSurfaceNormalFunctions,
//                    std::vector< double >{1},
//                    std::vector< double >{1},
//                    std::vector< double >{0},
//                    []() { return Eigen::Quaterniond( Eigen::Matrix3d::Identity()); },
//                    occultingBodyPositions,
//                    std::vector< double >{});
//
//    auto radiationPressureAcceleration = std::make_shared<PanelledRadiationPressureAcceleration>(
//            radiationPressureInterface,
//            []() { return 1; }
//            );
//
//    radiationPressureInterface->updateInterface();
//    radiationPressureAcceleration->updateMembers();
//
//    std::cout << radiationPressureAcceleration->getAcceleration() << " m/s^2" << std::endl;


    return 0;
}

void loadLROSpiceKernels()
{
    std::string path = "/home/dominik/dev/tudat-bundle/spice/lro/data/spk";
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
        if (entry.path().extension() == ".bsp")
        {
            spice_interface::loadSpiceKernelInTudat(entry.path().string());
        }
    }
}

void saveLROEphemeris(double coverageSpiceLROStart, std::shared_ptr<ephemerides::SpiceEphemeris> ephemeris)
{
    std::map<double, Eigen::Vector3d> positions;
    for (int t = 1 * 24 * 3600; t < 10 * 365 * 24 * 3600; t += 60*60)
    {
        auto timestampInJD = coverageSpiceLROStart + t / physical_constants::JULIAN_DAY;
        positions[timestampInJD] = ephemeris->getCartesianPosition(t);
    }

    tudat::input_output::writeDataMapToTextFile( positions,
                                                 "lroPositionHourly.dat",
                                                 boost::filesystem::current_path(),
                                                 "",
                                                 std::numeric_limits< double >::digits10,
                                                 std::numeric_limits< double >::digits10,
                                                 "," );
}
