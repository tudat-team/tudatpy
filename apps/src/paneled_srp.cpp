#include <Eigen/Core>
#include <iostream>

#include "tudat/astro/electromagnetism/panelledRadiationPressure.h"

using namespace tudat::electromagnetism;

int main()
{
//    tudat::spice_interface::loadStandardSpiceKernels();
//
//    auto ephemericLRO = std::make_shared<tudat::ephemerides::SpiceEphemeris>(
//            "LRO", "Earth"
//            );
//    std::cout << ephemericLRO->getCartesianPosition(0) << std::endl;

    std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormalFunctions{
            []( const double ) { return Eigen::Vector3d(-1, 0, 0); }
    };
    std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions{};

    auto radiationPressureInterface = std::make_shared<PanelledRadiationPressureInterface >(
                    []() { return 1; },
                    []() { return Eigen::Vector3d(0, 0, 0); },
                    []() { return Eigen::Vector3d(1, 0, 0); },
                    localFrameSurfaceNormalFunctions,
                    std::vector< double >{1},
                    std::vector< double >{1},
                    std::vector< double >{0},
                    []() { return Eigen::Quaterniond( Eigen::Matrix3d::Identity()); },
                    occultingBodyPositions,
                    std::vector< double >{});

    auto radiationPressureAcceleration = std::make_shared<PanelledRadiationPressureAcceleration>(
            radiationPressureInterface,
            []() { return 1; }
            );

    radiationPressureInterface->updateInterface();
    radiationPressureAcceleration->updateMembers();

    std::cout << radiationPressureAcceleration->getAcceleration() << " m/s^2" << std::endl;

    return 0;
}
