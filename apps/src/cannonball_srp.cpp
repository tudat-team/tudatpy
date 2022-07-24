#include <Eigen/Core>
#include <iostream>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"
#include "tudat/astro/electromagnetism/cannonBallRadiationPressureAcceleration.h"

using namespace tudat::electromagnetism;

int main()
{
    std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions{};

    auto radiationPressureInterface = std::make_shared<RadiationPressureInterface>(
            []() { return 3.839e26; },
            []() { return Eigen::Vector3d(0, 0, 0); },
            []() { return Eigen::Vector3d(tudat::physical_constants::ASTRONOMICAL_UNIT, 0, 0); },
            1,
            1,
            occultingBodyPositions,
            std::vector< double >{}
            );

    auto radiationPressureAcceleration = std::make_shared< CannonBallRadiationPressureAcceleration >(
            radiationPressureInterface->getSourcePositionFunction(),
            radiationPressureInterface->getTargetPositionFunction(),
            std::bind( &RadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface ),
            std::bind( &RadiationPressureInterface::getRadiationPressureCoefficient, radiationPressureInterface ),
            std::bind( &RadiationPressureInterface::getArea, radiationPressureInterface ),
            []() { return 2; }
            );

    radiationPressureInterface->updateInterface();
    radiationPressureAcceleration->updateMembers();

    std::cout << radiationPressureInterface->getCurrentRadiationPressure() << " N/m^2" << std::endl;
    std::cout << radiationPressureAcceleration->getAcceleration() << " m/s^2" << std::endl;

    return 0;
}
