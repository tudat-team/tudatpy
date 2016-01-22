#ifndef CREATERADIATIONPRESSUREINTERFACE_H
#define CREATERADIATIONPRESSUREINTERFACE_H

#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace simulation_setup
{

static const std::map< std::string, double > defaultRadiatedPowerValues =
{ { "Sun",  3.839E26 } };

enum RadiationPressureType
{
    cannon_ball
};

class RadiationPressureInterfaceSettings
{
public:
    RadiationPressureInterfaceSettings(
            const RadiationPressureType radiationPressureType,
            const std::string& sourceBody,
            const std::vector< std::string > occultingBodies = std::vector< std::string >( ) ):
        radiationPressureType_( radiationPressureType ), sourceBody_( sourceBody ),
        occultingBodies_( occultingBodies ){  }

    virtual ~RadiationPressureInterfaceSettings( ){ }

    RadiationPressureType getRadiationPressureType( ){ return radiationPressureType_; }

    std::string getSourceBody( ){ return sourceBody_; }

    std::vector< std::string > getOccultingBodies( ){ return occultingBodies_; }

protected:
    RadiationPressureType radiationPressureType_;

    std::string sourceBody_;

    std::vector< std::string > occultingBodies_;
};

class CannonBallRadiationPressureInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:
    CannonBallRadiationPressureInterfaceSettings(
            const std::string& sourceBody, const double area, const double radiationPressureCoefficient,
            const std::vector< std::string >& occultingBodies = std::vector< std::string >( ) ):
        RadiationPressureInterfaceSettings( cannon_ball, sourceBody, occultingBodies ),
        area_( area ), radiationPressureCoefficient_( radiationPressureCoefficient ){ }

    double getArea( ){ return area_; }

    double getRadiationPressureCoefficient( ){ return radiationPressureCoefficient_; }


private:
    double area_;

    double radiationPressureCoefficient_;
};

void getOccultingBodiesInformation(
        const NamedBodyMap& bodyMap, const std::vector< std::string >& occultingBodies,
        std::vector< boost::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii );

boost::shared_ptr< electro_magnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const NamedBodyMap& bodyMap );



}

}

#endif // CREATERADIATIONPRESSUREINTERFACE_H
