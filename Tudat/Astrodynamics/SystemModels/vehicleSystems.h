#ifndef VEHICLESYSTEMS_H
#define VEHICLESYSTEMS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace system_models
{

class VehicleSystems
{

    virtual double getCurrentMass( ) = 0;

};

class SingleStageVehicleSystems
{
public:

    double updateMass( const double mass )
    {
        if( dryMass_ > mass )
        {

        }
        else
        {
            propellantMass_ = mass - dryMass_;
        }
    }

    double getCurrentMass( )
    {
        return propellantMass_ + dryMass_;
    }

private:

    boost::shared_ptr< EngineModel > engineModel_;

    double propellantMass_;

    double dryMass_;

};

class MutltiStageVehicleSystems
{
public:
    double getCurrentMass( )
    {

    }

private:

    std::vector< boost::shared_ptr< SingleStageVehicleSystems > > singleStageModels_;

    int activeStage_;
};

}

}

#endif // VEHICLESYSTEMS_H
