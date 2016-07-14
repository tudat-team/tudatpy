#ifndef VEHICLESYSTEMS_H
#define VEHICLESYSTEMS_H

#include <map>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace system_models
{

class VehicleSystems
{
public:

    VehicleSystems( const double dryMass = TUDAT_NAN ):
        dryMass_( dryMass ){ }

    void updateMass( const double mass )
    {
        if( dryMass_ > mass )
        {
            throw std::runtime_error( "Error when setting vehicle mass, dry mass is larger than actual mass" );
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

    std::map< std::string, boost::shared_ptr< EngineModel > > getEngineModels( )
    {
        return engineModels_;
    }

    void setEngineModel(
            const boost::shared_ptr< EngineModel > engineModel, const std::string engineName = "" )
    {
        if( engineModels_.count( engineName ) )
        {
            std::cerr<<"Warning, engine model of name "<<engineModel<<" already exists, overriding old model"<<std::endl;
        }

        engineModels_[ engineName ] = engineModel;
    }

private:

    std::map< std::string, boost::shared_ptr< EngineModel > > engineModels_;

    double propellantMass_;

    double dryMass_;

};

//class MutltiStageVehicleSystems
//{
//public:
//    double getCurrentMass( )
//    {

//    }

//private:

//    std::vector< boost::shared_ptr< SingleStageVehicleSystems > > singleStageModels_;

//    int activeStage_;
//};

}

}

#endif // VEHICLESYSTEMS_H
