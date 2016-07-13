#ifndef CREATEMASSRATEMODELS_H
#define CREATEMASSRATEMODELS_H


#include <vector>
#include <string>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

class MassRateModelSettings
{
public:

    MassRateModelSettings(
            const basic_astrodynamics::AvailableMassRateModels massRateType ):
        massRateType_( massRateType ){ }

    basic_astrodynamics::AvailableMassRateModels massRateType_;

};

class CustomMassRateModelSettings: public MassRateModelSettings
{
public:

    CusromMassRateModelSettings(
            const boost::function< double( const double ) > massRateFunction ):
        MassRateModelSettings( basic_astrodynamics::custom_mass_rate_model ),
    massRateFunction_( massRateFunction ){ }

    boost::function< double( const double ) > massRateFunction_;

};

class FromThrustMassModelSettings: public MassRateModelSettings
{
public:

    FromThrustMassModelSettings(
            const bool useAllThrustModels = 1,
            const std::string& associatedThroustSource = "" ):
        MassRateModelSettings( basic_astrodynamics::from_thrust_mass_rate_model ),
    associatedThroustSource_( associatedThroustSource ), useAllThrustModels_( useAllThrustModels ){ }

    std::string associatedThroustSource_;

    bool useAllThrustModels_;

};


boost::shared_ptr< basic_astrodynamics::MassRateModel >
createMassRateModel(
        const std::string& bodyWithMassRate,
        const boost::shared_ptr< MassRateModelSettings > massRateModelSettings,
        const NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModels = basic_astrodynamics::AccelerationMap( ) );


std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > > createMassRateModelsMap(
        const NamedBodyMap& bodyMap,
        const std::map< std::string, std::vector< boost::shared_ptr< MassRateModelSettings > > >& massRateModelSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModels = AccelerationMap( ) );

} // namespace simulation_setup

} // namespace tudat

#endif // CREATEMASSRATEMODELS_H
