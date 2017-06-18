#ifndef CREATETORQUEMODEL_H
#define CREATETORQUEMODEL_H


#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr<
basic_astrodynamics::TorqueSettings > > > > SelectedTorqueMap;

boost::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const boost::shared_ptr< basic_astrodynamics::TorqueSettings > torqueType,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque );

basic_astrodynamics::TorqueModelMap createTorqueModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedTorqueMap& selectedTorquePerBody );

}

}

#endif // CREATETORQUEMODEL_H
