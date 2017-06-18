#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< basic_astrodynamics::TorqueModel > createTorqueModel(
        const boost::shared_ptr< simulation_setup::Body > bodyUndergoingTorque,
        const boost::shared_ptr< simulation_setup::Body > bodyExertingTorque,
        const boost::shared_ptr< basic_astrodynamics::TorqueSettings > torqueType,
        const std::string& nameOfBodyUndergoingTorque,
        const std::string& nameOfBodyExertingTorque )
{
    boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel;

    switch( torqueType->torqueType_ )
    {
    default:
        std::cerr<<"Error, did not recognize type "<<torqueType->torqueType_<<" when making torque model"<<std::endl;
    }

    return torqueModel;
}

basic_astrodynamics::TorqueModelMap createTorqueModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedTorqueMap& selectedTorquePerBody )
{
    basic_astrodynamics::TorqueModelMap torqueModelMap;

    for( SelectedTorqueMap::const_iterator acceleratedBodyIterator = selectedTorquePerBody.begin( );
         acceleratedBodyIterator != selectedTorquePerBody.end( ); acceleratedBodyIterator++ )
    {
        if( bodyMap.count( acceleratedBodyIterator->first ) == 0 )
        {
            std::cerr<<"Error A, could not find body "<<acceleratedBodyIterator->first<<" when making torque model map."<<std::endl;
        }
        else
        {
            for( std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::TorqueSettings > > >::const_iterator
                 acceleratingBodyIterator = acceleratedBodyIterator->second.begin( );
                 acceleratingBodyIterator != acceleratedBodyIterator->second.end( ); acceleratingBodyIterator++ )
            {
                if( bodyMap.count( acceleratingBodyIterator->first ) == 0 )
                {
                    std::cerr<<"Error B, could not find body "<<acceleratingBodyIterator->first<<" when making torque model map."<<std::endl;
                }
                for( unsigned int i = 0; i < acceleratingBodyIterator->second.size( ); i++ )
                {
                    torqueModelMap[ acceleratedBodyIterator->first ][ acceleratingBodyIterator->first ].push_back(
                                createTorqueModel(
                                bodyMap.at( acceleratedBodyIterator->first ), bodyMap.at( acceleratingBodyIterator->first ),
                                acceleratingBodyIterator->second.at( i ),
                                acceleratedBodyIterator->first, acceleratingBodyIterator->first ) );
                }
            }
        }
    }

    return torqueModelMap;
}

}

}

