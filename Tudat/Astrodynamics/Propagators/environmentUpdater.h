#ifndef ENVIRONMENTUPDATER_H
#define ENVIRONMENTUPDATER_H

#include <vector>
#include <string>
#include <map>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

enum EnvironmentModelsToUpdate
{
    body_transational_state_update = 0,
    body_rotational_state_update = 1,
    body_mass_update = 2,
    vehicle_flight_conditions_update = 3,
    radiation_pressure_interface_update = 4,
};

template< typename StateScalarType, typename TimeType >
class EnvironmentUpdater
{
public:

    EnvironmentUpdater( const simulation_setup::NamedBodyMap& bodyList,
                        const std::map< EnvironmentModelsToUpdate, std::vector< std::string > >& updateSettings,
                        const std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > >& integratedStates =
            ( std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > >( ) ) ):
        bodyList_( bodyList ), integratedStates_( integratedStates )
    {
        // Set update function to be evaluated as dependent variables of state and time during each integration time step.
        setUpdateFunctions( updateSettings );
    }

    void updateEnvironment(
            const TimeType currentTime,
            const std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& integratedStatesToSet,
            const std::vector< IntegratedStateType >& setIntegratedStatesFromEnvironment = std::vector< IntegratedStateType >( ) )
    {
        if( integratedStatesToSet.size( ) + setIntegratedStatesFromEnvironment.size( ) != integratedStates_.size( ) )
        {
            std::cerr<<"Error when updating environment, input size is inconsistent "<<
                       integratedStatesToSet.size( )<<" "<<setIntegratedStatesFromEnvironment.size( )<<" "<<integratedStates_.size( )<<std::endl;
        }

        setIntegratedStatesInEnvironment( integratedStatesToSet );

        setStatesFromEnvironment( setIntegratedStatesFromEnvironment, currentTime );

        for( outerCurrentStateFromEnvironmentIterator_ = currentStateFromEnvironmentList_.begin( );
             outerCurrentStateFromEnvironmentIterator_ != currentStateFromEnvironmentList_.end( );
             outerCurrentStateFromEnvironmentIterator_++ )
        {
            for( currentStateFromEnvironmentIterator_ = outerCurrentStateFromEnvironmentIterator_->second.begin( );
                 currentStateFromEnvironmentIterator_ != outerCurrentStateFromEnvironmentIterator_->second.end( );
                 currentStateFromEnvironmentIterator_++ )
            {
                currentStateFromEnvironmentIterator_->second( currentTime );
            }
        }

        // Evaluate update functions (dependent variables of state and time) determined by setUpdateFunctions
        for( updateFunctionIterator = updateFunctionList_.begin( ); updateFunctionIterator != updateFunctionList_.end( );
             updateFunctionIterator++ )
        {
            for( unsigned int i = 0; i < updateFunctionIterator->second.size( ); i++ )
            {
                updateFunctionIterator->second.at( i ).second( );
            }
        }

        // Evaluate time-dependent update functions (dependent variables of state and time) determined by setUpdateFunctions
        for( updateTimeIterator = updateTimeFunctionList_.begin( ); updateTimeIterator != updateTimeFunctionList_.end( );
             updateTimeIterator++ )
        {
            for( unsigned int i = 0; i < updateTimeIterator->second.size( ); i++ )
            {
                updateTimeIterator->second.at( i ).second( static_cast< double >( currentTime ) );
            }
        }
    }

private:

    void setIntegratedStatesInEnvironment(
            const std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& integratedStatesToSet )
    {
        for( typename std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator integratedStateIterator =
             integratedStatesToSet.begin( ); integratedStateIterator != integratedStatesToSet.end( ); integratedStateIterator++ )
        {
            switch( integratedStateIterator->first )
            {
            case transational_state:
            {
                std::vector< std::pair< std::string, std::string > > bodiesWithIntegratedStates =
                        integratedStates_.at( transational_state );
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_[ bodiesWithIntegratedStates[ i ].first ]->template setTemplatedState< StateScalarType >(
                                integratedStateIterator->second.segment( i * 6, 6 ) );
                }
                break;
            };
            default:
                std::cerr<<"Error, could not find integrated state settings for "<<integratedStateIterator->first<<std::endl;
            }
        }
    }

    void setStatesFromEnvironment(
            const std::vector< IntegratedStateType >& statesToSet,
            const TimeType currentTime )
    {
        for( unsigned int i = 0; i < statesToSet.size( ); i++ )
        {
            switch( statesToSet.at( i ) )
            {
            case transational_state:
            {
                std::vector< std::pair< std::string, std::string > > bodiesWithIntegratedStates =
                        integratedStates_.at( transational_state );
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_[ bodiesWithIntegratedStates[ i ].first ]->template setTemplatedStateFromEphemeris< StateScalarType, TimeType >(
                                currentTime );

                }
                break;
            }
            default:
                std::cerr<<"Error, could not find  state settings for "<<statesToSet.at( i )<<std::endl;
            }
        }
    }

    void setTranslationalStateUpdateFunctions( std::vector< std::string > bodiesWithIntegratedState )
    {
        std::vector< std::string >::iterator ephemerisBodyLookup;

        // Iterate over all bodies.
        for( simulation_setup::NamedBodyMap::iterator bodyIterator = bodyList_.begin( ); bodyIterator != bodyList_.end( ); bodyIterator++ )
        {
            ephemerisBodyLookup =
                    std::find( bodiesWithIntegratedState.begin( ), bodiesWithIntegratedState.end( ), bodyIterator->first );
            if( ephemerisBodyLookup == bodiesWithIntegratedState.end( ) )
            {
                boost::function< void( const TimeType ) > stateSetFunction =
                        boost::bind( &simulation_setup::Body::setTemplatedStateFromEphemeris< StateScalarType, TimeType >, bodyIterator->second, _1 );
                currentStateFromEnvironmentList_[ body_transational_state_update ].insert( std::make_pair( bodyIterator->first, stateSetFunction ) );
            }
        }
    }

    void setRotationalStateUpdateFunctions( std::vector< std::string > bodiesWithIntegratedState )
    {
        std::vector< std::string >::iterator ephemerisBodyLookup;

        // Iterate over all bodies.
        for( simulation_setup::NamedBodyMap::iterator bodyIterator = bodyList_.begin( ); bodyIterator != bodyList_.end( ); bodyIterator++ )
        {
            ephemerisBodyLookup =
                    std::find( bodiesWithIntegratedState.begin( ), bodiesWithIntegratedState.end( ), bodyIterator->first );
            if( ephemerisBodyLookup == bodiesWithIntegratedState.end( ) )
            {
                boost::function< void( const TimeType ) > stateSetFunction =
                        boost::bind( &simulation_setup::Body::setCurrentRotationalStateToLocalFrameFromEphemeris, bodyIterator->second, _1 );

                currentStateFromEnvironmentList_[ body_rotational_state_update ].insert( std::make_pair( bodyIterator->first, stateSetFunction ) );
            }
        }
    }

    void setUpdateFunctions( const std::map< EnvironmentModelsToUpdate, std::vector< std::string > >& updateSettings )
    {
        if( integratedStates_.count( transational_state ) > 0 )
        {
            std::vector< std::string > integratedBodies;
            for( unsigned int i = 0; i < integratedStates_.at( transational_state ).size( ); i++ )
            {
                integratedBodies.push_back( integratedStates_.at( transational_state ).at( i ).first );
            }
            setTranslationalStateUpdateFunctions( integratedBodies );
        }

        for( std::map< EnvironmentModelsToUpdate, std::vector< std::string > >::const_iterator updateIterator =
             updateSettings.begin( ); updateIterator != updateSettings.end( ); updateIterator++ )
        {
            std::vector< std::string > currentBodies = updateIterator->second;
            for( unsigned int i = 0; i < currentBodies.size( ); i++ )
            {
                if( currentBodies.at( i ) != "" )
                {
                    if( bodyList_.count( currentBodies.at( i ) ) == 0 )
                    {
                        std::cerr<<"Error when setting environment update functions, could not find body "<<currentBodies.at( i )<<std::endl;
                    }
                    switch( updateIterator->first )
                    {
                    case body_transational_state_update:
                    {
                        bool addUpdate = 1;
                        if( integratedStates_.count( transational_state ) > 0 )
                        {
                            std::pair< std::string, std::string > bodyToCheck = std::make_pair( currentBodies.at( i ), "" );
                            std::vector< std::pair< std::string, std::string > > integratedTranslationalStates =
                                    integratedStates_.at( transational_state );
                            if( std::find( integratedTranslationalStates.begin( ), integratedTranslationalStates.end( ), bodyToCheck ) !=
                                    integratedTranslationalStates.end( ) )
                            {
                                addUpdate = 0;
                            }
                        }
                        break;
                    }
                    case body_rotational_state_update:
                    {

                        boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris =
                                bodyList_.at( currentBodies.at( i ) )->getRotationalEphemeris( );
                        if( rotationalEphemeris != NULL )
                        {
                            boost::function< void( const TimeType ) > rotationalStateSetFunction =
                                    boost::bind( &simulation_setup::Body::setCurrentRotationalStateToLocalFrameFromEphemeris,
                                                 bodyList_.at( currentBodies.at( i ) ), _1 );
                            currentStateFromEnvironmentList_[ body_rotational_state_update ].insert( std::make_pair( currentBodies.at( i ), rotationalStateSetFunction ) );
                        }

                        break;

                    }
                    case body_mass_update:
                    {
                        updateTimeFunctionList_[ body_mass_update ].push_back(
                                    std::make_pair( currentBodies.at( i ),
                                                    boost::bind( &simulation_setup::Body::updateMass, bodyList_.at( currentBodies.at( i ) ), _1  ) ) );
                        break;
                    }
                    case vehicle_flight_conditions_update:
                    {
                        // Check if current body has flight conditions set.
                        if( bodyList_.at( currentBodies.at( i ) )->getFlightConditions( ) != NULL )
                        {
                            // If vehicle has flight conditions, add flight conditions update function to update list.
                            updateFunctionList_[ vehicle_flight_conditions_update ].push_back(
                                        std::make_pair(
                                            currentBodies.at( i ), boost::bind( &aerodynamics::FlightConditions::updateConditions,
                                                                                bodyList_.at( currentBodies.at( i ) )->getFlightConditions( ) ) ) );
                        }
                        else
                        {
                            std::cerr<<"Request flight condition update of "<<currentBodies.at( i )<<", but body ihas no flight conditions"<<std::endl;
                        }
                        break;
                    }

                    case radiation_pressure_interface_update:
                    {
                        // Get body radiation pressure interface(s) (one per source)
                        std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > > radiationPressureInterfaces =
                                bodyList_.at( currentBodies.at( i ) )->getRadiationPressureInterfaces( );

                        if( radiationPressureInterfaces.size( ) == 0 )
                        {
                            std::cerr<<"Request radiation pressure update of "<<currentBodies.at( i )<<", but body has no radiation pressure interfaces"<<std::endl;
                        }
                        else if( radiationPressureInterfaces.size( ) > 1 )
                        {
                            std::cerr<<"Request radiation pressure update of "<<currentBodies.at( i )<<
                                       ", but body has multiple radiation pressure interfaces: updating all."<<std::endl;
                        }

                        // Add each interface update function to update list.
                        for( std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > > ::iterator iterator =
                             radiationPressureInterfaces.begin( ); iterator != radiationPressureInterfaces.end( ); iterator++ )
                        {
                            updateTimeFunctionList_[ radiation_pressure_interface_update ].push_back(
                                        std::make_pair( currentBodies.at( i ),
                                                        boost::bind( &electro_magnetism::RadiationPressureInterface::updateInterface,
                                                                     iterator->second, _1 ) ) );
                        }
                        break;
                    }
                    }
                }
            }
        }
    }


    simulation_setup::NamedBodyMap bodyList_;


    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStates_;


    std::map< EnvironmentModelsToUpdate, std::multimap< std::string, boost::function< void( const TimeType ) > > > currentStateFromEnvironmentList_;

    typename std::map< EnvironmentModelsToUpdate, std::multimap< std::string, boost::function< void( const TimeType ) > > >::iterator outerCurrentStateFromEnvironmentIterator_;

    typename std::multimap< std::string, boost::function< void( const TimeType ) > >::iterator currentStateFromEnvironmentIterator_;


    std::map< EnvironmentModelsToUpdate, std::vector< std::pair< std::string, boost::function< void( ) > > > > updateFunctionList_;

    std::map< EnvironmentModelsToUpdate, std::vector< std::pair< std::string, boost::function< void( ) > > > >::iterator updateFunctionIterator;


    std::map< EnvironmentModelsToUpdate, std::vector< std::pair< std::string, boost::function< void( const double ) > > > > updateTimeFunctionList_;

    std::map< EnvironmentModelsToUpdate, std::vector< std::pair< std::string, boost::function< void( const double ) > > > >::iterator updateTimeIterator;


};

}

}

#endif // ENVIRONMENTUPDATER_H
