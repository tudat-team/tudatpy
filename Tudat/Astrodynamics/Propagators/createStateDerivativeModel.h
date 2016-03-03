#ifndef CREATESTATEDERIVATIVEMODEL_H
#define CREATESTATEDERIVATIVEMODEL_H

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

template< typename ScalarStateType, typename TimeType >
boost::shared_ptr< CentralBodyData< ScalarStateType, TimeType > > createCentralBodyData(
        const std::vector< std::string >& centralBodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    // Check whether the bodies that are to e integrated exist in bodyMap
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        if( bodyMap.count( bodiesToIntegrate.at( i ) ) == 0 )
        {
            std::cerr<<"Warning when creating CentralBodyData, body "<<bodiesToIntegrate.at( i )<<std::endl;
            std::cerr<<" not present in provided body map, cannot integrate body!"<<std::endl;
        }
    }

    std::vector< std::string > centralBodiesToUse;
    // Generate central body data for requested settings; if no central bodies provided, get inertial SSB as central body
    if( centralBodies.size( ) == 0 )
    {
        for( int i = 0; i < bodiesToIntegrate.size( ); i++ )
        {
            centralBodiesToUse.push_back( "SSB" );
        }
    }
    else
    {
        centralBodiesToUse = centralBodies;
    }

    std::map< std::string, boost::function< Eigen::Matrix< ScalarStateType, 6, 1 >( const TimeType ) > > bodyStateFunctions;

    for( unsigned int i = 0; i < centralBodiesToUse.size( ); i++ )
    {
        if( centralBodiesToUse.at( i )  != "SSB" )
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] =
                    boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris< ScalarStateType, TimeType >,
                                 bodyMap.at( centralBodiesToUse.at( i ) ), _1 );
        }
        else
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] = boost::lambda::constant( Eigen::Matrix< ScalarStateType, 6, 1 >::Zero( ) );
        }
    }
    return boost::make_shared< CentralBodyData< ScalarStateType, TimeType > >( centralBodiesToUse, bodiesToIntegrate, bodyStateFunctions );
}

template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createTranslationalStateDerivativeModel(
        const boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationPropagatorSettings,
        const  simulation_setup::NamedBodyMap& bodyMap, const TimeType startTime )
{
    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;\

    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData =
            createCentralBodyData< StateScalarType, TimeType >(
                translationPropagatorSettings->centralBodies_,
                translationPropagatorSettings->bodiesToIntegrate_,
                bodyMap );

    switch( translationPropagatorSettings->propagator_ )
    {
    case cowell:
    {
        stateDerivativeModel = boost::make_shared< NBodyCowellStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->accelerationsMap_, centralBodyData, translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    default:
        std::cerr<<"Error, did not recognize translational state propagation type: "<<translationPropagatorSettings->propagator_<<std::endl;
    }
    return stateDerivativeModel;
}

template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createStateDerivativeModel(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings, const  simulation_setup::NamedBodyMap& bodyMap, const TimeType startTime )
{
    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;
    switch( propagatorSettings->stateType_ )
    {
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationPropagatorSettings == NULL )
        {
            std::cerr<<"Error, expected translational state propagation settings when making state derivative model"<<std::endl;
        }
        else
        {
            stateDerivativeModel = createTranslationalStateDerivativeModel< StateScalarType, TimeType >(
                        translationPropagatorSettings, bodyMap, startTime );
        }
        break;
    }    
    default:
        std::cerr<<"Error, could not process state type "<<propagatorSettings->stateType_<<" when making state derivative model"<<std::endl;
    }
    return stateDerivativeModel;
}

template< typename StateScalarType = double, typename TimeType = double >
std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > createStateDerivativeModels(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings, const  simulation_setup::NamedBodyMap& bodyMap, const TimeType startTime )
{
    std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels;
    switch( propagatorSettings->stateType_ )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             propagatorIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                if( propagatorIterator->first != hybrid )
                {
                    stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                                         propagatorIterator->second.at( i ), bodyMap, startTime ) );
                }
                else
                {
                    std::cerr<<"Error when making state derivative model, cannot process nested hybrid propagators"<<std::endl;
                }
            }
        }
        break;
    }
    default:
        stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                             propagatorSettings, bodyMap, startTime ) );
    }

    return stateDerivativeModels;
}

}

}

#endif // CREATESTATEDERIVATIVEMODEL_H
