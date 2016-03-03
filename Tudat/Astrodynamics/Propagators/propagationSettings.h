#ifndef PROPAGATIONSETTINGS_H
#define PROPAGATIONSETTINGS_H

#include <vector>
#include <string>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Basics/utilities.h"

namespace tudat
{

namespace propagators
{

enum IntegratedStateType
{
    hybrid,
    transational_state,
    rotational_state,
    proper_time
};

int getSingleIntegrationSize( const IntegratedStateType stateType );

int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType );

enum TranslationalPropagatorType
{
    cowell = 0
};


template< typename StateScalarType >
class PropagatorSettings
{
public:
    PropagatorSettings( const IntegratedStateType stateType,
                        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates ):
        stateType_( stateType ), initialStates_( initialBodyStates ), singleArcStateSize_( initialBodyStates.rows( ) ),
        bodyArcStartAndEndTimes_( std::vector< std::pair< double, double > >( ) ), isPropagatorMultiArc_( 0 ){ }

    virtual ~PropagatorSettings( ){ }

    IntegratedStateType stateType_;

    double getNumberOfIntegrationArcs( )
    {
        return bodyArcStartAndEndTimes_.size( );
    }

    virtual void setIntegrationArcs( const std::vector< std::pair< double, double > >& bodyArcStartAndEndTimes )
    {
        bodyArcStartAndEndTimes_ = bodyArcStartAndEndTimes;
        isPropagatorMultiArc_ = 1;
        singleArcStateSize_ = initialStates_.rows( ) / bodyArcStartAndEndTimes_.size( );
    }

    bool isPropagatorMultiArc_;

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( )
    {
        return initialStates_;
    }

    virtual void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        initialStates_ = initialBodyStates;
        if( bodyArcStartAndEndTimes_.size( ) > 0 )
        {
            singleArcStateSize_ = initialStates_.rows( ) / bodyArcStartAndEndTimes_.size( );
        }
        else
        {
            singleArcStateSize_ = initialStates_.rows( );
        }
    }


    int getSingleArcStateSize( )
    {
        return singleArcStateSize_;
    }

    std::vector< std::pair< double, double > > getBodyArcStartAndEndTimes( )
    {
        return bodyArcStartAndEndTimes_;
    }

protected:
    std::vector< std::pair< double, double > > bodyArcStartAndEndTimes_;

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    int singleArcStateSize_;

};

template< typename StateScalarType >
int getMultiTypePropagatorSingleArcStateSize(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorSettingsMap )
{
    int singleArcStateSize = 0;
    for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsMap.begin( ); typeIterator != propagatorSettingsMap.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            singleArcStateSize += typeIterator->second.at( i )->getSingleArcStateSize( );
        }
    }
    return singleArcStateSize;
}

template< typename StateScalarType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > createCombinedInitialState(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >& propagatorSettingsList )
{
    int totalSize = getMultiTypePropagatorSingleArcStateSize( propagatorSettingsList );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > combinedInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalSize, 1 );

    int currentIndex = 0;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState;
    for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            currentInitialState = typeIterator->second.at( i )->getInitialStates( );
            combinedInitialState.segment( currentIndex, currentInitialState.rows( ) ) = currentInitialState;
            currentIndex += currentInitialState.rows( );
        }
    }
    return combinedInitialState;
}

template< typename StateScalarType >
std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > createPerTypePropagatorSettingsMap(
        const std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > >& propagatorsSettingsList )
{
    std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorsSettingsMap;
    for( unsigned int i = 0; i < propagatorsSettingsList.size( ); i++ )
    {
        propagatorsSettingsMap[ propagatorsSettingsList.at( i )->stateType_ ].push_back( propagatorsSettingsList.at( i ) );
    }
    return propagatorsSettingsMap;
}


template< typename StateScalarType >
class MultiTypePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:
    using PropagatorSettings< StateScalarType >::singleArcStateSize_;
    using PropagatorSettings< StateScalarType >::bodyArcStartAndEndTimes_;

    MultiTypePropagatorSettings(
            const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorSettingsMap ):
        PropagatorSettings< StateScalarType >( hybrid, createCombinedInitialState< StateScalarType >( propagatorSettingsMap ) ),
        propagatorSettingsMap_( propagatorSettingsMap )
    {
        singleArcStateSize_ = getMultiTypePropagatorSingleArcStateSize( propagatorSettingsMap_ );
    }


    ~MultiTypePropagatorSettings( ){ }

    void setIntegrationArcs( const std::vector< std::pair< double, double > >& bodyArcStartAndEndTimes )
    {
        bodyArcStartAndEndTimes_ = bodyArcStartAndEndTimes;
        this->isPropagatorMultiArc_ = 1;
        singleArcStateSize_ = this->initialStates_.rows( ) / bodyArcStartAndEndTimes_.size( );

        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             propagatorIterator = propagatorSettingsMap_.begin( );
             propagatorIterator!= propagatorSettingsMap_.end( ); propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                propagatorIterator->second.at( i )->setIntegrationArcs( bodyArcStartAndEndTimes );
            }
        }
    }

    void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        int currentStartIndex = 0;
        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = propagatorSettingsMap_.begin( ); propagatorIterator != propagatorSettingsMap_.end( );
             propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                int currentParameterSize = propagatorIterator->second.at( i )->getInitialStates( ).rows( );

                if( currentParameterSize + currentStartIndex > initialBodyStates.rows( ) )
                {
                    std::cerr<<"Error when resetting multi-arc state, sizes are incompatible "<<std::endl;
                }
                propagatorIterator->second.at( i )->resetInitialStates(
                            initialBodyStates.block( currentStartIndex, 0, currentParameterSize, 1 ) );
                currentStartIndex += currentParameterSize;

            }
        }

        if( currentStartIndex != initialBodyStates.rows( ) )
        {
            std::cerr<<"Error when resetting multi-arc state, sizes are incompatible B "<<
                       currentStartIndex<<" "<<initialBodyStates.rows( )<<std::endl;
        }

        updateInitialStates( );
    }

    void updateInitialStates( )
    {
        this->initialStates_ = createCombinedInitialState< StateScalarType >( propagatorSettingsMap_ );
    }

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorSettingsMap_;

};

//! Settings for propagator to use in numerical integration of equations of motion.
/*!
 *  Settings for propagator to use in numerical integration of equations of motion. The propagator defines the form of the equations of
 *  motion (i.e. Cowell, Encke, Gauss etc.). This base class can be used for Cowell propagator. Other propagators have dedicated derived class.
 */
template< typename StateScalarType = double >
class TranslationalStatePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor for user-defined central bodies
    /*!
     *  Constructor for user-defined central bodies.
     *  \param numberOfBodies Central bodies to use for each of the numerically integrated bodies (in same order as bodiesToIntegrate_)
     *  in DynamicsSimulator
     *  \param propagator Type of propagator to use.
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const TranslationalPropagatorType propagator = cowell):
        PropagatorSettings< StateScalarType >( transational_state, initialBodyStates ), centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    ~TranslationalStatePropagatorSettings( ){ }

    //! Type of propagator to use.
    /*!
     *  Type of propagator to use.
     */
    TranslationalPropagatorType propagator_;

    AccelerationMap accelerationsMap_;

    std::vector< std::string > bodiesToIntegrate_;

    //! Central bodies to use for each of the numerically integrated bodies.
    /*!
     *  Central bodies to use for each of the numerically integrated bodies (in same order as bodiesToIntegrate_)
     *  in DynamicsSimulator
     */
    std::vector< std::string > centralBodies_;
};

template< typename StateScalarType >
std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStateList;
    switch( propagatorSettings->stateType_ )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > singleTypeIntegratedStateList;

        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            if( typeIterator->first != hybrid )
            {
                for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
                {
                    singleTypeIntegratedStateList = getIntegratedTypeAndBodyList< StateScalarType >( typeIterator->second.at( i ) );

                    if( singleTypeIntegratedStateList.begin( )->first != typeIterator->first || singleTypeIntegratedStateList.size( ) != 1 )
                    {
                        std::cerr<<"Error when making integrated state list for hybrid propagator, inconsistency encountered "<<
                                   singleTypeIntegratedStateList.begin( )->first<<" "<<typeIterator->first<<" "<<
                                   singleTypeIntegratedStateList.size( )<<" "<<singleTypeIntegratedStateList.begin( )->second.size( )<<std::endl;
                    }
                    else
                    {
                        for( unsigned int j = 0; j < singleTypeIntegratedStateList[ typeIterator->first ].size( ); j++ )
                        {
                            integratedStateList[ typeIterator->first ].push_back( singleTypeIntegratedStateList.begin( )->second.at( j ) );
                        }

                    }
                }
            }
            else
            {
                std::cerr<<"Error when making integrated state list, cannot handle hybrid propagator inside hybrid propagator"<<std::endl;
            }
        }

        break;
    }
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < translationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair( translationalPropagatorSettings->bodiesToIntegrate_.at( i ), "" ) );
        }

        integratedStateList[ transational_state ] = integratedBodies;

        break;
    }
    default:
        std::cerr<<"Error, could not process integrated state type "<<propagatorSettings->stateType_<<std::endl;
    }

    return integratedStateList;
}

}

}

#endif // PROPAGATIONSETTINGS_H
