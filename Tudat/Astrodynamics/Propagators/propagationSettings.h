#ifndef PROPAGATIONSETTINGS_H
#define PROPAGATIONSETTINGS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{

namespace propagators
{

enum IntegratedStateType
{
    transational_state
};

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
        stateType_( stateType ), initialStates_( initialBodyStates ){ }

    virtual ~PropagatorSettings( ){ }

    IntegratedStateType stateType_;

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( )
    {
        return initialStates_;
    }

    virtual void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        initialStates_ = initialBodyStates;
    }

protected:
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

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


    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
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

    //! Central bodies to use for each of the numerically integrated bodies.
    /*!
     *  Central bodies to use for each of the numerically integrated bodies (in same order as bodiesToIntegrate_)
     *  in DynamicsSimulator
     */
    std::vector< std::string > centralBodies_;

    basic_astrodynamics::AccelerationMap accelerationsMap_;

    std::vector< std::string > bodiesToIntegrate_;


    //! Type of propagator to use.
    /*!
     *  Type of propagator to use.
     */
    TranslationalPropagatorType propagator_;

};

template< typename StateScalarType >
std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStateList;
    switch( propagatorSettings->stateType_ )
    {    
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
        throw std::runtime_error( "Error, could not process integrated state type " +
                                  boost::lexical_cast< std::string >( propagatorSettings->stateType_ ) );
    }

    return integratedStateList;
}

}

}

#endif // PROPAGATIONSETTINGS_H
