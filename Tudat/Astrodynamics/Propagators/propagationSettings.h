/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONSETTINGS_H
#define TUDAT_PROPAGATIONSETTINGS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <unordered_map>

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{

namespace propagators
{


//! Enum listing types of dynamics that can be numerically integrated
enum IntegratedStateType
{
    transational_state
};


//! Enum listing propagator types for translational dynamics that can be used.
enum TranslationalPropagatorType
{
    cowell = 0
};

enum PropagationDependentVariables
{
    mach_number_dependent_variable,
    altitude_dependent_variable,
    airspeed_dependent_variable,
    local_density_dependent_variable,
    relative_speed_dependent_variable,
    relative_distance_dependent_variable,
    radiation_pressure_dependent_variable,    
    total_acceleration_norm_dependent_variable
};

enum PropagationTerminationTypes
{
    time_stopping_condition,
    dependent_variable_stopping_condition,
    hybrid_stopping_condition
};


class PropagationTerminationSettings
{
public:
    PropagationTerminationSettings( const PropagationTerminationTypes terminationType ):
        terminationType_( terminationType ){ }

    virtual ~PropagationTerminationSettings( ){ }

    PropagationTerminationTypes terminationType_;
};

class PropagationTimeTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationTimeTerminationSettings( const double terminationTime ):
        PropagationTerminationSettings( time_stopping_condition ),
        terminationTime_( terminationTime ){ }

    ~PropagationTimeTerminationSettings( ){ }

    double terminationTime_;
};

class PropagationDependentVariableTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationDependentVariableTerminationSettings( const PropagationDependentVariables variableType,
                                                     const std::string associatedBody,
                                                     const double limitValue,
                                                     const bool useAsLowerLimit,
                                                     const std::string secondaryBody = "" ):
        PropagationTerminationSettings( dependent_variable_stopping_condition ),
        variableType_( variableType ), associatedBody_( associatedBody ),
        limitValue_( limitValue ), useAsLowerLimit_( useAsLowerLimit ), secondaryBody_( secondaryBody ){ }

    ~PropagationDependentVariableTerminationSettings( ){ }

    PropagationDependentVariables variableType_;

    std::string associatedBody_;

    double limitValue_;

    bool useAsLowerLimit_;

    std::string secondaryBody_;
};

class PropagationHybridTerminationSettings: public PropagationTerminationSettings
{
public:
    PropagationHybridTerminationSettings(
            const std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings,
            const bool fulFillSingleCondition = 0 ):
        PropagationTerminationSettings( hybrid_stopping_condition ),
        terminationSettings_( terminationSettings ),
        fulFillSingleCondition_( fulFillSingleCondition ){ }

    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettings_;

    bool fulFillSingleCondition_;

};

//! Get size of state for single propagated state of given type.
/*!
 * Get size of state for single propagated state of given type (i.e. 6 for translational state).
 * \param stateType Type of state
 * \return Size of single state.
 */
int getSingleIntegrationSize( const IntegratedStateType stateType );

//! Get order of differential equation for governing equations of dynamics of given type.
/*!
 * Get order of differential equation for governing equations of dynamics of given type (i.e. 2 for translational state).
 * \param stateType Type of state
 * \return Order of differential equations.
 */
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType );

//! Base class for defining setting of a propagator
/*!
 *  Base class for defining setting of a propagator. This class is non-functional, and each state type requires its
 *  own derived class (which may have multiple derived classes of its own).
 */
template< typename StateScalarType >
class PropagatorSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state being propagated
     * \param initialBodyStates Initial state used as input for numerical integration
     */
    PropagatorSettings( const IntegratedStateType stateType,
                        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
                        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings ):
        stateType_( stateType ), initialStates_( initialBodyStates ), stateSize_( initialBodyStates.rows( ) ),
    terminationSettings_( terminationSettings ){ }

    //! Virtual destructor.
    virtual ~PropagatorSettings( ){ }

    //!T ype of state being propagated
    IntegratedStateType stateType_;

    //! Function to retrieve the initial state used as input for numerical integration
    /*!
     * Function to retrieve the initial state used as input for numerical integration
     * \return Initial state used as input for numerical integration
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( )
    {
        return initialStates_;
    }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration
     */
    virtual void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        initialStates_ = initialBodyStates;
        stateSize_ = initialStates_.rows( );
    }

    //! Get total size of the propagated state.
    /*!
     * Get total size of the propagated state.
     * \return Total size of the propagated state.
     */
    int getStateSize( )
    {
        return stateSize_;
    }

    boost::shared_ptr< PropagationTerminationSettings > getTerminationSettings( )
    {
        return terminationSettings_;
    }

protected:

    //!  Initial state used as input for numerical integration
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    //! Total size of the propagated state.
    int stateSize_;

    boost::shared_ptr< PropagationTerminationSettings > terminationSettings_;

};

//! Class for defining settings for propagating translational dynamics.
/*!
 *  Class for defining settings for propagating translational dynamics. The propagator defines the form of the equations of
 *  motion (i.e. Cowell, Encke, Gauss etc.). This base class can be used for Cowell propagator.
 *  Other propagators have dedicated derived class.
 */
template< typename StateScalarType = double >
class TranslationalStatePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param propagator Type of translational state propagator to be used
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                          const TranslationalPropagatorType propagator = cowell):
        PropagatorSettings< StateScalarType >( transational_state, initialBodyStates, terminationSettings ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }

    //! Constructor
    /*!
     * Constructor
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param propagator Type of translational state propagator to be used
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const double endTime,
                                          const TranslationalPropagatorType propagator = cowell):
        PropagatorSettings< StateScalarType >(
            transational_state, initialBodyStates,  boost::make_shared< PropagationTimeTerminationSettings >( endTime ) ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }


    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    ~TranslationalStatePropagatorSettings( ){ }

    //! List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
    std::vector< std::string > centralBodies_;

    //! A map containing the list of accelerations acting on each body
    /*!
     *  A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationsMap_;

    //! List of bodies for which the translational state is to be propagated.
    std::vector< std::string > bodiesToIntegrate_;

    //! Type of translational state propagator to be used
    TranslationalPropagatorType propagator_;

};

template< typename StateScalarType >
//! Function to retrieve the list of integrated state types and reference ids
/*!
* Function to retrieve the list of integrated state types and reference ids. For translational and rotational dynamics,
* the id refers only to the body being propagated (and the second entry of the pair is empty: ""). For proper time
* propagation, a body and a reference point may be provided, resulting in non-empty first and second pair entries.
* \param propagatorSettings Settings that are to be used for the propagation.
* \return List of integrated state types and reference ids
*/
std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStateList;

    // Identify propagator type
    switch( propagatorSettings->stateType_ )
    {    
    case transational_state:
    {

        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationalPropagatorSettings == NULL )
        {
            throw std::runtime_error( "Error getting integrated state type list, translational state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
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

namespace std
{

//! Hash for IntegratedStateType enum.
template< >
struct hash< tudat::propagators::IntegratedStateType >
{
   typedef tudat::propagators::IntegratedStateType argument_type;
   typedef size_t result_type;

   result_type operator () (const argument_type& x) const
   {
      using type = typename std::underlying_type<argument_type>::type;
      return std::hash< type >( )( static_cast< type >( x ) );
   }
};

}

#endif // TUDAT_PROPAGATIONSETTINGS_H
