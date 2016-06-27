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

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Propagators/propagationOutputSettings.h"
#include "Tudat/Astrodynamics/Propagators/propagationTerminationSettings.h"

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
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    PropagatorSettings( const IntegratedStateType stateType,
                        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
                        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
                        const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                        const double printInterval = TUDAT_NAN ):
        stateType_( stateType ), initialStates_( initialBodyStates ), stateSize_( initialBodyStates.rows( ) ),
        terminationSettings_( terminationSettings ), dependentVariablesToSave_( dependentVariablesToSave ),
        printInterval_( printInterval){ }

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

    //! Function to retrieve settings for creating the object that checks whether the propagation is finished.
    /*!
     * Function to retrieve settings for creating the object that checks whether the propagation is finished.
     * \return Settings for creating the object that checks whether the propagation is finished.
     */
    boost::shared_ptr< PropagationTerminationSettings > getTerminationSettings( )
    {
        return terminationSettings_;
    }

    //! Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
    /*!
     * Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
     * \return Settings for the dependent variables that are to be saved during propagation (default none).
     */
    boost::shared_ptr< DependentVariableSaveSettings > getDependentVariablesToSave( )
    {
        return dependentVariablesToSave_;
    }

    //! Function to retrieve how often the current state and time are to be printed to console
    /*!
     * Function to retrieve how often the current state and time are to be printed to console
     * \return Time intercal with which the current state and time are to be printed to console (default NaN, meaning never).
     */
    double getPrintInterval( )
    {
        return printInterval_;
    }


protected:

    //!  Initial state used as input for numerical integration
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    //! Total size of the propagated state.
    int stateSize_;

    //! Settings for creating the object that checks whether the propagation is finished.
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings_;

    //! Settings for the dependent variables that are to be saved during propagation (default none).
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;

    //! Variable indicating how often (once per printInterval_ seconds or propagation independenty variable) the
    //! current state and time are to be printed to console (default never).
    double printInterval_;

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

    //! Constructor for generic stopping conditions.
    /*!
     * Constructor for generic stopping conditions.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param propagator Type of translational state propagator to be used
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( transational_state, initialBodyStates, terminationSettings,
                                               dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }

    //! Constructor for fixed propagation time stopping conditions.
    /*!
     * Constructor for fixed propagation time stopping conditions.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param endTime Time at which to stop the numerical propagation
     * \param propagator Type of translational state propagator to be used
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const double endTime,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >(
            transational_state, initialBodyStates,  boost::make_shared< PropagationTimeTerminationSettings >( endTime ),
            dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }


    //! Destructor
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


} // namespace propagators

} // namespace tudat

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

} // namespace std

#endif // TUDAT_PROPAGATIONSETTINGS_H
