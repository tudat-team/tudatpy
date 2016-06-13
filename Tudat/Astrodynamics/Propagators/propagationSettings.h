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

#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"

namespace tudat
{

namespace propagators
{


//! Enum listing types of dynamics that can be numerically integrated
enum IntegratedStateType
{
    hybrid = 0,
    transational_state = 1,
    body_mass_state = 2
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
 *  Base class for defining setting of a propagator. This class is non-functional, and each state
 *  type requires its own derived class (which may have multiple derived classes of its own).
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
                        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates ):
        stateType_( stateType ), initialStates_( initialBodyStates ), stateSize_( initialBodyStates.rows( ) ){ }

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
    virtual void resetInitialStates( const Eigen::Matrix< StateScalarType,
                                     Eigen::Dynamic, 1 >& initialBodyStates )
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

protected:

    //!  Initial state used as input for numerical integration
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    //! Total size of the propagated state.
    int stateSize_;

};

//! Class for defining settings for propagating translational dynamics.
/*!
 *  Class for defining settings for propagating translational dynamics. The propagator defines the
 *  form of the equations of motion (i.e. Cowell, Encke, Gauss etc.). This base class can be used
 *  for Cowell propagator.  Other propagators have dedicated derived class.
 */
template< typename StateScalarType = double >
class TranslationalStatePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor of translational state propagator settings
    /*!
     * Constructor creating translational state propagator settings object
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are
     *  propagated.     
     * \param accelerationsMap A map containing the list of accelerations acting on each body,
     *  identifying the body being acted on and the body acted on by an acceleration. The map has as
     *  key a string denoting the name of the body the list of accelerations, provided as the value
     *  corresponding to a key, is acting on.  This map-value is again a map with string as key,
     *  denoting the body exerting the acceleration, and as value a pointer to an acceleration model.     
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param propagator Type of translational state propagator to be used
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType,
                                                               Eigen::Dynamic, 1 >& initialBodyStates,
                                          const TranslationalPropagatorType propagator = cowell):
        PropagatorSettings< StateScalarType >( transational_state, initialBodyStates ),
        centralBodies_( centralBodies ), accelerationsMap_( accelerationsMap ),
        bodiesToIntegrate_( bodiesToIntegrate ), propagator_( propagator ){ }

    //! Destructor
    ~TranslationalStatePropagatorSettings( ){ }

    //! List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
    std::vector< std::string > centralBodies_;

    //! A map containing the list of accelerations acting on each body
    /*!
     *  A map containing the list of accelerations acting on each body, identifying the body being
     *  acted on and the body acted on by an acceleration. The map has as key a string denoting the
     *  name of the body the list of accelerations, provided as the value corresponding to a key, is
     *  acting on.  This map-value is again a map with string as key, denoting the body exerting the
     *  acceleration, and as value a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationsMap_;

    //! List of bodies for which the translational state is to be propagated.
    std::vector< std::string > bodiesToIntegrate_;

    //! Type of translational state propagator to be used
    TranslationalPropagatorType propagator_;

};


//! Class for defining settings for propagating the mass of a body
/*!
 *  Class for defining settings for propagating the mass of a body. The body masses are propagated in their natural
 *  form (i.e. no choice of equations of motion as is the case for translational dynamics)l
 */
template< typename StateScalarType >
class MassPropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor of mass state propagator settings
    /*!
     * Constructor  of mass state propagator settings
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateModels List of mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses ):
        PropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses ),
      bodiesWithMassToPropagate_( bodiesWithMassToPropagate ), massRateModels_( massRateModels )
    { }

    //! List of bodies for which the mass is to be propagated.
    std::vector< std::string > bodiesWithMassToPropagate_;

    //! List of mass rate models per propagated body.
    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels_;
};

//! Function to retrieve the state size for a list of propagator settings.
/*!
 *  Function to retrieve the initial state for a list of propagator settings.
 *  \param propagatorSettingsList List of propagator settings (sorted by type as key). Map value provides list
 *  of propagator settings for given type.
 *  \return Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the
 *  vector of PropagatorSettings of given type.
 */
template< typename StateScalarType >
int getMultiTypePropagatorStateSize(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    int stateSize = 0;

    // Iterate over all propagation settings and add size to list
    for( typename std::map< IntegratedStateType,
         std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            stateSize += typeIterator->second.at( i )->getStateSize( );
        }
    }
    return stateSize;
}

//! Function to retrieve the initial state for a list of propagator settings.
/*!
 *  Function to retrieve the initial state for a list of propagator settings.
 *  \param propagatorSettingsList List of propagator settings (sorted by type as key). Map value provides list
 *  of propagator settings for given type.
 *  \return Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the
 *  vector of PropagatorSettings of given type.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > createCombinedInitialState(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    // Get total size of propagated state
    int totalSize = getMultiTypePropagatorStateSize( propagatorSettingsList );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > combinedInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalSize, 1 );

    // Iterate over all propagation settings and add to total list
    int currentIndex = 0;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState;
    for( typename std::map< IntegratedStateType,
         std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
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

//! Class for defining settings for propagating multiple types of dynamics concurrently.
template< typename StateScalarType >
class MultiTypePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsMap List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     */
    MultiTypePropagatorSettings(
            const std::map< IntegratedStateType,
            std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorSettingsMap ):
        PropagatorSettings< StateScalarType >(
            hybrid, createCombinedInitialState< StateScalarType >( propagatorSettingsMap ) ),
        propagatorSettingsMap_( propagatorSettingsMap )
    { }

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsVector Vector of propagator settings to use.
     */
    MultiTypePropagatorSettings(
            const std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > propagatorSettingsVector ):
        PropagatorSettings< StateScalarType >(
            hybrid, Eigen::VectorXd::Zero( 0 ) )
    {
        for( unsigned int i = 0; i < propagatorSettingsVector.size( ); i++ )
        {
            propagatorSettingsMap_[ propagatorSettingsVector.at( i )->stateType_ ].push_back(
                        propagatorSettingsVector.at( i ) );
        }

        this->initialStates_ = createCombinedInitialState< StateScalarType >( propagatorSettingsMap_ );
    }

    //! Destructor
    ~MultiTypePropagatorSettings( ){ }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration, sorted in order of
     * IntegratedStateType, and then in the order of the vector of PropagatorSettings of given type.
     */
    void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        // Iterate over all propagator settings.
        int currentStartIndex = 0;
        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = propagatorSettingsMap_.begin( ); propagatorIterator != propagatorSettingsMap_.end( );
             propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                // Get current state size
                int currentParameterSize = propagatorIterator->second.at( i )->getInitialStates( ).rows( );

                // Check consistency
                if( currentParameterSize + currentStartIndex > initialBodyStates.rows( ) )
                {
                    throw std::runtime_error(
                                "Error when resetting multi-type state, sizes are incompatible " );
                }

                // Reset state for current settings
                propagatorIterator->second.at( i )->resetInitialStates(
                            initialBodyStates.block( currentStartIndex, 0, currentParameterSize, 1 ) );
                currentStartIndex += currentParameterSize;

            }
        }

        // Check consistency
        if( currentStartIndex != initialBodyStates.rows( ) )
        {
            std::string errorMessage = "Error when resetting multi-type state, total size is incompatible "+
                    boost::lexical_cast< std::string >( currentStartIndex ) +
                    boost::lexical_cast< std::string >( initialBodyStates.rows( ) );
            throw std::runtime_error( errorMessage );
        }
    }

    //! List of propagator settings to use
    /*!
     * List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     */

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >
    propagatorSettingsMap_;

};

template< typename StateScalarType >
//! Function to retrieve the list of integrated state types and reference ids
/*!
 * Function to retrieve the list of integrated state types and reference ids. For translational and
 * rotational dynamics, the id refers only to the body being propagated (and the second entry of the
 * pair is empty: ""). For proper time propagation, a body and a reference point may be provided,
 * resulting in non-empty first and second pair entries.
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
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > singleTypeIntegratedStateList;

        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            if( typeIterator->first != hybrid )
            {
                for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
                {
                    singleTypeIntegratedStateList = getIntegratedTypeAndBodyList< StateScalarType >(
                                typeIterator->second.at( i ) );

                    if( singleTypeIntegratedStateList.begin( )->first != typeIterator->first
                            || singleTypeIntegratedStateList.size( ) != 1 )
                    {
                        std::cerr<<"Error when making integrated state list for hybrid propagator, inconsistency encountered "<<
                                   singleTypeIntegratedStateList.begin( )->first<<" "<<typeIterator->first<<" "<<
                                   singleTypeIntegratedStateList.size( )<<" "<<
                                   singleTypeIntegratedStateList.begin( )->second.size( )<<std::endl;
                    }
                    else
                    {
                        for( unsigned int j = 0; j < singleTypeIntegratedStateList[ typeIterator->first ].size( ); j++ )
                        {
                            integratedStateList[ typeIterator->first ].push_back(
                                        singleTypeIntegratedStateList.begin( )->second.at( j ) );
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
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationalPropagatorSettings = boost::dynamic_pointer_cast<
                     TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationalPropagatorSettings == NULL )
        {
            throw std::runtime_error( "Error getting integrated state type list, translational state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < translationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair(
                translationalPropagatorSettings->bodiesToIntegrate_.at( i ), "" ) );
        }
        integratedStateList[ transational_state ] = integratedBodies;

        break;
    }
    case body_mass_state:
    {
        boost::shared_ptr< MassPropagatorSettings< StateScalarType > >
                massPropagatorSettings = boost::dynamic_pointer_cast<
                     MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == NULL )
        {
            throw std::runtime_error( "Error getting integrated state type list, mass state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < massPropagatorSettings->bodiesWithMassToPropagate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair(
                massPropagatorSettings->bodiesWithMassToPropagate_.at( i ), "" ) );
        }
        integratedStateList[ body_mass_state ] = integratedBodies;

        break;
    }
    default:
        throw std::runtime_error( "Error, could not process integrated state type ingetIntegratedTypeAndBodyList " +
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

}

#endif // TUDAT_PROPAGATIONSETTINGS_H
