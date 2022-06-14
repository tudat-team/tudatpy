/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"

namespace tudat
{

namespace propagators
{

//! Base class for defining propagation settings, derived classes split into settings for single- and multi-arc dynamics
template< typename StateScalarType = double >
class PropagatorSettings
{

public:
    //! Constructor
    /*!
     * Constructor
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param isMultiArc Boolean denoting whether the propagation settings are multi-arc (if true) or single arc (if false).
     */
    PropagatorSettings( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
                        const bool isMultiArc ):
        initialStates_( initialBodyStates ), stateSize_( initialBodyStates.rows( ) ), isMultiArc_( isMultiArc ){ }

    //! Destructor
    virtual ~PropagatorSettings( ){ }

    //! Function to retrieve the initial state used as input for numerical integration
    /*!
     * Function to retrieve the initial state used as input for numerical integration
     * \return Initial state used as input for numerical integration
     */
    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( )
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
    virtual int getPropagatedStateSize( )
    {
        return stateSize_;
    }

    //! Get total size of the conventional state.
    /*!
     * Get total size of the conventional state.
     * \return Total size of the conventional state.
     */
    int getConventionalStateSize( )
    {
        return stateSize_;
    }

    //! Function to get boolean denoting whether the propagation settings are multi-arc (if true) or single arc (if false).
    /*!
     *  Function to get boolean denoting whether the propagation settings are multi-arc (if true) or single arc (if false).
     *  \return Boolean denoting whether the propagation settings are multi-arc (if true) or single arc (if false).
     */
    bool getIsMultiArc( )
    {
        return isMultiArc_;
    }

    //! Function to create the integrated state models (e.g. acceleration/torque/mass-rate models).
    /*!
     * Function to create the integrated state models (e.g. acceleration/torque/mass-rate models).
     *
     * Derived classes must provide an implementation for this method. This function will use the provided
     * system of bodies and the member containing the acceleration/torque/mass-rate settings to create the
     * actual models and assign them to the corresponding model map members.
     *
     * The implementation for MultiArc and MultiType (hybrid state) propagations will call the
     * `resetIntegratedStateModels` method of each of the fundamental propagators that they contain.
     *
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies ) = 0;

protected:

    //!  Initial state used as input for numerical integration
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    //! Total size of the propagated state.
    int stateSize_;

    //! Boolean denoting whether the propagation settings are multi-arc (if true) or single arc (if false).
    bool isMultiArc_;
};


//! Base class for defining setting of a propagator for single-arc dynamics
/*!
 *  Base class for defining setting of a propagator for single-arc dynamics. This class is non-functional, and each state type
 *  requires its own derived class (which may have multiple derived classes of its own).
 */
template< typename StateScalarType = double >
class SingleArcPropagatorSettings: public PropagatorSettings< StateScalarType >
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
    SingleArcPropagatorSettings( const IntegratedStateType stateType,
                                 const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
                                 const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                 const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                 const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( initialBodyStates, false ),
        stateType_( stateType ), terminationSettings_( terminationSettings ),
        dependentVariablesToSave_( dependentVariablesToSave ), printInterval_( printInterval)
    { }

    //! Virtual destructor.
    virtual ~SingleArcPropagatorSettings( ){ }

    //! Type of state being propagated.
    IntegratedStateType getStateType( )
    {
        return stateType_;
    }

    //! Function to retrieve settings for creating the object that checks whether the propagation is finished.
    /*!
     * Function to retrieve settings for creating the object that checks whether the propagation is finished.
     * \return Settings for creating the object that checks whether the propagation is finished.
     */
    std::shared_ptr< PropagationTerminationSettings > getTerminationSettings( )
    {
        return terminationSettings_;
    }

    //! Function to reset settings for creating the object that checks whether the propagation is finished.
    /*!
     * Function to reset settings for creating the object that checks whether the propagation is finished.
     * \param terminationSettings New settings for creating the object that checks whether the propagation is finished.
     */
    void resetTerminationSettings( const std::shared_ptr< PropagationTerminationSettings > terminationSettings )
    {
        terminationSettings_ = terminationSettings;
    }

    //! Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
    /*!
     * Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
     * \return Settings for the dependent variables that are to be saved during propagation (default none).
     */
    std::shared_ptr< DependentVariableSaveSettings > getDependentVariablesToSave( )
    {
        return dependentVariablesToSave_;
    }

    //! Function to reset settings for the dependent variables that are to be saved during propagation
    /*!
     * Function to reset settings for the dependent variables that are to be saved during propagation.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation.
     */
    void resetDependentVariablesToSave(
            const std::shared_ptr< DependentVariableSaveSettings >& dependentVariablesToSave )
    {
        dependentVariablesToSave_ = dependentVariablesToSave;
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

    //! Function to modify settings for creating the object that checks whether the propagation is finished.
    /*!
     * Function to modify settings for creating the object that checks whether the propagation is finished.
     * \param terminationSettings new settings for creating the object that checks whether the propagation is finished.
     */
    void setTerminationSettings(
            std::shared_ptr< PropagationTerminationSettings > terminationSettings )
    {
        terminationSettings_ = terminationSettings;
    }

protected:

    //!Type of state being propagated
    IntegratedStateType stateType_;

    //! Settings for creating the object that checks whether the propagation is finished.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings_;

    //! Settings for the dependent variables that are to be saved during propagation (default none).
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;

    //! Variable indicating how often (once per printInterval_ seconds or propagation independenty variable) the
    //! current state and time are to be printed to console (default never).
    double printInterval_;

};

//! Function to get the total size of multi-arc initial state vector
/*!
 *  Function to get the total size of multi-arc initial state vector, e.g. the size of the single-arc initial states, concatenated
 *  into a single vector
 *  \param singleArcPropagatorSettings ist of single-arc propagation settings for which the concatenated initial state size is to
 *  be determined.
 *  \return Total size of multi-arc initial state vector
 */
template< typename StateScalarType = double >
int getConcatenatedStateSize(
        const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >& singleArcPropagatorSettings )
{
    int vectorSize = 0;

    for( unsigned int i = 0; i < singleArcPropagatorSettings.size( ); i++ )
    {
        vectorSize += singleArcPropagatorSettings.at( i )->getConventionalStateSize( );
    }

    return vectorSize;
}

//! Function to concatenate the initial states for a list of single-arc propagations into a single list.
/*!
 *  Function to concatenate the initial states for a list of single-arc propagations into a single list.
 *  \param singleArcPropagatorSettings List of single-arc propagation settings for which the initial states are to be
 *  concatenated into a single vector.
 *  \return Vector with concatenated initial states from singleArcPropagatorSettings.
 */
template< typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getConcatenatedInitialStates(
        const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >& singleArcPropagatorSettings )
{
    // Define size of return vector
    int vectorSize = getConcatenatedStateSize( singleArcPropagatorSettings );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( vectorSize );

    // Retrieve single-arc initial states arc-by-arac
    int currentIndex = 0;
    int currentBlockSize = 0;
    for( unsigned int i = 0; i < singleArcPropagatorSettings.size( ); i++ )
    {
        currentBlockSize = singleArcPropagatorSettings.at( i )->getConventionalStateSize( );
        initialStates.segment( currentIndex, currentBlockSize ) = singleArcPropagatorSettings.at( i )->getInitialStates( );
        currentIndex += currentBlockSize;
    }

    return initialStates;
}

//! Class for defining setting of a propagator for multi-arc dynamics
/*!
 *  Base class for defining setting of a propagator for multi-arc dynamics. This class contains single-arc propagator settings
 *  for each arc.
 */
template< typename StateScalarType = double >
class MultiArcPropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //!  Constructor
    /*!
     * Constructor
     * \param singleArcSettings List of propagator settings for each arc in propagation.
     * \param transferInitialStateInformationPerArc Boolean denoting whether the initial state of arc N+1 is to be taken from
     * arc N (for N>0)
     */
    MultiArcPropagatorSettings(
            const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >& singleArcSettings,
            const bool transferInitialStateInformationPerArc = 0 ):
        PropagatorSettings< StateScalarType >( getConcatenatedInitialStates( singleArcSettings ), true )
    {
        singleArcSettings_ = singleArcSettings;
        for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
        {
            // If information is to be transferred between arcs, set arc N>0 initial state as NaN
            if( transferInitialStateInformationPerArc )
            {
                if( i != 0 )
                {
                    singleArcSettings_.at( i )->resetInitialStates(
                                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Constant(
                                    singleArcSettings_.at( i )->getConventionalStateSize( ), TUDAT_NAN ) );
                }
            }

            initialStateList_.push_back( singleArcSettings_.at( i )->getInitialStates( ) );
        }
        if( transferInitialStateInformationPerArc )
        {
            this->initialStates_ = getConcatenatedInitialStates( singleArcSettings );
        }
    }

    //! Destructor
    virtual ~MultiArcPropagatorSettings( ){ }

    //! Function get the list of propagator settings for each arc in propagation.
    /*!
     * Function get the list of propagator settings for each arc in propagation.
     * \return List of propagator settings for each arc in propagation.
     */
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > getSingleArcSettings( )
    {
        return singleArcSettings_;
    }

    //! Function to retrieve the number of arcs
    /*!
     * Function to retrieve the number of arcs
     * \return Number of arcs
     */
    int getNmberOfArcs( )
    {
        return singleArcSettings_.size( );
    }

    //! Function get the list of initial states for each arc in propagation.
    /*!
     * Function get the list of initial states for each arc in propagation.
     * \return List of initial states for each arc in propagation.
     */
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getInitialStateList( )
    {
        return initialStateList_;
    }

    //! Function to create the integrated state models (e.g. acceleration/torque/mass-rate models).
    /*!
     * Function to create the integrated state models (e.g. acceleration/torque/mass-rate models) for
     * each fo the propagators existing in each propagation arc.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        for ( std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcSettings :
              singleArcSettings_ )
        {
            if ( singleArcSettings )
            {
                singleArcSettings->resetIntegratedStateModels( bodies );
            }
        }
    }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration
     */
    void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        if( this->stateSize_ != this->initialStates_.rows( ) )
        {
            std::cerr << "Warning when resetting multi-arc initial states, size is incomparible with original size." << std::endl;
        }

        this->initialStates_ = initialBodyStates;

        int currentIndex = 0;
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            initialStateList_[ i ] = this->initialStates_.segment( currentIndex, singleArcSettings_.at( i )->getConventionalStateSize( ) );
            singleArcSettings_.at( i )->resetInitialStates( initialStateList_[ i ] );
            currentIndex += singleArcSettings_.at( i )->getConventionalStateSize( );
        }

    }

    //! Function to reset the initial state used as input for numerical integration as a vector of Eigen Vectors
    /*!
     * Function to reset the initial state used as input for numerical integration as a vector of Eigen Vectors
     * \param initialStateList New initial states used as input for numerical integration
     */
    void resetInitialStatesList( const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStateList )
    {
        if( initialStateList_.size( ) != initialStateList.size( ) )
        {
            std::cerr << "Warning when resetting multi-arc initial state list, size is incomparible with original size." << std::endl;
        }

        initialStateList_ = initialStateList;

        int currentIndex = 0;
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            this->initialStates_.segment( currentIndex, singleArcSettings_.at( i )->getConventionalStateSize( ) ) =  initialStateList_[ i ];
            singleArcSettings_.at( i )->resetInitialStates( initialStateList_[ i ] );
            currentIndex += singleArcSettings_.at( i )->getConventionalStateSize( );
        }
    }

    void updateInitialStateFromConsituentSettings( )
    {
        int currentIndex = 0;
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            this->initialStates_.segment( currentIndex, singleArcSettings_.at( i )->getConventionalStateSize( ) ) =
                    singleArcSettings_.at( i )->getInitialStates( );
            currentIndex += singleArcSettings_.at( i )->getConventionalStateSize( );
        }
    }

protected:

    //! List of propagator settings for each arc in propagation.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > singleArcSettings_;

    //! List of initial states for each arc in propagation.
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > initialStateList_;
};

template< typename StateScalarType = double >
std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings(
        const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >& singleArcSettings,
        const bool transferInitialStateInformationPerArc = 0 )
{
        return std::make_shared< MultiArcPropagatorSettings< StateScalarType > >(
                    singleArcSettings, transferInitialStateInformationPerArc );
}

//! Class for defining setting of a propagator for a combination of single- and multi-arc dynamics
template< typename StateScalarType = double >
class HybridArcPropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param singleArcPropagatorSettings Settings for single-arc propagation component
     * \param multiArcPropagatorSettings Settings for multi-arc propagation component
     */
    HybridArcPropagatorSettings(
            const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcPropagatorSettings,
            const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings ):
        PropagatorSettings< StateScalarType >(
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 0 ), false ),
        singleArcPropagatorSettings_( singleArcPropagatorSettings ),
        multiArcPropagatorSettings_( multiArcPropagatorSettings )
    {
        // Set initial states
        this->initialStates_ =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >(
                    singleArcPropagatorSettings->getPropagatedStateSize( ) +
                    multiArcPropagatorSettings->getPropagatedStateSize( ) );
        setInitialStatesFromConstituents( );

        // Set initial state sizes
        this->stateSize_ = this->initialStates_.rows( );
        singleArcStateSize_ = singleArcPropagatorSettings_->getPropagatedStateSize( );
        multiArcStateSize_ = multiArcPropagatorSettings_->getPropagatedStateSize( );
    }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration
     */
    void resetInitialStates(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        this->initialStates_ = initialBodyStates;
        this->stateSize_ = this->initialStates_.rows( );

        singleArcPropagatorSettings_->resetInitialStates(
                    this->initialStates_.segment( 0, singleArcStateSize_ ) );
        multiArcPropagatorSettings_->resetInitialStates(
                    this->initialStates_.segment( singleArcStateSize_, multiArcStateSize_ ) );
    }

    //! Function to retrieve settings for single-arc propagation component
    /*!
     * Function to retrieve settings for single-arc propagation component
     * \return Settings for single-arc propagation component
     */
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > getSingleArcPropagatorSettings( )
    {
        return singleArcPropagatorSettings_;
    }

    //! Function to retrieve settings for multi-arc propagation component
    /*!
     * Function to retrieve settings for multi-arc propagation component
     * \return Settings for multi-arc propagation component
     */
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > getMultiArcPropagatorSettings( )
    {
        return multiArcPropagatorSettings_;
    }

    //! Function that sets initial states from single- and multi-arc initial states
    void setInitialStatesFromConstituents( )
    {
        this->initialStates_.segment( 0, singleArcPropagatorSettings_->getPropagatedStateSize( ) ) =
                singleArcPropagatorSettings_->getInitialStates( );
        this->initialStates_.segment(
                    singleArcPropagatorSettings_->getPropagatedStateSize( ), multiArcPropagatorSettings_->getPropagatedStateSize( ) ) =
                multiArcPropagatorSettings_->getInitialStates( );
    }

    //! Function to create the integrated state models (e.g. acceleration/torque/mass-rate models).
    /*!
     * Function to create the integrated state models (e.g. acceleration/torque/mass-rate models) from associated settings objects.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        singleArcPropagatorSettings_->resetIntegratedStateModels( bodies );
        multiArcPropagatorSettings_->resetIntegratedStateModels( bodies );
    }

protected:

    //! Settings for single-arc propagation component
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcPropagatorSettings_;

    //! Settings for multi-arc propagation component
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings_;

    //! Size of total single-arc initial state
    int singleArcStateSize_;

    //! Size of total multi-arc initial state
    int multiArcStateSize_;
};


template< typename StateScalarType = double >
std::shared_ptr< HybridArcPropagatorSettings< StateScalarType > > hybridArcPropagatorSettings(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcPropagatorSettings,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings )
{
        return std::make_shared< HybridArcPropagatorSettings< StateScalarType > >(
                    singleArcPropagatorSettings, multiArcPropagatorSettings );
}

//! Class for defining settings for propagating translational dynamics.
/*!
 *  Class for defining settings for propagating translational dynamics. The propagator defines the form of the equations of
 *  motion (i.e. Cowell, Encke, Gauss etc.). This base class can be used for Cowell propagator.
 *  Other propagators have dedicated derived class.
 */
template< typename StateScalarType = double >
class TranslationalStatePropagatorSettings: public SingleArcPropagatorSettings< StateScalarType >
{
public:

    //! Constructor for generic stopping conditions, providing an alreay-created accelerations map.
    /*!
     * Constructor for generic stopping conditions, providing an alreay-created accelerations map.
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
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independent
     * variable) the current state and time are to be printed to console (default never).
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( translational_state, initialBodyStates, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ),
        accelerationsMap_( accelerationsMap ) { verifyInput( ); }

    //! Constructor for generic stopping conditions, providing settings to create accelerations map.
    /*!
     * Constructor for generic stopping conditions, providing settings to create accelerations map.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationSettingsMap A map containing settings for the accelerations acting on each body, identifying
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
                                          const simulation_setup::SelectedAccelerationMap& accelerationSettingsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( translational_state, initialBodyStates, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ),
        accelerationSettingsMap_( accelerationSettingsMap ) { verifyInput( ); }

    //! Constructor for fixed propagation time stopping conditions, providing an alreay-created accelerations map.
    /*!
     * Constructor for fixed propagation time stopping conditions, providing an alreay-created accelerations map.
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
                                          const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >(
            translational_state, initialBodyStates,  std::make_shared< PropagationTimeTerminationSettings >( endTime ),
            dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ),
        accelerationsMap_( accelerationsMap ) { verifyInput( ); }

    //! Constructor for fixed propagation time stopping conditions, providing settings to create accelerations map.
    /*!
     * Constructor for fixed propagation time stopping conditions, providing settings to create accelerations map.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationSettingsMap A map containing settings for the accelerations acting on each body, identifying
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
                                          const simulation_setup::SelectedAccelerationMap& accelerationSettingsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const double endTime,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >(
            translational_state, initialBodyStates,  std::make_shared< PropagationTimeTerminationSettings >( endTime ),
            dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ),
        accelerationSettingsMap_( accelerationSettingsMap ) { verifyInput( ); }

    //! Destructor
    ~TranslationalStatePropagatorSettings( ){ }

    //! List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
    std::vector< std::string > centralBodies_;

    //! List of bodies for which the translational state is to be propagated.
    std::vector< std::string > bodiesToIntegrate_;

    //! Type of translational state propagator to be used
    TranslationalPropagatorType propagator_;

    //! Function to create the acceleration models.
    /*!
     * Function to create the acceleration models.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        accelerationsMap_ = simulation_setup::createAccelerationModelsMap(
                    bodies, accelerationSettingsMap_, bodiesToIntegrate_, centralBodies_ );
    }

    void resetAccelerationModelsMap(
            const simulation_setup::SelectedAccelerationMap accelerationSettingsMap,
            const simulation_setup::SystemOfBodies& bodies )
    {
        accelerationSettingsMap_ = accelerationSettingsMap;
        accelerationsMap_ = simulation_setup::createAccelerationModelsMap(
                    bodies, accelerationSettingsMap_, bodiesToIntegrate_, centralBodies_ );
    }

    //! Function to get the acceleration settings map.
    /*!
     * Function to get the acceleration settings map.
     * \return The acceleration settings map.
     */
    simulation_setup::SelectedAccelerationMap getAccelerationSettingsMap( ) const
    {
        return accelerationSettingsMap_;
    }

    //! Function to get the accelerations map.
    /*!
     * Function to get the accelerations map.
     * \return The accelerations map.
     */
    basic_astrodynamics::AccelerationMap getAccelerationsMap( ) const
    {
        if ( accelerationsMap_.size( ) == 0 && accelerationSettingsMap_.size( ) != 0 )
        {
            std::cerr << "Inconsistent sizes for map of aceleration settings and map of acceleration models. "
                      << "Did you forget to call resetIntegratedStateModels on the propagator?" << std::endl;
        }
        return accelerationsMap_;
    }

    int getPropagatedStateSize( )
    {
        return ( propagator_ == unified_state_model_quaternions ||
                 propagator_ == unified_state_model_modified_rodrigues_parameters ) ?
                    ( this->stateSize_ + bodiesToIntegrate_.size( ) ) : this->stateSize_;
    }

private:

    void verifyInput( )
    {
        if( this->initialStates_.rows( ) != static_cast< int >( 6 * bodiesToIntegrate_.size( ) ) )
        {
            throw std::runtime_error( "Error when defining body translational propagator settings, provided initial state size (" +
                                      std::to_string( this->initialStates_.rows( ) ) +
                                      ") is incompatible with list of bodies for which translational state is to be propagated (size " +
                                      std::to_string( bodiesToIntegrate_.size( ) ) +
                                      ")");
        }
    }

    //! A map containing the list of settings for the accelerations acting on each body
    /*!
     *  A map containing the list of settings for the accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration settings object.
     */
    simulation_setup::SelectedAccelerationMap accelerationSettingsMap_;

    //! A map containing the list of accelerations acting on each body
    /*!
     *  A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationsMap_;

};

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const basic_astrodynamics::AccelerationMap& accelerationsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const TranslationalPropagatorType propagator = cowell,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
        std::shared_ptr< DependentVariableSaveSettings >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationsMap, bodiesToIntegrate, initialBodyStates,
                terminationSettings, propagator, dependentVariablesToSave, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const basic_astrodynamics::AccelerationMap& accelerationsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const TranslationalPropagatorType propagator = cowell,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariablesToSave =
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesTotal = nullptr;
    if( dependentVariablesToSave.size( ) > 0 )
    {
        dependentVariablesTotal =
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave );
    }
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationsMap, bodiesToIntegrate, initialBodyStates,
                terminationSettings, propagator, dependentVariablesTotal, printInterval );
}


template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SelectedAccelerationMap& accelerationSettingsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const TranslationalPropagatorType propagator = cowell,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
        std::shared_ptr< DependentVariableSaveSettings >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationSettingsMap, bodiesToIntegrate, initialBodyStates,
                terminationSettings, propagator, dependentVariablesToSave, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SelectedAccelerationMap& accelerationSettingsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const TranslationalPropagatorType propagator = cowell,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariablesToSave =
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesTotal = nullptr;
    if( dependentVariablesToSave.size( ) > 0 )
    {
        dependentVariablesTotal =
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave );
    }
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationSettingsMap, bodiesToIntegrate, initialBodyStates,
                terminationSettings, propagator, dependentVariablesTotal, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const basic_astrodynamics::AccelerationMap& accelerationsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const double endTime,
        const TranslationalPropagatorType propagator = cowell,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
        std::shared_ptr< DependentVariableSaveSettings >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationsMap, bodiesToIntegrate, initialBodyStates,
                endTime, propagator, dependentVariablesToSave, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const basic_astrodynamics::AccelerationMap& accelerationsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const double endTime,
        const TranslationalPropagatorType propagator = cowell,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariablesToSave =
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesTotal = nullptr;
    if( dependentVariablesToSave.size( ) > 0 )
    {
        dependentVariablesTotal =
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave );
    }
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationsMap, bodiesToIntegrate, initialBodyStates,
                endTime, propagator, dependentVariablesTotal, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
translationalStatePropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SelectedAccelerationMap& accelerationSettingsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const double endTime,
        const TranslationalPropagatorType propagator = cowell,
        const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
        std::shared_ptr< DependentVariableSaveSettings >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< TranslationalStatePropagatorSettings < StateScalarType > >(
                centralBodies, accelerationSettingsMap, bodiesToIntegrate, initialBodyStates,
                endTime, propagator, dependentVariablesToSave, printInterval );
}

//! Class for defining settings for propagating rotational dynamics.
template< typename StateScalarType = double >
class RotationalStatePropagatorSettings: public SingleArcPropagatorSettings< StateScalarType >
{
public:

    //! Constructor with already-created torque models.
    /*!
     * Constructor with already-created torque models.
     * \param torqueModelMap List of torque models that are to be used in propagation
     * \param bodiesToIntegrate List of bodies that are to be propagated numerically.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param propagator Type of rotational state propagator to be used.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     *
     */
    RotationalStatePropagatorSettings( const basic_astrodynamics::TorqueModelMap& torqueModelMap,
                                       const std::vector< std::string >& bodiesToIntegrate,
                                       const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                       const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                       const RotationalPropagatorType propagator = quaternions,
                                       const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                       const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( rotational_state, initialBodyStates, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        bodiesToIntegrate_( bodiesToIntegrate ), propagator_( propagator ), torqueModelMap_( torqueModelMap )
    { verifyInput( ); }

    //! Constructor with settings for torque models.
    /*!
     * Constructor with settings for torque models.
     * \param torqueSettingsMap List of settings for torque models that are to be used in propagation
     * \param bodiesToIntegrate List of bodies that are to be propagated numerically.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param propagator Type of rotational state propagator to be used.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     *
     */
    RotationalStatePropagatorSettings( const simulation_setup::SelectedTorqueMap& torqueSettingsMap,
                                       const std::vector< std::string >& bodiesToIntegrate,
                                       const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                       const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                       const RotationalPropagatorType propagator = quaternions,
                                       const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
                                       const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( rotational_state, initialBodyStates, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        bodiesToIntegrate_( bodiesToIntegrate ), propagator_( propagator ), torqueSettingsMap_( torqueSettingsMap )
    { verifyInput( ); }

    //! Destructor
    ~RotationalStatePropagatorSettings( ){ }

    //! List of bodies that are to be propagated numerically.
    std::vector< std::string > bodiesToIntegrate_;

    //! Type of translational state propagator to be used
    RotationalPropagatorType propagator_;

    //! Function to create the torque models.
    /*!
     * Function to create the torque models.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        torqueModelMap_ = simulation_setup::createTorqueModelsMap( bodies, torqueSettingsMap_, bodiesToIntegrate_ );
    }

    //! Function to get the torque settings map.
    /*!
     * Function to get the torque settings map.
     * \return The torque settings map.
     */
    simulation_setup::SelectedTorqueMap getTorqueSettingsMap( ) const
    {
        return torqueSettingsMap_;
    }

    //! Function to get the torque models map.
    /*!
     * Function to get the torque models map.
     * \return The torque models map.
     */
    basic_astrodynamics::TorqueModelMap getTorqueModelsMap( ) const
    {
        if ( torqueModelMap_.size( ) == 0 && torqueSettingsMap_.size( ) != 0 )
        {
            std::cerr << "Inconsistent sizes for map of torque settings and map of torque models. "
                      << "Did you forget to call resetIntegratedStateModels on the propagator?" << std::endl;
        }
        return torqueModelMap_;
    }


    void resetTorqueModelsMap( const basic_astrodynamics::TorqueModelMap& torqueModelMap )
    {
        torqueModelMap_ = torqueModelMap;
    }

private:

    void verifyInput( )
    {
        if( this->initialStates_.rows( ) != static_cast< int >( 7 * bodiesToIntegrate_.size( ) ) )
        {
            throw std::runtime_error( "Error when defining body rotational propagator settings, provided initial state size (" +
                                      std::to_string( this->initialStates_.rows( ) ) +
                                      ") is incompatible with list of bodies for which rotational state is to be propagated (size " +
                                      std::to_string( bodiesToIntegrate_.size( ) ) +
                                      ")");
        }
    }

    //! List of torque settings that are to be used to create the torque models
    simulation_setup::SelectedTorqueMap torqueSettingsMap_;

    //! List of torque models that are to be used in propagation
    basic_astrodynamics::TorqueModelMap torqueModelMap_;

};


template< typename StateScalarType = double >
inline std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModelMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const RotationalPropagatorType propagator = quaternions,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> >(),
        const double printInterval = TUDAT_NAN)
{
    return std::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                torqueModelMap, bodiesToIntegrate, initialBodyStates, terminationSettings, propagator,
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ), printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModelMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const double finalTime,
        const RotationalPropagatorType propagator = quaternions,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                torqueModelMap, bodiesToIntegrate, initialBodyStates,
                std::make_shared< PropagationTimeTerminationSettings >( finalTime ), propagator,
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ), printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings(
        const simulation_setup::SelectedTorqueMap& torqueSettingsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const RotationalPropagatorType propagator = quaternions,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> >(),
        const double printInterval = TUDAT_NAN)
{
    return std::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                torqueSettingsMap, bodiesToIntegrate, initialBodyStates, terminationSettings, propagator,
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ), printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings(
        const simulation_setup::SelectedTorqueMap& torqueSettingsMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
        const double finalTime,
        const RotationalPropagatorType propagator = quaternions,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                torqueSettingsMap, bodiesToIntegrate, initialBodyStates,
                std::make_shared< PropagationTimeTerminationSettings >( finalTime ), propagator,
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ), printInterval );
}


//! Class for defining settings for propagating the mass of a body
/*!
 *  Class for defining settings for propagating the mass of a body. The body masses are propagated in their natural
 *  form (i.e. no choice of equations of motion as is the case for translational dynamics)l
 */
template< typename StateScalarType >
class MassPropagatorSettings: public SingleArcPropagatorSettings< StateScalarType >
{
public:

    //! Constructor of mass state propagator settings, with single mass rate model per body.
    /*!
     * Constructor  of mass state propagator settings, with single mass rate model per body.
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateModels List of mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > >& massRateModels,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        bodiesWithMassToPropagate_( bodiesWithMassToPropagate )
    {
        for( std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > >::const_iterator
             massRateIterator = massRateModels.begin( ); massRateIterator != massRateModels.end( ); massRateIterator++ )
        {
            massRateModels_[ massRateIterator->first ].push_back( massRateIterator->second );
        }
    }

    //! Constructor of mass state propagator settings, with already-created mass rate models.
    /*!
     * Constructor  of mass state propagator settings, with already-created mass rate models.
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateModels List of mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > >&
            massRateModels,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        bodiesWithMassToPropagate_( bodiesWithMassToPropagate ), massRateModels_( massRateModels )
    { verifyInput( ); }

    //! Constructor of mass state propagator settings, with settings for mass rate models.
    /*!
     * Constructor  of mass state propagator settings, with settings for mass rate models.
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateSettings List of settings for mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const simulation_setup::SelectedMassRateModelMap& massRateSettings,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        bodiesWithMassToPropagate_( bodiesWithMassToPropagate ), massRateSettingsMap_( massRateSettings )
    { verifyInput( ); }

    //! List of bodies for which the mass is to be propagated.
    std::vector< std::string > bodiesWithMassToPropagate_;

    //! Function to create the mass-rate models with support for thrust-acceleration-based mass-rate models.
    /*!
     * Function to create the mass-rate models, with the possibility to specify an acceleration map for setting up
     * mass-rate models determined from thrust accelerations.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     * \param accelerationMap Map of accelerations in the propagation.
     */
    void resetIntegratedStateModels(
            const simulation_setup::SystemOfBodies& bodies,
            const basic_astrodynamics::AccelerationMap& accelerationMap )
    {
        massRateModels_ = simulation_setup::createMassRateModelsMap( bodies, massRateSettingsMap_, accelerationMap );
    }

    //! Function to create the mass-rate models.
    /*!
     * Function to create the mass-rate models.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        resetIntegratedStateModels( bodies, basic_astrodynamics::AccelerationMap( ) );
    }

    //! Function to get the mass-rate settings map.
    /*!
     * Function to get the mass-rate settings map.
     * \return The mass-rate settings map.
     */
    simulation_setup::SelectedMassRateModelMap getMassRateSettingsMap( ) const
    {
        return massRateSettingsMap_;
    }

    //! Function to get the mass-rate models map.
    /*!
     * Function to get the mass-rate models map.
     * \return The mass-rate models map.
     */
    basic_astrodynamics::MassRateModelMap getMassRateModelsMap( ) const
    {
        if ( massRateModels_.size( ) == 0 && massRateSettingsMap_.size( ) != 0 )
        {
            std::cerr << "Inconsistent sizes for map of mass-rate settings and map of mass-rate models. "
                      << "Did you forget to call resetIntegratedStateModels on the propagator?" << std::endl;
        }
        return massRateModels_;
    }


private:

    void verifyInput( )
    {
        if( this->initialStates_.rows( ) != static_cast< int >( bodiesWithMassToPropagate_.size( ) ) )
        {
            throw std::runtime_error( "Error when defining body mass propagator settings, provided initial state size (" +
                                      std::to_string( this->initialStates_.rows( ) ) +
                                      ") is incompatible with list of bodies for which mass is to be propagated (size " +
                                      std::to_string( bodiesWithMassToPropagate_.size( ) ) +
                                      ")");
        }
    }

    //! List of mass rate settings per propagated body.
    simulation_setup::SelectedMassRateModelMap massRateSettingsMap_;

    //! List of mass rate models per propagated body.
    basic_astrodynamics::MassRateModelMap massRateModels_;

};

template< typename StateScalarType = double >
inline std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings(
		const std::vector< std::string > bodiesWithMassToPropagate,
		const std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > >& massRateModels,
		const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
		const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
		const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
		std::shared_ptr< DependentVariableSaveSettings >( ),
		const double printInterval = TUDAT_NAN )
{
	return std::make_shared< MassPropagatorSettings< StateScalarType > >(bodiesWithMassToPropagate,
			  massRateModels, initialBodyMasses, terminationSettings, dependentVariablesToSave,
			  printInterval);
}

template< typename StateScalarType = double >
inline std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings(
		const std::vector< std::string > bodiesWithMassToPropagate,
		const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > >&
				massRateModels,
		const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
		const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
		const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
				std::shared_ptr< DependentVariableSaveSettings >( ),
		const double printInterval = TUDAT_NAN )
{
	return std::make_shared< MassPropagatorSettings< StateScalarType > >(bodiesWithMassToPropagate,
			 massRateModels, initialBodyMasses, terminationSettings, dependentVariablesToSave,
			 printInterval);
}

template< typename StateScalarType = double >
inline std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings(
		const std::vector< std::string > bodiesWithMassToPropagate,
		const simulation_setup::SelectedMassRateModelMap& massRateSettings,
		const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
		const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
		const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
				std::shared_ptr< DependentVariableSaveSettings >( ),
		const double printInterval = TUDAT_NAN )
{
	return std::make_shared< MassPropagatorSettings< StateScalarType > >(bodiesWithMassToPropagate,
			 massRateSettings, initialBodyMasses, terminationSettings, dependentVariablesToSave,
			 printInterval);
}

template< typename StateScalarType = double >
inline std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings(
        const std::vector< std::string > bodiesWithMassToPropagate,
        const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > >&
                massRateModels,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< MassPropagatorSettings< StateScalarType > >(bodiesWithMassToPropagate,
             massRateModels, initialBodyMasses, terminationSettings, std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ),
             printInterval);
}

template< typename StateScalarType = double >
inline std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings(
        const std::vector< std::string > bodiesWithMassToPropagate,
        const simulation_setup::SelectedMassRateModelMap& massRateSettings,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< MassPropagatorSettings< StateScalarType > >(bodiesWithMassToPropagate,
             massRateSettings, initialBodyMasses, terminationSettings, std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ),
             printInterval);
}

//! Function to evaluate a floating point state-derivative function as though it was a vector state function
/*!
 *  Function to evaluate a floating point state-derivative function as though it was a vector state function.
 *  \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
 *  \param currentTime Time at which to evaluate the state derivative function
 *  \param currentStateVector Vector of size 1 containing the current state
 *  \return Current state derivative (as vector of size 1).
 */
template< typename StateScalarType = double, typename TimeType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertScalarToVectorStateFunction(
        const std::function< StateScalarType( const TimeType, const StateScalarType ) > stateDerivativeFunction,
        const TimeType currentTime,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& currentStateVector )
{
    if( currentStateVector.rows( ) != 1 )
    {
        throw std::runtime_error( "Error, expected vector of size one when converting scalar to vector state function" );
    }
    return ( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( 1 )
             << stateDerivativeFunction( currentTime, currentStateVector( 0 ) ) ).finished( );
}


//! Class used to provide settings for a custom state derivative model
/*!
 *  Class used to provide settings for a custom state derivative model. The custom state derivative function has to be
 *  defined by the user and provided as input here. It may depend on the current state, and any dependent variables may
 *  be computed for other state derivative models.
 */
template< typename StateScalarType = double, typename TimeType = double >
class CustomStatePropagatorSettings: public SingleArcPropagatorSettings< StateScalarType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateVectorType;

    //! Constructor for scalar custom state
    /*!
     * Constructor for scalar custom state
     * \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
     * \param initialState Initial state value of custom state
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    CustomStatePropagatorSettings(
            const std::function< StateScalarType( const TimeType, const StateScalarType ) > stateDerivativeFunction,
            const StateScalarType initialState,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >(
            custom_state, ( StateVectorType( 1 ) << initialState ).finished( ), terminationSettings,
            dependentVariablesToSave, printInterval ),
        stateDerivativeFunction_( std::bind( &convertScalarToVectorStateFunction< StateScalarType, TimeType >,
                                             stateDerivativeFunction, std::placeholders::_1, std::placeholders::_2 ) ), stateSize_( 1 )
    { }

    //! Constructor for vector custom state
    /*!
     * Constructor for vector custom state
     * \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
     * \param initialState Initial state value of custom state
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    CustomStatePropagatorSettings(
            const std::function< StateVectorType( const TimeType, const StateVectorType& ) > stateDerivativeFunction,
            const StateVectorType initialState,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >( custom_state, initialState, terminationSettings,
                                                        dependentVariablesToSave, printInterval ),
        stateDerivativeFunction_( stateDerivativeFunction ), stateSize_( initialState.rows( ) ){ }

    //! Destructor
    ~CustomStatePropagatorSettings( ){ }

    //! Function to compute the state derivative, as a function of current time and state.
    std::function< StateVectorType( const TimeType, const StateVectorType& ) > stateDerivativeFunction_;

    //! Size of the state that is propagated.
    int stateSize_;

    //! Function to create the integrated state models. Always throws an error.
    /*!
     * Function to create the integrated state models. This method always throws a `runtime_error` because
     * the integrated state models cannot be created automatically for `CustomStatePropagatorSettings`.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        throw std::runtime_error( "Could not create integrated state models for custom state propagator." );
    }
};

template< typename StateScalarType = double, typename TimeType = double >
inline std::shared_ptr< CustomStatePropagatorSettings< StateScalarType > >
customStatePropagatorSettings(
        const std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialState,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN  )
{
    return std::make_shared< CustomStatePropagatorSettings < StateScalarType > >(
                stateDerivativeFunction, initialState, terminationSettings, std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ), printInterval );
}


//! Function to create multi-arc propagator settings by merging an existing multi-arc with single-arc settings
/*!
 *  Function to create multi-arc propagator settings by merging an existing multi-arc with single-arc settings. The single-arc
 *  settings are converted to multi-arc, and the single-arc propagated bodies are appended at the beginning of the vector
 *  of propagated bodies. Currently, only translational dynamics is supported.
 *  \param singleArcSettings Single-arc settings that are to be added (to the head of the list of propagated bodies) of the input
 *  multi-arc settings.
 *  \param multiArcSettings Multi-arc settings that are to be extended
 *  \param numberofArcs Number of arcs in which the single-arc dynamics is to be split
 *  \return Multi-arc propagator settings by merging an existing multi-arc with single-arc settings
 */
template< typename StateScalarType = double >
std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > getExtendedMultiPropagatorSettings(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcSettings,
        const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcSettings,
        const int numberofArcs )
{
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > constituentSingleArcSettings;

    // Check parameter type
    switch( singleArcSettings->getStateType( ) )
    {
    case translational_state:
    {
        // Check single-arc consistency
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > singleArcTranslationalSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( singleArcSettings );
        if( singleArcTranslationalSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error when making multi-arc propagator settings from single arc. Translational input not consistent." );
        }

        // Create list of single-arc settings
        for( int i = 0; i < numberofArcs; i++ )
        {
            // Check multi-arc consistency
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > currentArcTranslationalSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                        multiArcSettings->getSingleArcSettings( ).at( i ) );
            if( currentArcTranslationalSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error when making multi-arc propagator settings from single arc. Multi-arc input not consistent with single-arc (translational)." );
            }

            // Create full list of central bodies
            std::vector< std::string > multiArcCentralBodies = currentArcTranslationalSettings->centralBodies_;
            std::vector< std::string > fullCentralBodies = singleArcTranslationalSettings->centralBodies_;
            fullCentralBodies.insert( fullCentralBodies.end( ), multiArcCentralBodies.begin( ), multiArcCentralBodies.end( ) );

            // Create full accelerations map
            basic_astrodynamics::AccelerationMap multiArcAccelerationsMap = currentArcTranslationalSettings->getAccelerationsMap( );
            basic_astrodynamics::AccelerationMap fullAccelerationsMap = singleArcTranslationalSettings->getAccelerationsMap( );
            fullAccelerationsMap.insert( multiArcAccelerationsMap.begin( ), multiArcAccelerationsMap.end( ) );

            // Create full list of propagated bodies
            std::vector< std::string > multiArcBodiesToIntegrate = currentArcTranslationalSettings->bodiesToIntegrate_;
            std::vector< std::string > fullBodiesToIntegrate = singleArcTranslationalSettings->bodiesToIntegrate_;
            fullBodiesToIntegrate.insert( fullBodiesToIntegrate.end( ), multiArcBodiesToIntegrate.begin( ), multiArcBodiesToIntegrate.end( ) );

            // Create full initial state list
            int fullSingleArcSize = 6 * fullCentralBodies.size( );
            int singleArcSize = 6 * singleArcTranslationalSettings->centralBodies_.size( );

            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  > currentArcInitialStates =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  >::Zero( fullSingleArcSize );

            // Get single arc initial states (NaN if not first arc)
            if( i == 0 )
            {
                currentArcInitialStates.segment( 0, singleArcSize ) = singleArcTranslationalSettings->getInitialStates( );
            }
            else
            {
                currentArcInitialStates.segment( 0, singleArcSize ) =
                        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  >::Constant( singleArcSize, TUDAT_NAN );
            }

            // Get existing multi-arc initial states
            currentArcInitialStates.segment( singleArcSize, fullSingleArcSize - singleArcSize ) =
                    multiArcSettings->getSingleArcSettings( ).at( i )->getInitialStates( );


            TranslationalPropagatorType propagatorToUse = currentArcTranslationalSettings->propagator_;

            // Retrieve dependent variables that are to be saved.
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > multiArcDependentVariablesToSave;
            if( currentArcTranslationalSettings->getDependentVariablesToSave( ) != nullptr )
            {
                multiArcDependentVariablesToSave  = currentArcTranslationalSettings->getDependentVariablesToSave( )->dependentVariables_;
            }
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > fullDependentVariablesToSave;
            if( singleArcTranslationalSettings->getDependentVariablesToSave( ) != nullptr )
            {
                fullDependentVariablesToSave = singleArcTranslationalSettings->getDependentVariablesToSave( )->dependentVariables_;
            }
            fullDependentVariablesToSave.insert(
                        fullDependentVariablesToSave.end( ), multiArcDependentVariablesToSave.begin( ),
                        multiArcDependentVariablesToSave.end( ) );

            // Create dependent variables object
            std::shared_ptr< DependentVariableSaveSettings > fullDependentVariablesObject;
            if( fullDependentVariablesToSave.size( ) > 0 )
            {
                fullDependentVariablesObject = std::make_shared< DependentVariableSaveSettings >(
                            fullDependentVariablesToSave, true );
            }


            constituentSingleArcSettings.push_back(
                        std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >(
                            fullCentralBodies, fullAccelerationsMap, fullBodiesToIntegrate,
                            currentArcInitialStates,
                            multiArcSettings->getSingleArcSettings( ).at( i )->getTerminationSettings( ), propagatorToUse,
                            fullDependentVariablesObject, singleArcTranslationalSettings->getPrintInterval( ) ) );
        }

        break;
    }
    default:
        throw std::runtime_error( "Error when making multi-arc propagator settings from single arc. Input not recognized." );
    }

    return std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( constituentSingleArcSettings );
}

//! Function to retrieve the state size for a list of propagator settings.
/*!
 *  Function to retrieve the initial state for a list of propagator settings.
 *  \param propagatorSettingsList List of propagator settings (sorted by type as key). Map value provides list
 *  of propagator settings for given type.
 *  \return Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the
 *  vector of SingleArcPropagatorSettings of given type.
 */
template< typename StateScalarType >
int getMultiTypePropagatorStateSize(
        const std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    int stateSize = 0;

    // Iterate over all propagation settings and add size to list
    for( typename std::map< IntegratedStateType,
         std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            stateSize += typeIterator->second.at( i )->getPropagatedStateSize( );
        }
    }
    return stateSize;
}

template< typename StateScalarType >
int getMultiTypePropagatorConventionalStateSize(
        const std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    int stateSize = 0;

    // Iterate over all propagation settings and add size to list
    for( typename std::map< IntegratedStateType,
         std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            stateSize += typeIterator->second.at( i )->getConventionalStateSize( );
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
 *  vector of SingleArcPropagatorSettings of given type.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > createCombinedInitialState(
        const std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    // Get total size of propagated state
    int totalSize = getMultiTypePropagatorConventionalStateSize( propagatorSettingsList );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > combinedInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalSize, 1 );

    // Iterate over all propagation settings and add to total list
    int currentIndex = 0;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState;
    for( typename std::map< IntegratedStateType,
         std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::const_iterator
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
template< typename StateScalarType = double >
class MultiTypePropagatorSettings: public SingleArcPropagatorSettings< StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsMap List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MultiTypePropagatorSettings(
            const std::map< IntegratedStateType,
            std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > > propagatorSettingsMap,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >(
            hybrid, createCombinedInitialState< StateScalarType >( propagatorSettingsMap ),
            terminationSettings, dependentVariablesToSave, printInterval ),
        propagatorSettingsMap_( propagatorSettingsMap )
    {
        this->stateSize_ = getMultiTypePropagatorStateSize( propagatorSettingsMap_ );
    }

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsVector Vector of propagator settings to use.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MultiTypePropagatorSettings(
            const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsVector,
            const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        SingleArcPropagatorSettings< StateScalarType >(
            hybrid, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 0 ),
            terminationSettings, dependentVariablesToSave, printInterval )
    {
        for( unsigned int i = 0; i < propagatorSettingsVector.size( ); i++ )
        {
            propagatorSettingsMap_[ propagatorSettingsVector.at( i )->getStateType( ) ].push_back(
                        propagatorSettingsVector.at( i ) );
        }

        this->initialStates_ = createCombinedInitialState< StateScalarType >( propagatorSettingsMap_ );
        this->stateSize_ = getMultiTypePropagatorStateSize( propagatorSettingsMap_ );

    }

    //! Destructor
    ~MultiTypePropagatorSettings( ){ }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration, sorted in order of
     * IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.
     */
    void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        // Iterate over all propagator settings.
        int currentStartIndex = 0;
        for( typename std::map< IntegratedStateType,
             std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::iterator
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
            std::string errorMessage = "Error when resetting multi-type state, total size is incompatible " +
                    std::to_string( currentStartIndex ) +
                    std::to_string( initialBodyStates.rows( ) );
            throw std::runtime_error( errorMessage );
        }

        this->initialStates_ = createCombinedInitialState< StateScalarType >( propagatorSettingsMap_ );

    }

    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > getSingleTypePropagatorSettings(
              const IntegratedStateType stateType )
      {
          if( propagatorSettingsMap_.count( stateType ) == 0 )
          {
              throw std::runtime_error( "Error when requesting propagator settings of type " +
                                        std::to_string( stateType ) + " from multi-type settings, no such settings found" );
          }

          if( propagatorSettingsMap_.count( stateType ) != 1 )
          {
              throw std::runtime_error( "Error when requesting propagator settings of type " +
                                        std::to_string( stateType ) + " from multi-type settings, multiple settings of given type found" );
          }

          return propagatorSettingsMap_.at( stateType ).at( 0 );
      }

    //! List of propagator settings to use
    /*!
     * List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     */
    std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >
    propagatorSettingsMap_;

    std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >
    getPropagatorSettingsMap( )
    {
        return propagatorSettingsMap_;
    }

    //! Function to create the integrated state models (e.g. acceleration/torque/mass-rate models).
    /*!
     * Function to create the integrated state models (e.g. acceleration/torque/mass-rate models) for
     * each fo the propagators state types contained in `propagatorSettingsMap_`.
     * \param bodies Map of bodies in the propagation, with keys the names of the bodies.
     */
    virtual void resetIntegratedStateModels( const simulation_setup::SystemOfBodies& bodies )
    {
        std::vector< std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > >
                vectorOfTranslationalSettings;
        if ( propagatorSettingsMap_.count( translational_state ) > 0 )
        {
            for ( std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > translationalSettings :
                  propagatorSettingsMap_.at( translational_state ) )
            {
                vectorOfTranslationalSettings.push_back(
                            std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                                translationalSettings ) );
            }
        }

        for ( auto entry : propagatorSettingsMap_ )
        {
            for ( unsigned int i = 0; i < entry.second.size( ); ++i )
            {
                std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcSettings =
                        entry.second.at( i );
                if ( singleArcSettings )
                {
                    std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings
                            = std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >(
                                singleArcSettings );
                    if ( massPropagatorSettings && vectorOfTranslationalSettings.size( ) > 0 )
                    {
                        if ( entry.second.size( ) != vectorOfTranslationalSettings.size( ) )
                        {
                            throw std::runtime_error( "Could not create integrated state model for mass "
                                                      "propagator settings because a one-to-one relationship "
                                                      "between the specified mass-rate propagators and the "
                                                      "tranlational propagators could not be inferred. Create "
                                                      "the models manually or provide one translational "
                                                      "propagator settings for each mass rate propagator settings, "
                                                      "or provide no translational propagator settings at all." );
                        }
                        massPropagatorSettings->resetIntegratedStateModels(
                                    bodies, vectorOfTranslationalSettings.at( i )->getAccelerationsMap( ) );
                    }
                    else
                    {
                        singleArcSettings->resetIntegratedStateModels( bodies );
                    }
                }
            }
        }
    }

};

template< typename StateScalarType = double >
inline std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings(
		const std::map< IntegratedStateType,
		std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > > propagatorSettingsMap,
		const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
		const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
		std::shared_ptr< DependentVariableSaveSettings >( ),
		const double printInterval = TUDAT_NAN )
{
	return std::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
			propagatorSettingsMap, terminationSettings, dependentVariablesToSave, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings(
		const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsVector,
		const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
		const std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
		std::shared_ptr< DependentVariableSaveSettings >( ),
		const double printInterval = TUDAT_NAN )
{
	return std::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
			propagatorSettingsVector, terminationSettings, dependentVariablesToSave, printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings(
        const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsVector,
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
            propagatorSettingsVector, terminationSettings, std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ),
                printInterval );
}

template< typename StateScalarType = double >
inline std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings(
        const std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsVector,
        const double finalTime,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ),
        const double printInterval = TUDAT_NAN )
{
    return std::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
            propagatorSettingsVector, std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
                std::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave ),
                printInterval );
}



extern template class PropagatorSettings< double >;
extern template class SingleArcPropagatorSettings< double >;
extern template class MultiArcPropagatorSettings< double >;
extern template class TranslationalStatePropagatorSettings< double >;
extern template class RotationalStatePropagatorSettings< double >;
extern template class MassPropagatorSettings< double >;
extern template class CustomStatePropagatorSettings< double >;
extern template class MultiTypePropagatorSettings< double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class PropagatorSettings< long double >;
extern template class SingleArcPropagatorSettings< long double >;
extern template class MultiArcPropagatorSettings< long double >;
extern template class TranslationalStatePropagatorSettings< long double >;
extern template class RotationalStatePropagatorSettings< long double >;
extern template class MassPropagatorSettings< long double >;
extern template class MultiTypePropagatorSettings< long double >;
extern template class CustomStatePropagatorSettings< long double >;
extern template class CustomStatePropagatorSettings< double, Time >;
extern template class CustomStatePropagatorSettings< long double, Time >;
#endif

//! Function to retrieve list of accelerations from propagator settings
/*!
 *  Function to retrieve list of accelerations from propagator settings. Extracts the translational dynamics elements, and
 *  the associated acceleration models
 *  \param singleArcPropagatorSettings Propagator settings
 *  \return List of acceleration models
 */
template< typename StateScalarType = double >
basic_astrodynamics::AccelerationMap getAccelerationMapFromPropagatorSettings(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcPropagatorSettings )
{
    basic_astrodynamics::AccelerationMap accelerationMap;

    switch( singleArcPropagatorSettings->getStateType( ) )
    {
    case hybrid:
    {
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
        if( multiTypePropagatorSettings->propagatorSettingsMap_.count( translational_state ) > 0 )
        {
            if( multiTypePropagatorSettings->propagatorSettingsMap_.at( translational_state ).size( ) == 1 )
            {
                accelerationMap = getAccelerationMapFromPropagatorSettings(
                            multiTypePropagatorSettings->propagatorSettingsMap_.at( translational_state ).at( 0 ) );
            }
            else
            {
                throw std::runtime_error( "Error when get accelerations map from propagator settings, multi-type cannot be parsed properly" );
            }
        }

        break;
    }
    case translational_state:
    {
        accelerationMap = std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                    singleArcPropagatorSettings )->getAccelerationsMap( );

        break;
    }
    default:
        break;
    }

    return accelerationMap;
}

//! Function to retrieve the list of integrated state types and reference ids
/*!
* Function to retrieve the list of integrated state types and reference ids. For translational and rotational dynamics,
* the id refers only to the body being propagated (and the second entry of the pair is empty: ""). For proper time
* propagation, a body and a reference point may be provided, resulting in non-empty first and second pair entries.
* \param propagatorSettings Settings that are to be used for the propagation.
* \return List of integrated state types and reference ids
*/
template< typename StateScalarType >
std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStateList;

    // Identify propagator type
    switch( propagatorSettings->getStateType( ) )
    {
    case hybrid:
    {
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > singleTypeIntegratedStateList;


        for( typename std::map< IntegratedStateType,
             std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::const_iterator
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
                        std::string errorMessage = "Error when making integrated state list for hybrid propagator, inconsistency encountered " +
                                std::to_string( singleTypeIntegratedStateList.begin( )->first ) + " " +
                                std::to_string( typeIterator->first ) + " " +
                                std::to_string( singleTypeIntegratedStateList.size( ) ) + " " +
                                std::to_string( singleTypeIntegratedStateList.begin( )->second.size( ) );
                        throw std::runtime_error( errorMessage );
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
                throw std::runtime_error( "Error when making integrated state list, cannot handle hybrid propagator inside hybrid propagator" );
            }
        }
        break;
    }
    case translational_state:
    {
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationalPropagatorSettings = std::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        if( translationalPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error getting integrated state type list, translational state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < translationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair( translationalPropagatorSettings->bodiesToIntegrate_.at( i ), "" ) );
        }
        integratedStateList[ translational_state ] = integratedBodies;

        break;
    }
    case rotational_state:
    {
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < rotationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair( rotationalPropagatorSettings->bodiesToIntegrate_.at( i ), "" ) );
        }

        integratedStateList[ rotational_state ] = integratedBodies;

        break;
    }
    case body_mass_state:
    {
        std::shared_ptr< MassPropagatorSettings< StateScalarType > >
                massPropagatorSettings = std::dynamic_pointer_cast<
                MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == nullptr )
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
    case custom_state:
    {
        std::vector< std::pair< std::string, std::string > > customList;
        customList.push_back( std::make_pair( "N/A", "N/A" ) );
        integratedStateList[ custom_state ] = customList;
        break;
    }
    default:
        throw std::runtime_error( "Error, could not process integrated state type in getIntegratedTypeAndBodyList " +
                                  std::to_string( propagatorSettings->getStateType( ) ) );
    }

    return integratedStateList;
}


// addition for thesis work (Jonas Hener), considered generally useful
template< typename StateScalarType >
void addDepedentVariableSettings(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToAdd,
        const std::shared_ptr< propagators::SingleArcPropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave = propagatorSettings->getDependentVariablesToSave( );
    if( dependentVariablesToSave != nullptr )
    {
        dependentVariablesToSave->dependentVariables_.insert(
                dependentVariablesToSave->dependentVariables_.end( ), dependentVariablesToAdd.begin( ), dependentVariablesToAdd.end( ) );
    }
    else
    {
        dependentVariablesToSave = std::make_shared< DependentVariableSaveSettings >(
                dependentVariablesToAdd, false );
        propagatorSettings->resetDependentVariablesToSave( dependentVariablesToSave );
    }
}




template< typename StateScalarType >
void resetSingleArcInitialStates(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const std::map< propagators::IntegratedStateType, std::map< std::pair< std::string, std::string >,
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > currentArcInitialStates )
{
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    switch( propagatorSettings->getStateType( ) )
    {
    case translational_state:
    {
        if( currentArcInitialStates.count( translational_state ) == 0 )
        {
            throw std::runtime_error( "Error when resetting initial translational state from sorted data, no data found." );
        }
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalStateSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        std::vector< std::string > propagatedBodies = translationalStateSettings->bodiesToIntegrate_;

        if( currentArcInitialStates.at( translational_state ).size( ) != propagatedBodies.size( ) )
        {
            throw std::runtime_error( "Error when resetting initial translational state from sorted data, body list size is incompatible." );
        }

        VectorType totalInitialState = VectorType( 6 * propagatedBodies.size( ) );

        for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
        {
            if( currentArcInitialStates.at( translational_state ).count( std::make_pair( propagatedBodies.at( i ), "" ) ) == 0 )
            {
                throw std::runtime_error(
                            "Error when resetting initial translational state from sorted data, did not find body " + propagatedBodies.at( i ) );
            }

            totalInitialState.segment( i * 6, 6 ) =
                    currentArcInitialStates.at( translational_state ).at( std::make_pair( propagatedBodies.at( i ), "" ) );
        }
        translationalStateSettings->resetInitialStates( totalInitialState );

        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize state type " + std::to_string( propagatorSettings->getStateType( ) ) +
                                  " when resetting initial states from parameter data " );
    }

}


extern template std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList< double >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONSETTINGS_H
