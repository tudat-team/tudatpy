/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NBODYSTATEDERIVATIVE_H
#define TUDAT_NBODYSTATEDERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/propagators/centralBodyData.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace propagators
{

// Enum listing propagator types for translational dynamics that can be used.
//! @get_docstring(TranslationalPropagatorType.__docstring__)
enum TranslationalPropagatorType
{
    undefined_translational_propagator = -1,
    cowell = 0,
    encke = 1,
    gauss_keplerian = 2,
    gauss_modified_equinoctial = 3,
    unified_state_model_quaternions = 4,
    unified_state_model_modified_rodrigues_parameters = 5,
    unified_state_model_exponential_map = 6
};

// Function to remove the central gravity acceleration from an AccelerationMap
/*
 * Function to remove the central gravity acceleration from an AccelerationMap. This is crucial for propagation methods in
 * which the deviation from a reference Kepler orbit is propagated. If the central gravity is a spherical harmonic
 * acceleration, the point mass term is removed by setting the C(0,0) coefficnet to 0
 *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
 *  \param centralBodies List of names of bodies of which the central terms are to be removed
 *  (per entry of bodiesToIntegrate)
 *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
 *  body, identifying the body being acted on and the body acted on by an acceleration. The map
 *  has as key a string denoting the name of the body the list of accelerations, provided as the
 *  value corresponding to a key, is acting on.  This map-value is again a map with string as
 *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
 *  model.
 * \return Functions returning the gravitational parameters of the central terms that were removed.
 */
std::vector< std::function< double( ) > > removeCentralGravityAccelerations(
        const std::vector< std::string >& centralBodies, const std::vector< std::string >& bodiesToIntegrate,
        basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
        std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > >& removedAcceleration  );

// Function to determine in which order the ephemerides are to be updated
/*
 * Function to determine in which order the ephemerides are to be updated. The order depends on the
 * dependencies between the ephemeris/integration origins.
 * \param integratedBodies List of bodies that are numerically integrated.
 * \param centralBodies List of origins w.r.t. the integratedBodies' translational dynamics is propagated.
 * \param ephemerisOrigins Origin of the Ephemeris objects of the integratedBodies.
 * \return
 */
std::vector< std::string > determineEphemerisUpdateorder( std::vector< std::string > integratedBodies,
                                                          std::vector< std::string > centralBodies,
                                                          std::vector< std::string > ephemerisOrigins );

// State derivative for the translational dynamics of N bodies
/*
 * This class calculates the trabnslational state derivative of any
 * number of bodies, each under the influence of any number of bodies,
 * both from the set being integrated and otherwise.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

    // Constructor from data for translational Cartesian state derivative creation.
    // It is assumed that all acceleration are exerted on bodies by bodies.
    /*
     *  From this constructor, the object for generating the state derivative is created. Required
     *  are the acceleration models, a map of all (named) bodies involved in the simulation and a
     *  list of body names, which must be a subset of the bodyList that are to be numerically
     *  integrated. Note that the state derivative model currently has 3 degrees of freedom (3
     *  translational) in Cartesian coordinates.
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param propagatorType Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                          const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                          const TranslationalPropagatorType propagatorType,
                          const std::vector< std::string >& bodiesToIntegrate,
                          const bool removeCentralTerm = false ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::translational_state ),
        accelerationModelsPerBody_( accelerationModelsPerBody ),
        removedCentralAccelerations_( std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > >( ) ),
        updateRemovedAccelerations_( std::vector< std::string >( ) ),
        centralBodyData_( centralBodyData ),
        propagatorType_( propagatorType ),
        bodiesToBeIntegratedNumerically_( bodiesToIntegrate ),
        removeCentralTerm_( removeCentralTerm )
    {
        originalAccelerationModelsPerBody_ = this->accelerationModelsPerBody_ ;

        // Add empty acceleration map if body is to be propagated with no accelerations.
        for( unsigned int i = 0; i < bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            if( accelerationModelsPerBody_.count( bodiesToBeIntegratedNumerically_.at( i ) ) == 0 )
            {
                accelerationModelsPerBody_[ bodiesToBeIntegratedNumerically_.at( i ) ] =
                        basic_astrodynamics::SingleBodyAccelerationMap( );
            }
        }

        // Correct order of propagated bodies.
        for( outerAccelerationIterator = accelerationModelsPerBody_.begin( );
             outerAccelerationIterator != accelerationModelsPerBody_.end( );
             outerAccelerationIterator++ )
        {
            std::vector< std::string >::iterator findIterator =
                    std::find( bodiesToBeIntegratedNumerically_.begin( ), bodiesToBeIntegratedNumerically_.end( ),
                               outerAccelerationIterator->first );
            bodyOrder_.push_back( std::distance( bodiesToBeIntegratedNumerically_.begin( ), findIterator ) );
        }

        createAccelerationModelList( );
        verifyInput( );
    }

    // Destructor
    virtual ~NBodyStateDerivative( ){ }

    // Function to clear any reference/cached values of state derivative model
    /*
     * Function to clear any reference/cached values of state derivative model, in addition to those performed in the
     * clearTranslationalStateDerivativeModel function. Default implementation is empty.
     */
    virtual void clearDerivedTranslationalStateDerivativeModel( ){ }

    // Function to clear reference/cached values of acceleration models
    /*
     * Function to clear reference/cached values of acceleration models, to ensure that they are all recalculated.
     */
    void clearTranslationalStateDerivativeModel( )
    {
        for( unsigned int i = 0; i < accelerationModelList_.size( ); i++ )
        {
            accelerationModelList_.at( i )->resetTime( TUDAT_NAN );
        }        

        for( unsigned int i = 0; i < updateRemovedAccelerations_.size( ); i++ )
        {
            if( removedCentralAccelerations_.count( updateRemovedAccelerations_.at( i  ) ) > 0 )
            {
                removedCentralAccelerations_[ updateRemovedAccelerations_.at( i ) ]->resetTime( TUDAT_NAN );
            }
        }
    }

    // Function to clear reference/cached values of translational state derivative model
    /*
     * Function to clear reference/cached values of translational state derivative model. For each derived class, this
     * entails resetting the current time in the acceleration models to NaN (see clearTranslationalStateDerivativeModel).
     * Every derived class requiring additional values to be cleared should implement the
     * clearDerivedTranslationalStateDerivativeModel function.
     */
    void clearStateDerivativeModel(  )
    {
        clearTranslationalStateDerivativeModel( );
        clearDerivedTranslationalStateDerivativeModel( );
    }

    // Function to update the state derivative model to the current time.
    /*
     * Function to update the state derivative model (i.e. acceleration models) to the
     * current time. Note that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function.
     * \param currentTime Time at which state derivative is to be calculated
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {
        for( unsigned int i = 0; i < accelerationModelList_.size( ); i++ )
        {
            accelerationModelList_.at( i )->updateMembers( currentTime );
        }

        for( unsigned int i = 0; i < updateRemovedAccelerations_.size( ); i++ )
        {
            if( removedCentralAccelerations_.count( updateRemovedAccelerations_.at( i  ) ) > 0 )
            {
                removedCentralAccelerations_[ updateRemovedAccelerations_.at( i ) ]->updateMembers( currentTime );
            }
        }
    }

    // Function to convert the propagator-specific form of the state to the conventional form in the global frame.
    /*
     * Function to convert the propagator-specific form of the state to the conventional form in the
     * global frame.  The conventional form for translational dynamics this is the Cartesian
     * position and velocity).  The inertial frame is typically the barycenter with J2000/ECLIPJ2000
     * orientation, but may differ depending on simulation settings.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration).
     * \param time Current time at which the state is valid.
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the Cartesian state in inertial coordinates
     * (returned by reference).
     */
    void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        this->convertToOutputSolution( internalSolution, time, currentCartesianLocalSoluton );

        centralBodyData_->getReferenceFrameOriginInertialStates(
                    currentCartesianLocalSoluton, time, centralBodyStatesWrtGlobalOrigin_, true );

        for( unsigned int i = 0; i < centralBodyStatesWrtGlobalOrigin_.size( ); i++ )
        {
            currentCartesianLocalSoluton.block( i * 6, 0, 6, 1 ) += centralBodyStatesWrtGlobalOrigin_[ i ];
        }
    }

    // Function to get list of names of bodies that are to be integrated numerically.
    /*
     * Function to get list of names of bodies that are to be integrated numerically.
     * \return List of names of bodies that are to be integrated numerically.
     */
    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToBeIntegratedNumerically_;
    }

    // Function to get object providing the current integration origins
    /*
     * Function to get object responsible for providing the current integration origins from the
     * global origins.
     * \return Object providing the current integration origins from the global origins.
     */
    std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > getCentralBodyData( )
    {
        return centralBodyData_;
    }

    // Function to get type of propagator that is to be used (i.e. Cowell, Encke, etc.)
    /*
     * Function to type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     * \return Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     */
    TranslationalPropagatorType getTranslationalPropagatorType( )
    {
        return propagatorType_;
    }

    // Function to return the size of the state handled by the object
    /*
     * Function to return the size of the state handled by the object
     * \return Size of the state under consideration (6 times the number if integrated bodies).
     */
    int getConventionalStateSize( )
    {
        return 6 * bodiesToBeIntegratedNumerically_.size( );
    }

    // Function to retrieve the total acceleration acting on a given body.
    /*
     * Function to retrieve the total acceleration acting on a given body. The environment
     * and acceleration models must have been updated to the current state before calling this
     * function. NOTE: This function is typically used to retrieve the acceleration for output purposes, not to compute the
     * translational state derivative.
     * \param bodyName Name of body for which accelerations are to be retrieved.
     * \return
     */
    Eigen::Vector3d getTotalAccelerationForBody(
            const std::string& bodyName )
    {
        // Check if body is propagated.
        Eigen::Vector3d totalAcceleration = Eigen::Vector3d::Zero( );
        if( std::find( bodiesToBeIntegratedNumerically_.begin( ),
                       bodiesToBeIntegratedNumerically_.end( ),
                       bodyName ) == bodiesToBeIntegratedNumerically_.end( ) )
        {
            std::string errorMessage = "Error when getting total acceleration for body " + bodyName +
                    ", no such acceleration is found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
//            if( removedCentralAcceleration_ != nullptr && updateRemovedAcceleration_ == false )
//            {
//                std::string errorMessage = "Error when getting total acceleration for body " + bodyName +
//                        ", central term is removed, but cannot be evaluated.";
//                throw std::runtime_error( errorMessage );
//            }

            if( originalAccelerationModelsPerBody_.count( bodyName ) != 0 )
            {
                basic_astrodynamics::SingleBodyAccelerationMap accelerationsOnBody =
                        originalAccelerationModelsPerBody_.at( bodyName );

                // Iterate over all accelerations acting on body
                for( innerAccelerationIterator  = accelerationsOnBody.begin( );
                     innerAccelerationIterator != accelerationsOnBody.end( );
                     innerAccelerationIterator++ )
                {
                    for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                    {
                        // Calculate acceleration and add to state derivative.
                         innerAccelerationIterator->second[ j ]->addCurrentAcceleration( totalAcceleration );
                    }
                }
            }
        }
        return totalAcceleration;
    }

    // Function to retrieve the map containing the list of accelerations acting on each body.
    /*
     * Function to retrieve the map containing the list of accelerations acting on each body.
     * \return Map containing the list of accelerations acting on each body,
     */
    basic_astrodynamics::AccelerationMap getAccelerationsMap( )
    {
        return accelerationModelsPerBody_;
    }

    basic_astrodynamics::AccelerationMap getFullAccelerationsMap( )
    {
        return originalAccelerationModelsPerBody_;
    }

    std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > > getRemovedCentralAcceleration( )
    {
        return removedCentralAccelerations_;
    }

    void setUpdateRemovedAcceleration( const std::string bodyName )
    {
        if( removedCentralAccelerations_.count( bodyName ) > 0 )
        {
            updateRemovedAccelerations_.push_back( bodyName );
        }
    }



protected:

    void verifyInput( )
    {
        for( unsigned int i = 0; i < bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            if( accelerationModelsPerBody_.count( bodiesToBeIntegratedNumerically_.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error, requested propagation of translational dynamics of body " +
                                          bodiesToBeIntegratedNumerically_.at( i ) +
                                          ", but no acceleration models provided" );
            }
        }

        for( auto it : accelerationModelsPerBody_ )
        {
            if( std::find( bodiesToBeIntegratedNumerically_.begin( ),
                           bodiesToBeIntegratedNumerically_.end( ),
                           it.first ) == bodiesToBeIntegratedNumerically_.end( ) )
            {
                throw std::runtime_error( "Error, provided acceleration models for body " +
                                          it.first +
                                          ", but this body is not included in list of bodies for which translational dynamics is to be propagated." );
            }
        }
    }


    // Function to set the vector of acceleration models (accelerationModelList_) form the map of map of
    // acceleration models (accelerationModelsPerBody_).
    void createAccelerationModelList( )
    {
        // Iterate over all accelerations and update their internal state.
        accelerationModelList_.clear( );
        for( outerAccelerationIterator = accelerationModelsPerBody_.begin( );
             outerAccelerationIterator != accelerationModelsPerBody_.end( ); outerAccelerationIterator++ )
        {
            // Iterate over all accelerations acting on body
            for( innerAccelerationIterator  = outerAccelerationIterator->second.begin( );
                 innerAccelerationIterator != outerAccelerationIterator->second.end( );
                 innerAccelerationIterator++ )
            {
                // Update accelerations
                for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                {
                    accelerationModelList_.push_back( innerAccelerationIterator->second.at( j ) );
                }
            }
        }
    }

    // Function to get the state derivative of the system in Cartesian coordinates.
    /*
     * Function to get the state derivative of the system in Cartesian coordinates. The environment
     * and acceleration models must have been updated to the current state before calling this
     * function.
     * \param stateOfSystemToBeIntegrated Current Cartesian state of the system.
     * \param stateDerivative State derivative of the system in Cartesian coordinates (returned by reference).
     * \param addPositionDerivatives Boolean denoting whether the derivatives of the position (e.g. velocity) are to be added
     * to the state derivative vector.
     */
    void sumStateDerivativeContributions(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative,
            const bool addPositionDerivatives = true )
    {
        using namespace basic_astrodynamics;

        stateDerivative.setZero( );

        int currentBodyIndex = 0;
        int currentAccelerationIndex = 0;

        // Iterate over all bodies with accelerations.
        for( outerAccelerationIterator = accelerationModelsPerBody_.begin( );
             outerAccelerationIterator != accelerationModelsPerBody_.end( );
             outerAccelerationIterator++ )
        {
            currentBodyIndex = bodyOrder_[ currentAccelerationIndex ];

            // Iterate over all accelerations acting on body
            for( innerAccelerationIterator  = outerAccelerationIterator->second.begin( );
                 innerAccelerationIterator != outerAccelerationIterator->second.end( );
                 innerAccelerationIterator++ )
            {
                for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                {
                    // Calculate acceleration and add to state derivative.
                    stateDerivative.block( currentBodyIndex * 6 + 3, 0, 3, 1 ) +=
                                                ( innerAccelerationIterator->second[ j ]->getAccelerationReference( ) ).
                                                template cast< StateScalarType >( );

                }
            }

            if( addPositionDerivatives )
            {
                // Add body velocity as derivative of its position.
                stateDerivative.block( currentBodyIndex * 6, 0, 3, 1 ) =
                        ( stateOfSystemToBeIntegrated.segment( currentBodyIndex * 6 + 3, 3 ) );
            }
            currentAccelerationIndex++;
        }
    }

    // Function to get the state derivative of the system in Cartesian coordinates.
    /*
     * Function to get the state derivative of the system in Cartesian coordinates. The environment
     * and acceleration models must have been updated to the current state before calling this
     * function.
     * \param stateOfSystemToBeIntegrated Current Cartesian state of the system.
     * \param stateDerivative State derivative of the system in Cartesian coordinates (returned by reference).
     * \param addPositionDerivatives Boolean denoting whether the derivatives of the position (e.g. velocity) are to be added
     * to the state derivative vector.
     */
    void sumStateDerivativeContributions(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& stateDerivative,
            const bool addPositionDerivatives = true )
    {
        return sumStateDerivativeContributions(
                    stateOfSystemToBeIntegrated,
                    stateDerivative.block( 0, 0, stateDerivative.rows( ), stateDerivative.cols( ) ),
                    addPositionDerivatives );
    }

    // A map containing the list of accelerations acting on each body,
    /*
     * A map containing the list of accelerations acting on each body, identifying the body being
     * acted on and the body acted on by an acceleration. The map has as key a string denoting the
     * name of the body the list of accelerations, provided as the value corresponding to a key, is
     * acting on.  This map-value is again a map with string as key, denoting the body exerting the
     * acceleration, and as value a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationModelsPerBody_;

    basic_astrodynamics::AccelerationMap originalAccelerationModelsPerBody_;

    std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > > removedCentralAccelerations_;

    std::vector< std::string > updateRemovedAccelerations_;

    // Vector of acceleration models, containing all entries of accelerationModelsPerBody_.
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > accelerationModelList_;

    // Object responsible for providing the current integration origins from the global origins.
    std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData_;

    // Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
    TranslationalPropagatorType propagatorType_;

    // List of names of bodies that are to be integrated numerically.
    std::vector< std::string > bodiesToBeIntegratedNumerically_;

    std::vector< int > bodyOrder_;

    // Predefined iterator to save (de-)allocation time.
    std::unordered_map< std::string, std::vector<
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >::iterator innerAccelerationIterator;

    // Predefined iterator to save (de-)allocation time.
    std::unordered_map< std::string, std::unordered_map< std::string, std::vector<
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > > >::iterator outerAccelerationIterator;

    // List of states of the central bodies of the propagated bodies.
    std::vector< Eigen::Matrix< StateScalarType, 6, 1 >  > centralBodyStatesWrtGlobalOrigin_;

    Eigen::Vector3d currentAccelerationComponent_;

    bool removeCentralTerm_;

};

extern template class NBodyStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyStateDerivative< long double, double >;
extern template class NBodyStateDerivative< double, Time >;
extern template class NBodyStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYSTATEDERIVATIVE_H
