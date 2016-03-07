#ifndef TUDAT_NBODYSTATEDERIVATIVE_H
#define TUDAT_NBODYSTATEDERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/centralBodyData.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"


namespace tudat
{

namespace propagators
{

//! State derivative for the translational dynamics of N bodies
/*!
 * This class calculates the trabnslational state derivative of any number of bodies, each under the influence of any
 * number of bodies, both from the set being integrated and otherwise.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;


    //! Constructor from data for translational Cartesian state derivative creation.
    //! It is assumed that all acceleration are exerted on bodies by bodies.
    /*!
     *  From this constructor, the object for generating the state derivative is created. Required are the
     *  acceleration models, a map of all (named) bodies involved in the simulation and  a list of body names,
     *  which must be a subset of the bodyList that are to be numerically integrated. Note that the state derivative
     *  model currently has 3 degrees of freedom (3 translational) in Cartesian coordinates.
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     *  \param centralBodyData Object responsible for providing the current integration origins from the global origins.
     *  \param propagatorType Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                          const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                          const TranslationalPropagatorType propagatorType,
                          const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::transational_state ),
        accelerationModelsPerBody_( accelerationModelsPerBody ),
        centralBodyData_( centralBodyData ),
        propagatorType_( propagatorType ),
        bodiesToBeIntegratedNumerically_( bodiesToIntegrate )
    { }

    //! Destructor
    virtual ~NBodyStateDerivative( ){ }

    //! Function to update the state derivative model to the current time.
    /*!
     * Function to update the state derivative model (i.e. acceleration, torque, etc. models) to the current time. Note that
     * this function only updates the state derivative model itself, the environment models must be updated before calling
     * this function
     * \param currentTime Time at which state derivative is to be calculated
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {
        // Reser all acceleration times (to allow multiple evaluations at same time, e.g. stage 2 and 3 in RK4 integrator)
        for( accelerationMapIterator = accelerationModelsPerBody_.begin( );
             accelerationMapIterator != accelerationModelsPerBody_.end( ); accelerationMapIterator++ )
        {
            // Iterate over all accelerations acting on body
            for( innerAccelerationIterator  = accelerationMapIterator->second.begin( ); innerAccelerationIterator !=
                 accelerationMapIterator->second.end( ); innerAccelerationIterator++ )
            {
                // Update accelerations
                for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                {

                    innerAccelerationIterator->second[ j ]->resetTime( TUDAT_NAN );
                }
            }
        }

        // Iterate over all accelerations and update their internal state.
        for( accelerationMapIterator = accelerationModelsPerBody_.begin( );
             accelerationMapIterator != accelerationModelsPerBody_.end( ); accelerationMapIterator++ )
        {
            // Iterate over all accelerations acting on body
            for( innerAccelerationIterator  = accelerationMapIterator->second.begin( ); innerAccelerationIterator !=
                 accelerationMapIterator->second.end( ); innerAccelerationIterator++ )
            {
                // Update accelerations
                for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                {

                    innerAccelerationIterator->second[ j ]->updateMembers( currentTime );
                }
            }
        }
    }

    //! Function to convert the propagator-specific form of the state to the conventional form in the global frame.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form in the global frame.
     * The conventional form for translational dynamics this is the Cartesian position and velocity).
     * The inertial frame is typically the barycenter with J2000/ECLIPJ2000 orientation, but may differ depending on
     * simulation settings.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in numerical integration).
     * \param time Current time at which the state is valid.
     * \return State (internalSolution), converted to the Cartesian state in inertial coordinates.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > cartesianLocalSolution = this->convertToOutputSolution( internalSolution, time );

        std::vector< Eigen::Matrix< StateScalarType, 6, 1 >  > centralBodyInertialStates =
                centralBodyData_->getReferenceFrameOriginInertialStates( cartesianLocalSolution, time, true );
        for( unsigned int i = 0; i < centralBodyInertialStates.size( ); i++ )
        {
            cartesianLocalSolution.segment( i * 6, 6 ) += centralBodyInertialStates[ i ];
        }
        return cartesianLocalSolution;
    }

    //! Function to get list of names of bodies that are to be integrated numerically.
    /*!
     * Function to get list of names of bodies that are to be integrated numerically.
     * \return List of names of bodies that are to be integrated numerically.
     */
    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToBeIntegratedNumerically_;
    }

    //! Function to get map containing the list of accelerations acting on each body,
    /*!
     * Function to get map containing the list of accelerations acting on each body,
     * \return A map containing the list of accelerations acting on each body,
     */
    basic_astrodynamics::AccelerationMap getAccelerationsMap( )
    {
        return accelerationModelsPerBody_;
    }

    //! Function to get object responsible for providing the current integration origins from the global origins.
    /*!
     * Function to get object responsible for providing the current integration origins from the global origins.
     * \return Object responsible for providing the current integration origins from the global origins.
     */
    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > getCentralBodyData( )
    {
        return centralBodyData_;
    }

    //! Function to get type of propagator that is to be used (i.e. Cowell, Encke, etc.)
    /*!
     * Function to type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     * \return Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
     */
    TranslationalPropagatorType getPropagatorType( )
    {
        return propagatorType_;
    }

    //! Function to return the size of the state handled by the object
    /*!
     * Function to return the size of the state handled by the object
     * \return Size of the state under consideration (6 times the number if integrated bodies).
     */
    int getStateSize( )
    {
        return 6 * bodiesToBeIntegratedNumerically_.size( );
    }

protected:

    //! Function to get the state derivative of the system in Cartesian coordinates.
    /*!
     * Function to get the state derivative of the system in Cartesian coordinates. The environment and acceleration models
     * must have been updated to the current state before calling this function.
     * \param stateOfSystemToBeIntegrated Current Cartesian state of the system.
     * \param time Time at which the state derivative is to be computed
     * \return State derivative of the system in Cartesian coordinates.
     */
    Eigen::VectorXd sumStateDerivativeContributions(
            const Eigen::VectorXd& stateOfSystemToBeIntegrated, const TimeType time )
    {
        using namespace basic_astrodynamics;

        // Declare and initialize to zero vector to be returned.
        Eigen::VectorXd stateDerivative = Eigen::VectorXd::Zero( stateOfSystemToBeIntegrated.size( ) );

        // Iterate over all bodies with accelerations.
        for( unsigned int i = 0; i < bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            // If body undergoes acceleration, calculate and add accelerations.
            if( accelerationModelsPerBody_.count( bodiesToBeIntegratedNumerically_[ i ] ) != 0 )
            {

                {
                    // Iterate over all accelerations acting on body
                    for( innerAccelerationIterator  = accelerationModelsPerBody_[ bodiesToBeIntegratedNumerically_[ i ] ].begin( );
                         innerAccelerationIterator != accelerationModelsPerBody_[ bodiesToBeIntegratedNumerically_[ i ] ].end( );
                         innerAccelerationIterator++ )
                    {
                        for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                        {
                            // Calculate acceleration and add to state derivative.
                            stateDerivative.segment( i * 6 + 3, 3 ) += (
                                        innerAccelerationIterator->second[ j ]->getAcceleration( ) );
                        }
                    }
                }
            }
            // Add body velocity as derivative of its position.
            stateDerivative.segment( i * 6, 3 ) = stateOfSystemToBeIntegrated.segment( i * 6 + 3, 3 );
        }

        return stateDerivative;
    }


    //! A map containing the list of accelerations acting on each body,
    /*!
     * A map containing the list of accelerations acting on each body, identifying
     * the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     * the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     * This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     * a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationModelsPerBody_;

    //! Object responsible for providing the current integration origins from the global origins.
    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData_;

    //! Type of propagator that is to be used (i.e. Cowell, Encke, etc.)
    TranslationalPropagatorType propagatorType_;

    //! List of names of bodies that are to be integrated numerically.
    std::vector< std::string > bodiesToBeIntegratedNumerically_;

    //! Predefined iterator to save (de-)allocation time.
    basic_astrodynamics::AccelerationMap::iterator accelerationMapIterator;

    //! Predefined iterator to save (de-)allocation time.
    std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >::iterator
    innerAccelerationIterator;

};

}

}
#endif // TUDAT_NBODYSTATEDERIVATIVE_H
