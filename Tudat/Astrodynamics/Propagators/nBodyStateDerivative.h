#ifndef NBODYSTATEDERIVATIVE_H
#define NBODYSTATEDERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/centralBodyData.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"


namespace tudat
{

namespace propagators
{

/*!
 * This class calculates the state derivative of any number of bodies each under the influence of any
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
     *  \param bodyList. Map of string and pointers to body objects, providing the named list of bodies in simulation.
     *  \param bodiesToIntegrate. List of names of bodies that are to be integrated numerically.
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

    virtual ~NBodyStateDerivative( ){ }

    void updateStateDerivativeModel( const TimeType currentTime )
    {
        // Reser all acceleration times (to allow multiple evaluations at same time, e.g. stage 2 and 3 in RK4 integrator)
        for( accelerationMapIterator = accelerationModelsPerBody_.begin( );
             accelerationMapIterator != accelerationModelsPerBody_.end( ); accelerationMapIterator++ )
        {
            // Retrieve list of accelerations
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
        }

        // Iterate over all accelerations and update their internal state.
        for( accelerationMapIterator = accelerationModelsPerBody_.begin( );
             accelerationMapIterator != accelerationModelsPerBody_.end( ); accelerationMapIterator++ )
        {
            // Retrieve list of accelerations
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
    }

    //! Calculates the state derivative of the translational motion of the system.
    /*!
     * Calculates the state derivative (velocity+acceleration of each body) of the translational motion of the system
     * at the given time and position/velocity of bodies.
     *  \param ephemerisTime Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Caartesian
     *  position/velocity of the bodies being integrated. The order of the values is defined by the order of bodies in
     *  bodiesToBeIntegratedNumerically_
     *  \return Current state derivative (velocty+acceleration) of system of bodies integrated numerically.
     */
//    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > calculateSystemStateDerivative(
//            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated ) = 0;


//    //! Returns the name of nth body from which the state is retrieved from the ephemeris.
//    /*!
//     *  Returns the name of nth body from which the state is retrieved from the ephemeris (i.e. which is involved in simulation,
//     *  but whose state is not integrated numerically).
//     *  \param bodyIndex Index from bodiesWithStateFromEphemeris_ of body of which name is to be retrieved.
//     *  \return Name of requested body.
//     */
//    std::string getNameOfEphemerisBody( int bodyIndex ) { return bodiesWithStateFromEphemeris_[ bodyIndex ]; }


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

    void setCurrentArcIndex( const int currentArcIndex )
    {
        for( accelerationMapIterator = accelerationModelsPerBody_.begin( ); accelerationMapIterator != accelerationModelsPerBody_.end( );
             accelerationMapIterator++ )
        {
            for( innerAccelerationIterator = accelerationMapIterator->second.begin( ); innerAccelerationIterator != accelerationMapIterator->second.end( );
                 innerAccelerationIterator++ )
            {
                for( unsigned int i = 0; i < innerAccelerationIterator->second.size( ); i++ )
                {
                    if( astrodynamics::acceleration_models::getAccelerationModelType( innerAccelerationIterator->second.at( i ) ) ==
                            astrodynamics::acceleration_models::power_spectrum_acceleration )
                    {
                        boost::dynamic_pointer_cast< basic_astrodynamics::PowerSpectrumAcceleration >(
                                    innerAccelerationIterator->second.at( i ) )->setCurrentArc( currentArcIndex );
                    }
                }
            }
        }
    }

    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToBeIntegratedNumerically_;
    }

    basic_astrodynamics::AccelerationMap getAccelerationsMap( )
    {
        return accelerationModelsPerBody_;
    }

    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > getCentralBodyData( )
    {
        return centralBodyData_;
    }

    int getStateSize( )
    {
        return 6 * bodiesToBeIntegratedNumerically_.size( );
    }

    TranslationalPropagatorType getPropagatorType( )
    {
        return propagatorType_;
    }

protected:

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


    basic_astrodynamics::AccelerationMap accelerationModelsPerBody_;

    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData_;

    TranslationalPropagatorType propagatorType_;

    std::vector< std::string > bodiesToBeIntegratedNumerically_;


    basic_astrodynamics::AccelerationMap::iterator accelerationMapIterator;

    std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > >::iterator
    innerAccelerationIterator;

};

}

/*
std::map< std::string, std::map< double, Eigen::VectorXd > > propagateBodies(
        AccelerationMap accelerationModelsPerBody,
        NamedBodyMap bodyMap,
        Eigen::VectorXd systemInitialState,
        double initialEphemerisTime,
        double timeStep,
        int numberOfTimeSteps,
        std::vector< std::string > bodiesToIntegrate);
*/

}
#endif // NBODYSTATEDERIVATIVE_H
