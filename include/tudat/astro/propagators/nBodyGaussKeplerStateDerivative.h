/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H
#define TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H

#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the Gauss planetary equations for Kepler elements
/*!
 * Function to evaluate the Gauss planetary equations for Kepler elements, providing the time-derivatives of the
 * Kepler elements from the accelerations expressed in an RSW frame (see Vallado, 2001). This function takes a number
 * of precomputed quantities as input, to reduce computational burden
 * \param currentOsculatingKeplerElements Current osculating Kepler elements of the body for which the Gauss equations are
 * to be evaluated
 * \param accelerationsInRswFrame Accelerations acting on body, expressed in RSW frame
 * \param semiLatusRectum Semi-latus rectum of osculating orbit
 * \param distance Distance from propagated body to its central body
 * \param meanMotion Instantaneous mean motion of body w.r.t. its central body
 * \param orbitalAngularMomentum Norm of orbital angular momentum vector of body, w.r.t. its central body
 * \return Time derivatives of osculating Kepler elements.
 */
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double semiLatusRectum,
        const double distance,
        const double meanMotion,
        const double orbitalAngularMomentum );

//! Function to evaluate the Gauss planetary equations for Kepler elements
/*!
 * Function to evaluate the Gauss planetary equations for Kepler elements, providing the time-derivatives of the
 * Kepler elements from the accelerations expressed in an RSW frame (see Vallado, 2001).
 * \param currentOsculatingKeplerElements Current osculating Kepler elements of the body for which the Gauss equations are
 * to be evaluated
 * \param accelerationsInRswFrame Accelerations acting on body, expressed in RSW frame
 * \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
 * \return Time derivatives of osculating Kepler elements.
 */
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter );

//! Function to evaluate the Gauss planetary equations for Kepler elements
/*!
 * Function to evaluate the Gauss planetary equations for Kepler elements, providing the time-derivatives of the
 * Kepler elements from the accelerations expressed in an RSW frame (see Vallado, 2001).  This function takes the accelerations
 * in the inertial frame, as well as the Cartesian inertial state, and converts the accelerations to the RSW frame.
 * \param currentOsculatingKeplerElements Current osculating Kepler elements of the body for which the Gauss equations are
 * to be evaluated
 * \param currentCartesianState Current Cartesian state of the body for which the Gauss equations are to be evaluated
 * \param accelerationsInInertialFrame Accelerations acting on body, expressed in inertial frame
 * \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
 * \return Time derivatives of osculating Kepler elements.
 */
Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter );

//! Class for computing the state derivative of translational motion of N bodies, using a Gauss method with Kepler elements.
/*!
 * Class for computing the state derivative of translational motion of N bodies, using a Gauss method with Kepler elements.
 * In this method, the derivative of the Kepler elements are computed from the total Cartesian accelerations, with the Kepler
 * elements of the bodies the states being numerically propagated.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyGaussKeplerStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on. This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyGaussKeplerStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                                     const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                                     const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, gauss_keplerian, bodiesToIntegrate, true )
    {
        currentTrueAnomalies_.resize( bodiesToIntegrate.size( ) );

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_, this->removedCentralAccelerations_ );
        this->createAccelerationModelList( );

    }

    //! Destructor
    ~NBodyGaussKeplerStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system, using the Gauss equations for Kepler elememts
    /*!
     *  Calculates the state derivative of the translational motion of the system, using the Gauss equations for
     *  Kepler elements. The input is the current state in Kepler elememts. The state derivate of this
     *  set is computed. TO do so the accelerations are internally transformed into the RSW frame, using the current
     *  Cartesian state as set by the last call to the convertToOutputSolution function
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Kepler
     *  elements of the bodies being integrated.
     *  The order of the values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current derivative of the Kepler elements of the
     *  system of bodies integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        // Get total inertial accelerations acting on bodies
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, false );

        // Compute RSW accelerations for each body, and evaluate Gauss equations for Keplerian elements.
        Eigen::Vector3d currentAccelerationInRswFrame;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentAccelerationInRswFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                        currentCartesianLocalSolution_.segment( i * 6, 6 ) ) *
                    stateDerivative.block( i * 6 + 3, 0, 3, 1 ).template cast< double >( );

            stateDerivative.block( i * 6, 0, 6, 1 ) = computeGaussPlanetaryEquationsForKeplerElements(
                        ( Eigen::Vector6d( ) << stateOfSystemToBeIntegrated.block( i * 6, 0, 5, 1 ).template cast< double >( ),
                          currentTrueAnomalies_.at( i ) ).finished( ), currentAccelerationInRswFrame,
                        centralBodyGravitationalParameters_.at( i )( ) ).template cast< StateScalarType >( );
        }
    }

    //! Function to convert the state in the conventional form to the Kepler Elements form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Gauss-Kepler propagator,
     * this transforms the Cartesian state w.r.t. the central body (conventional form) to the Kepler Elements
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid, used to computed Kepler orbits
     * \return State (outputSolution), converted to the Kepler Elements
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( cartesianSolution.rows( ) );

        // Subtract frame origin and Keplerian states from inertial state.
        Eigen::Matrix< StateScalarType, 6, 1 > currentCartesianState;
        Eigen::Matrix< StateScalarType, 6, 1 > currentKeplerianState;

        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianState = cartesianSolution.block( i * 6, 0, 6, 1 );
            currentKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements (
                        currentCartesianState, static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ) );
            StateScalarType eccentricAnomaly = orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                        currentKeplerianState( 5 ), currentKeplerianState( 1 ) );
            StateScalarType meanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomaly, currentKeplerianState( 1 ) );
            currentKeplerianState( 5 ) = meanAnomaly;
            currentState.segment( i * 6, 6 ) = currentKeplerianState;
        }

        return currentState;

    }

    //! Function to convert the Kepler states of the bodies to the conventional form.
    /*!
     * Function to convert the Kepler Elements state to the conventional form. For the Gauss-Kepler
     * propagator, this transforms Kepler Elements w.r.t. the central bodies to the Cartesian states w.r.t. these
     * same central bodies: In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in Kepler Elemements (i.e. form that is used in
     * numerical integration)
     * \param time Current time at which the state is valid
     * \param currentCartesianLocalSolution State (internalSolution, which is Encke-formulation),
     *  converted to the 'conventional form' (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution )
    {
        // Add Keplerian state to perturbation from Encke algorithm to get Cartesian state in local frames.
        Eigen::Matrix< StateScalarType, 6, 1 > currentKeplerianState;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            try
            {
                currentKeplerianState = internalSolution.block( i * 6, 0, 6, 1 );
                StateScalarType currentTrueAnomaly;
                if( currentKeplerianState( 1 ) < 1.0 )
                {
                    StateScalarType currentEccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(
                                currentKeplerianState( 1 ), currentKeplerianState( 5 ) );
                    currentTrueAnomaly = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                                currentEccentricAnomaly, currentKeplerianState( 1 ) );
                    currentKeplerianState( 5 ) = currentTrueAnomaly;
                }
                else
                {
                    StateScalarType currentEccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToHyperbolicEccentricAnomaly(
                                currentKeplerianState( 1 ), currentKeplerianState( 5 ) );
                    currentTrueAnomaly = orbital_element_conversions::convertHyperbolicEccentricAnomalyToTrueAnomaly(
                                currentEccentricAnomaly, currentKeplerianState( 1 ) );
                    currentKeplerianState( 5 ) = currentTrueAnomaly;
                }

                currentTrueAnomalies_[ i ] = currentTrueAnomaly;
                currentCartesianLocalSolution.block( i * 6, 0, 6, 1 ) = orbital_element_conversions::convertKeplerianToCartesianElements(
                            currentKeplerianState, static_cast< StateScalarType >( centralBodyGravitationalParameters_.at( i )( ) ) );
            }
            catch( std::runtime_error const& )
            {
                currentTrueAnomalies_[ i ] = TUDAT_NAN;
                currentCartesianLocalSolution.block( i * 6, 0, 6, 1 ) =
                        Eigen::Matrix< StateScalarType, 6, 1 >::Constant( TUDAT_NAN );
            }

        }

        currentCartesianLocalSolution_ = currentCartesianLocalSolution.template cast< double >( );
    }

private:

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< std::function< double( ) > > centralBodyGravitationalParameters_;

    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    //! Current full Cartesian state of the propagated bodies, w.r.t. trhe central bodies
    /*!
     *  Current full Cartesian state of the propagated bodies, w.r.t. trhe central bodies. These variables are set when calling
     *  the convertToOutputSolution function.
     */
    Eigen::VectorXd currentCartesianLocalSolution_;

    //! List of true anomalies of the bodies being propagated, computed by the last call to the convertToOutputSolution function
    std::vector< double > currentTrueAnomalies_;

};

extern template class NBodyGaussKeplerStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyGaussKeplerStateDerivative< long double, double >;
extern template class NBodyGaussKeplerStateDerivative< double, Time >;
extern template class NBodyGaussKeplerStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H
