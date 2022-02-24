/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H
#define TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H

#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the Gauss planetary equations for modified equinictial elements
/*!
 * Function to evaluate the Gauss planetary equations for modified equinictial elements, providing the time-derivatives of the
 * modified equinictial elements from the accelerations expressed in an RSW frame (see Vallado, 2001).
 * \param osculatingModifiedEquinoctialElements Current osculating modified equinoctial elements of the body for which the Gauss
 * equations are  to be evaluated
 * \param accelerationsInRswFrame Accelerations acting on body, expressed in RSW frame
 * \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
 * \return Time derivatives of osculating dified equinictial elements.
 */
Eigen::Vector6d computeGaussPlanetaryEquationsForModifiedEquinoctialElements(
        const Eigen::Vector6d& osculatingModifiedEquinoctialElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter );

//! Class for computing the state derivative of translational motion of N bodies, using a Gauss method with MEE.
/*!
 * Class for computing the state derivative of translational motion of N bodies, using a Gauss method with Modified Equinoctial
 * elements (MEE). In this method, the derivative of the MEE are computed from the total Cartesian accelerations, with the MEE
 * of the bodies the states being numerically propagated.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyGaussModifiedEquinictialStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     *  \param initialKeplerElements Kepler elements of bodiesToIntegrate, at initial propagation time
     */
    NBodyGaussModifiedEquinictialStateDerivative(
            const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
            const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const std::vector< std::string >& bodiesToIntegrate ,
            const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, gauss_modified_equinoctial, bodiesToIntegrate, true )
    {
        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_, this->removedCentralAccelerations_ );
        this->createAccelerationModelList( );

        // Check if singularity in equations of motion should be palced at 0 or 180 degrees inclination
        flipSingularities_.resize( bodiesToIntegrate.size( ) );
        for( unsigned int i = 0; i < initialKeplerElements.size( ); i++ )
        {
            if( initialKeplerElements.at( i )( orbital_element_conversions::inclinationIndex ) >  mathematical_constants::PI / 2.0 )
            {
                flipSingularities_[ i ] = false;
                std::cerr << "Warning when using Gauss-MEE propagation, body " << bodiesToIntegrate.at( i ) << " has inclination of "
                          << initialKeplerElements.at( i )( orbital_element_conversions::inclinationIndex ) * 180.0 /
                           mathematical_constants::PI << " degrees, but propagator has singularity at i=180 degrees" << std::endl;
            }
            else
            {
                flipSingularities_[ i ] = false;
            }

            if( initialKeplerElements.at( i )( orbital_element_conversions::inclinationIndex ) >
                    mathematical_constants::PI  )
            {
                throw std::runtime_error( "Error when setting N Body Gauss MEE propagator, initial inclination of body is larger than pi." );
            }
        }

    }

    //! Destructor
    ~NBodyGaussModifiedEquinictialStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system, using the Gauss equations for MEE
    /*!
     *  Calculates the state derivative of the translational motion of the system, using the Gauss equations for modified
     *  equinoctial elements. The input is the current state in modified equinoctial elememts. The state derivate of this
     *  set is computed. TO do so the accelerations are internally transformed into the RSW frame, using the current
     *  Cartesian state as set by the last call to the convertToOutputSolution function
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing modified equinoctial
     *  elements of the bodies being integrated.
     *  The order of the values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current derivative of the modified equinoctial elements of the
     *  system of bodies integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        // Get total inertial accelerations acting on bodies
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, false );

        // Compute RSW accelerations for each body, and evaluate MEE Gauss equations.
        Eigen::Vector3d currentAccelerationInRswFrame;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentAccelerationInRswFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                        currentCartesianLocalSolution_.segment( i * 6, 6 ).template cast< double >( ) ) *
                    stateDerivative.block( i * 6 + 3, 0, 3, 1 ).template cast< double >( );

            stateDerivative.block( i * 6, 0, 6, 1 ) = computeGaussPlanetaryEquationsForModifiedEquinoctialElements(
                        stateOfSystemToBeIntegrated.block( i * 6, 0, 6, 1 ).template cast< double >( ), currentAccelerationInRswFrame,
                        centralBodyGravitationalParameters_.at( i )( ) ).template cast< StateScalarType >( );
        }

    }

    //! Function to convert the state in the conventional form to the Modified Equinoctial Elements form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Gauss-MEE propagator,
     * this transforms the Cartesian state w.r.t. the central body (conventional form) to the Modified Equinoctial Elements
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid, used to computed Kepler orbits
     * \return State (outputSolution), converted to the Modified Equinoctial Elements
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( cartesianSolution.rows( ) );

        // Convert state to MEE for each body
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentState.segment( i * 6, 6 ) =
                    orbital_element_conversions::convertCartesianToModifiedEquinoctialElements< StateScalarType >(
                        cartesianSolution.block( i * 6, 0, 6, 1 ), static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ), flipSingularities_.at( i ) );
        }

        return currentState;

    }

    //! Function to convert the MEE states of the bodies to the conventional form.
    /*!
     * Function to convert the Modified Equinoctial Elements state to the conventional form. For the Gauss-MEE
     * propagator, this transforms Modfied Equinoctial Elements w.r.t. the central bodies to the Cartesian states w.r.t. these
     * same central bodies: In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in Modified Equinoctial Elemements (i.e. form that is used in
     * numerical integration)
     * \param time Current time at which the state is valid
     * \param currentCartesianLocalSolution State (internalSolution, which is Encke-formulation),
     *  converted to the 'conventional form' (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution )
    {
        // Convert state to Cartesian for each body
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianLocalSolution.block( i * 6, 0, 6, 1 ) =
                    orbital_element_conversions::convertModifiedEquinoctialToCartesianElements< StateScalarType >(
                        internalSolution.block( i * 6, 0, 6, 1 ), static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ), flipSingularities_.at( i ) );
        }

        currentCartesianLocalSolution_ = currentCartesianLocalSolution;
    }


private:

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< std::function< double( ) > > centralBodyGravitationalParameters_;

    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    //! Current full Cartesian state of the propagated bodies, w.r.t. the central bodies
    /*!
     *  Current full Cartesian state of the propagated bodies, w.r.t. the central bodies. These variables are set when calling
     *  the convertToOutputSolution function.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentCartesianLocalSolution_;

    //! List of booleans to denote, per propagated body, if the singularity in the MEE equations it o be flipped to 0 inclination
    //! (if true) or if it remains at 180 degree inclination (if false)
    std::vector< bool > flipSingularities_;

};

extern template class NBodyGaussModifiedEquinictialStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class NBodyGaussModifiedEquinictialStateDerivative< long double, double >;
extern template class NBodyGaussModifiedEquinictialStateDerivative< double, Time >;
extern template class NBodyGaussModifiedEquinictialStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H
