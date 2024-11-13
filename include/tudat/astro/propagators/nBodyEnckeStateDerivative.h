/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NBODYENCKESTATEDERIVATIVE_H
#define TUDAT_NBODYENCKESTATEDERIVATIVE_H

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"

#include "tudat/math/root_finders/rootFinder.h"
#include "tudat/astro/gravitation/centralGravityModel.h"

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to calculate Encke's function, to be used during propagation using Encke's method
/*!
 *  Function to calculate Encke's function, to be used during propagation using Encke's method
 *  \param qValue Value of free parameter in Encke's function (typically denoted as q)
 *  \return value of Encke's function for given free parameter.
 */
template< typename StateScalarType = double >
StateScalarType calculateEnckeQFunction( const StateScalarType qValue )
{
    StateScalarType powerTerm =  mathematical_constants::getFloatingInteger< StateScalarType >( 1 )  +
            mathematical_constants::getFloatingInteger< StateScalarType >( 2 ) * qValue;
    return mathematical_constants::getFloatingInteger< StateScalarType >( 1 ) - 1.0 / ( powerTerm * std::sqrt( powerTerm ) );
}

//! Class for computing the state derivative of translational motion of N bodies, using an Encke propagator.
/*!
 * Class for computing the state derivative of translational motion of N bodies, using an Encke propagator.
 * The Encke propagator propagates the Cartesian deviation from an ideal (pre-defined) Keplerian orbit.
 * See e.g. Wakker, astro II for mathematical details.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyEnckeStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor, computes required reference quantities, and removes central gravity from acceleration list.
    /*!
     * Constructor, computes required reference quantities, and removes central gravity from acceleration list. For
     * a spherical harmonic central gravity, the C(0,0) coefficient is set to zero.
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     *  \param initialKeplerElements Kepler elements of bodiesToIntegrate, valid at initialTime.
     *  \param initialTime Time at which the initialKeplerElements provide the orbital state.
     */
    NBodyEnckeStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                               const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                               const std::vector< std::string >& bodiesToIntegrate,
                               const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements,
                               const TimeType& initialTime ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, encke, bodiesToIntegrate, true ),
        initialKeplerElements_( initialKeplerElements ),
        initialTime_( initialTime ),
        currentKeplerOrbitTime_( TUDAT_NAN )
    {
        currentKeplerianOrbitCartesianState_.resize( bodiesToIntegrate.size( ) );

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_, this->removedCentralAccelerations_ );

        // Create root-finder for Kepler orbit propagation
        rootFinder_ = root_finders::createRootFinder< StateScalarType >(
                    root_finders::newtonRaphsonRootFinderSettings(
                        TUDAT_NAN, 5.0 * std::numeric_limits< StateScalarType >::epsilon( ),
                        TUDAT_NAN, 50, root_finders::accept_result ) );


        this->createAccelerationModelList( );
    }

    //! Function to clear reference values of Encke state derivative model
    /*!
     * Function to clear reference values of Encke state derivative model, in addition to those performed in the
     * clearTranslationalStateDerivativeModel function. It resets the currentKeplerOrbitTime_ to ensure that
     * the reference orbit is recomputed.
     */
    void clearDerivedStateDerivativeModel( )
    {
        currentKeplerOrbitTime_ = TUDAT_NAN;
    }

    //! Calculates the state derivative of the translational motion of the system, using the Encke algorithm
    /*!
     *  Calculates the state derivative the translational motion of the system
     *  at the given time and state of bodies. The velocity and acceleration of each body w.r.t. their reference Kepler
     *  orbits are computed by this function, i.e. the derivative of the Encke state.
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Cartesian
     *  position/velocity deviations from the reference Kepler orbits  of the bodies being integrated.
     *  The order of the values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current derivative of Encke state (velocity + acceleration w.r.t. reference Kepler orbit) of
     *  system of bodies integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );

        // Retrieve Keplerian orbit state for each body.
        calculateKeplerTrajectoryCartesianStates( time );

        // Get Cartesian state derivative for all bodies of Encke state (excluding central gravitational accelerations).
        this->sumStateDerivativeContributions(
                    stateOfSystemToBeIntegrated, stateDerivative, true );

        // Initialize Encke algorithm variables.
        StateScalarType qValue = 0.0;
        StateScalarType qFunction = 0.0;
        StateScalarType keplerianRadius = 0.0;
        Eigen::Matrix< StateScalarType, 3, 1 > positionPerturbation = Eigen::Matrix< StateScalarType, 3, 1 >::Zero( );

        // Update state derivative for each body.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            // Get position perturbation.
            positionPerturbation = stateOfSystemToBeIntegrated.segment( i * 6, 3 );

            // Get distance from central body, assuming purely Keplerian orbit.
            keplerianRadius = currentKeplerianOrbitCartesianState_[ i ].segment( 0, 3 ).norm( );

            // Calculate Encke algorithm variables.
            qValue = positionPerturbation.dot( currentKeplerianOrbitCartesianState_[ i ].segment( 0, 3 ) +
                                               0.5 * positionPerturbation ) / ( keplerianRadius * keplerianRadius );
            qFunction = calculateEnckeQFunction( qValue );

            // Update state derivative with Encke term.
            stateDerivative.block( i * 6 + 3, 0, 3, 1 ) += static_cast< StateScalarType >(
                        centralBodyGravitationalParameters_[ i ]( ) ) /
                    ( keplerianRadius * keplerianRadius * keplerianRadius ) * (
                        ( positionPerturbation + currentKeplerianOrbitCartesianState_[ i ].segment( 0, 3 ) ) * qFunction -
                        positionPerturbation );
        }

    }

    //! Function to convert the state in the conventional form to the Encke-propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Encke propagator,
     * this transforms the Cartesian state w.r.t. the central body (conventional form) to the Cartesian deviation
     * from the Kepler orbit w.r.t. this central body (Encke form).
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid, used to computed Kepler orbits
     * \return State (outputSolution), converted to the Encke-propagator-specific form
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState = cartesianSolution;

        // Calculate Keplerian orbit states in local frames.
        calculateKeplerTrajectoryCartesianStates( time );

        // Subtract frame origin and Keplerian states from inertial state.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentState.segment( i * 6, 6 ) -= ( currentKeplerianOrbitCartesianState_[ i ] );
        }

        return currentState;

    }

    //! Function to convert the Encke-propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the Encle-propagator-specific form of the state to the conventional form. For the Encke
     * propagator, this transforms the Cartesian deviation from the Kepler orbit w.r.t. this central body (Encke form) to
     * Cartesian state w.r.t. the central body (conventional form).
     * In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in Encke-propagator-specific form (i.e. form that is used in
     * numerical integration)/
     * \param time Current time at which the state is valid
     * \param currentCartesianLocalSoluton State (internalSolution, which is Encke-formulation),
     *  converted to the 'conventional form'.
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        // Calculate Keplerian orbit state around centeal bodies.
        calculateKeplerTrajectoryCartesianStates( time );

        // Add Keplerian state to perturbation from Encke algorithm to get Cartesian state in local frames.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianLocalSoluton.block( i * 6, 0, 6, 1 ) = currentKeplerianOrbitCartesianState_[ i ] +
                    internalSolution.block( i * 6, 0, 6, 1 );
        }
    }


private:

    //! Function to calculate and set the reference Kepler orbit in Cartesian coordinates for given body.
    /*!
     * Function to calculate and set the reference Kepler orbit in Cartesian coordinates for given body.
     * \param time Time at which Kepler orbit is to be computed.
     * \param bodyIndex Index in list of bodies for which Kepler orbit is to be computed.
     */
    void calculateKeplerTrajectoryCartesianState(
            const TimeType time,
            const int bodyIndex )
    {
        // Propagate Kepler orbit to current time and set.
        currentKeplerianOrbitCartesianState_[ bodyIndex ] =
                orbital_element_conversions::convertKeplerianToCartesianElements< StateScalarType >(
                    orbital_element_conversions::propagateKeplerOrbit< StateScalarType >(
                        initialKeplerElements_.at( bodyIndex ), static_cast< StateScalarType >( time - initialTime_ ),
                        static_cast< StateScalarType >( centralBodyGravitationalParameters_.at( bodyIndex )( ) ),
                        rootFinder_ ),
                    static_cast< StateScalarType >( centralBodyGravitationalParameters_.at( bodyIndex )( ) ) );
    }

    //! Function to calculate and set the reference Kepler orbit in Cartesian coordinates for all bodies.
    /*!
     * Function to calculate and set the reference Kepler orbit in Cartesian coordinates for all bodies.
     * \param time Time at which Kepler orbits are to be computed.
     */
    void calculateKeplerTrajectoryCartesianStates(
            const TimeType time )
    {
        // Check if update is neede.
        if( !( currentKeplerOrbitTime_ == time ) )
        {
            // Iterate over bodies and calculate Cartesian state of associated Kepler orbit at current time.
            for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
            {
                calculateKeplerTrajectoryCartesianState( time, i );
            }
            currentKeplerOrbitTime_ = time;
        }
    }

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< std::function< double( ) > > centralBodyGravitationalParameters_;

    //!  Kepler elements of bodiesToIntegrate, valid at initialTime_.
    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements_;

    //! Time at which the initialKeplerElements provide the reference Keper orbit.
    TimeType initialTime_ ;

    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_/
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    //! Root finder used to propagate Kepler orbit.
    std::shared_ptr< root_finders::RootFinder< StateScalarType > > rootFinder_;

    //! Current Cartesian states of reference Kepler orbits, valid at currentKeplerOrbitTime_, computed by
    //! calculateKeplerTrajectoryCartesianStates
    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > currentKeplerianOrbitCartesianState_;

    //! Time at which the currentKeplerianOrbitCartesianState_ provide the Cartesian representation of the
    //! referfence Kepler state.
    TimeType currentKeplerOrbitTime_;


};

//extern template class NBodyEnckeStateDerivative< double, double >;


} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYENCKESTATEDERIVATIVE_H
