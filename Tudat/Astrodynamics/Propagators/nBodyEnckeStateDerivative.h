#ifndef TUDAT_NBODYENCKESTATEDERIVATIVE_H
#define TUDAT_NBODYENCKESTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"

#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType = double >
StateScalarType calculateEnckeQFunction( const StateScalarType qValue )
{
    return mathematical_constants::getFloatingInteger< StateScalarType >( 1 ) -
            std::pow( mathematical_constants::getFloatingInteger< StateScalarType >( 1 )  +
                     mathematical_constants::getFloatingInteger< StateScalarType >( 2 ) * qValue,
                      mathematical_constants::getFloatingFraction< StateScalarType >( -3, 2 ) );
}


std::vector< boost::function< double( ) > > removeCentralGravityAccelerations(
        const std::vector< std::string >& centralBodies, const std::vector< std::string >& bodiesToBeIntegratedNumerically,
        basic_astrodynamics::AccelerationMap& accelerationModelsPerBody );

template< typename StateScalarType = double, typename TimeType = double >
class NBodyEnckeStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:
    NBodyEnckeStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                               const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                               const std::vector< std::string >& bodiesToIntegrate,
                               const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements,
                               const TimeType& initialTime ):
        NBodyStateDerivative< StateScalarType, TimeType >( accelerationModelsPerBody, centralBodyData, encke, bodiesToIntegrate ),
        initialKeplerElements_( initialKeplerElements ),
        initialTime_( initialTime )
    {
        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_ );

        rootFinder_ = boost::make_shared< root_finders::NewtonRaphsonCore< StateScalarType > >(
                    boost::bind( &root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition< StateScalarType >::
                                 checkTerminationCondition,
                                 boost::make_shared< root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition< StateScalarType > >(
                                     20.0 * std::numeric_limits< StateScalarType >::epsilon( ), 1000 ),
                                 _1, _2, _3, _4, _5 ) );
    }

    virtual void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );

        // Retrieve Keplerian orbit state for each body.
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > keplerianOrbitCartesianState =
                calculateKeplerTrajectoryCartesianStates( time );;

        // Get Cartesian state derivative for all bodies of Encke state (excluding central gravitational accelerations).
        this->sumStateDerivativeContributions(
                    stateOfSystemToBeIntegrated, stateDerivative );

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
            keplerianRadius = keplerianOrbitCartesianState[ i ].segment( 0, 3 ).norm( );

            // Calculate Encke algorithm variables.
            qValue = positionPerturbation.dot( keplerianOrbitCartesianState[ i ].segment( 0, 3 ) + 0.5 * positionPerturbation ) /
                    ( keplerianRadius * keplerianRadius );
            qFunction = calculateEnckeQFunction( qValue );

            // Update state derivative with Encke term.
            stateDerivative.block( i * 6 + 3, 0, 3, 1 ) += static_cast< StateScalarType >(
                        centralBodyGravitationalParameters_[ i ]( ) ) /
                    ( keplerianRadius * keplerianRadius * keplerianRadius ) * (
                        ( positionPerturbation + keplerianOrbitCartesianState[ i ].segment( 0, 3 ) ) * qFunction - positionPerturbation );
        }

    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution, const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState = cartesianSolution;

        // Calculate Keplerian orbit states in local frames.
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > keplerianOrbitCartesianState =
                calculateKeplerTrajectoryCartesianStates( time );

        // Subtract frame origin and Keplerian states from inertial state.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentState.segment( i * 6, 6 ) -= ( keplerianOrbitCartesianState[ i ] );
        }


        return currentState;

    }

    void convertToOutputSolution(
                const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
                Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        // Calculate Keplerian orbit state around centeal bodies.
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > keplerianOrbitCartesianState =
                calculateKeplerTrajectoryCartesianStates( time );

        // Add Keplerian state to perturbation from Encke algorithm to get Cartesian state in local frames.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianLocalSoluton.segment( i * 6, 6 ) = keplerianOrbitCartesianState[ i ] +
                    internalSolution.block( i * 6, 0, 6, 1 );
        }
    }

    double getInitialTime( )
    {
        return initialTime_;
    }

    void updateInitialKeplerElements(
            const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements )
    {
        initialKeplerElements_ = initialKeplerElements;
    }

private:
    Eigen::Matrix< StateScalarType, 6, 1 >  calculateKeplerTrajectoryCartesianState(
            const TimeType time,
            const int bodyIndex )
    {
        // Propagate Kepler orbit to current time and return.
        return orbital_element_conversions::convertKeplerianToCartesianElements< StateScalarType >(
                    orbital_element_conversions::propagateKeplerOrbit< StateScalarType >(
                        initialKeplerElements_[ bodyIndex ], static_cast< StateScalarType >( time - initialTime_ ),
                        static_cast< StateScalarType >( centralBodyGravitationalParameters_[ bodyIndex ]( ) ),
                        rootFinder_ ), static_cast< StateScalarType >( centralBodyGravitationalParameters_[ bodyIndex ]( ) ) );
    }

    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > calculateKeplerTrajectoryCartesianStates(
            const TimeType time )
    {
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > keplerianOrbitCartesianState;
        keplerianOrbitCartesianState.resize( this->bodiesToBeIntegratedNumerically_.size( ) );

        // Iterate over bodies and calculate Cartesian state of associated Kepler orbit at current time.
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            keplerianOrbitCartesianState[ i ] = calculateKeplerTrajectoryCartesianState( time, i );
        }

        return keplerianOrbitCartesianState;
    }


    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters_;

    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements_;

    TimeType initialTime_ ;

    std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    boost::shared_ptr< root_finders::RootFinderCore< StateScalarType > > rootFinder_;
};

template< typename StateScalarType = double, typename TimeType = double >
void updateEnckePropagator(
        const boost::shared_ptr< NBodyEnckeStateDerivative< StateScalarType, TimeType > > enckeStateDerivative,
        const simulation_setup::NamedBodyMap& bodyMap,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newInitialBodyStates )
{
    if( newInitialBodyStates.rows( )  != enckeStateDerivative->getStateSize( ) )
    {
        std::cerr<<"Error when updating Encke state derivative, state size is inconsistent"<<std::endl;
    }

    // Get central bodies of bodies that are integrated.
    std::vector< std::string > centralBodies = enckeStateDerivative->getCentralBodyData( )->getCentralBodies( );

    // Intialize vector of initial kepler elements
    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements;
    initialKeplerElements.resize( centralBodies.size( ) );

    for( unsigned int i = 0; i < centralBodies.size( ); i++ )
    {
        // Get and verify central body.
        boost::shared_ptr< simulation_setup::Body > currentCentralBody = bodyMap.at( centralBodies[ i ] );
        if( currentCentralBody == NULL )
        {
            std::cerr<<"Error when updating encke propagator, central body "<<centralBodies[ i ]<<" is not a celestial body."<<std::endl;
        }

        // Recalculate initial Kepler elements of bodies w.r.t. their respective centraln bodies.
        initialKeplerElements[ i ] = orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                    newInitialBodyStates.segment( i * 6, 6 ), static_cast< StateScalarType >(
                        bodyMap.at( centralBodies[ i ] )->getGravityFieldModel( )->
                        getGravitationalParameter( ) ) );

    }

    // Reset initial kepler elements in propagator
    enckeStateDerivative->updateInitialKeplerElements( initialKeplerElements );
}


}

}

#endif // TUDAT_NBODYENCKESTATEDERIVATIVE_H
