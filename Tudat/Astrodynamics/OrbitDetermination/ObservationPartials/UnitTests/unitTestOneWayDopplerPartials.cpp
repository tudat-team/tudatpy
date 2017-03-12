/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayDopplerObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_one_way_observation_partials)

Eigen::Vector3d computeUnitVectorToReceiverFromTransmitterState(
        const Eigen::Vector3d receiverPosition,
        const boost::function< Eigen::Vector6d( const double ) > transmitterStateFunction,
        const double evaluationTime )
{
    return ( receiverPosition - transmitterStateFunction( evaluationTime ).segment( 0, 3 ) ).normalized( );
}


Eigen::Vector3d computeUnitVectorToReceiverFromReceiverState(
        const boost::function< Eigen::Vector6d( const double ) > receiverStateFunction,
        const Eigen::Vector3d transmitterPosition,
        const double evaluationTime )
{
    return ( receiverStateFunction( evaluationTime ).segment( 0, 3 ) - transmitterPosition ).normalized( );
}

//! Test partial derivatives of one-way doppler observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testOneWayDopplerPartials )
{

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    std::cout<<"*******************************************************"<<std::endl;

    // Test ancilliary functions
    {
        double nominalEvaluationTime = 1.1E7;

        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Create transmitter/receriver state functions
        boost::function< Eigen::Vector6d( const double ) > transmitterStateFunction =
                getLinkEndCompleteEphemerisFunction< double, double >( linkEnds[ transmitter ], bodyMap );
        boost::function< Eigen::Vector6d( const double ) > receiverStateFunction =
                getLinkEndCompleteEphemerisFunction< double, double >( linkEnds[ receiver ], bodyMap );

        // Define (independent!) transmission/reception times
        double transmissionTime = nominalEvaluationTime;
        double receptionTime = nominalEvaluationTime + 1.0E3;

        // Compute associated states
         Eigen::Vector6d nominalTransmitterState = transmitterStateFunction( transmissionTime );
         Eigen::Vector6d nominalReceiverState = receiverStateFunction( receptionTime );
         Eigen::Vector3d nominalVectorToReceiver = ( nominalReceiverState - nominalTransmitterState ).segment( 0, 3 );

         double timePerturbation = 100.0;

         // Partials for fixed receiver
         {
             // Compute numerical derivative of transmitter state for acceleration)
             Eigen::Vector6d numericalStateDerivative = numerical_derivatives::computeCentralDifference(
                         transmitterStateFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute unit vector derivative numerically
             boost::function< Eigen::Vector3d( const double ) > unitVectorFunction =
                     boost::bind( &computeUnitVectorToReceiverFromTransmitterState,
                                  nominalReceiverState.segment( 0, 3 ), transmitterStateFunction, _1 );
             Eigen::Vector3d numericalUnitVectorDerivative = numerical_derivatives::computeCentralDifference(
                         unitVectorFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute projected velocoty vector derivative numerically
             boost::function< double( const double) > projectedVelocityFunction =
                     boost::bind( &calculateLineOfSightVelocityAsCFractionFromTransmitterStateFunction< double, double >,
                                  nominalReceiverState.segment( 0, 3 ), transmitterStateFunction, _1 );
             double numericalProjectedVelocityDerivative =
                     numerical_derivatives::computeCentralDifference(
                         projectedVelocityFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute analytical partial derivatives
             Eigen::Vector3d analyticalUnitVectorDerivative =
                     -computePartialOfUnitVectorWrtLinkEndTime(
                         nominalVectorToReceiver, nominalVectorToReceiver.normalized( ),
                         nominalVectorToReceiver.norm( ), nominalTransmitterState.segment( 3, 3 ) );
             double analyticalProjectedVelocityDerivative = computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                         nominalVectorToReceiver, nominalTransmitterState.segment( 3, 3 ), \
                         numericalStateDerivative.segment( 3, 3 ), false );


             for( unsigned int i = 0; i < 3; i++ )
             {
                 BOOST_CHECK_SMALL( std::fabs( analyticalUnitVectorDerivative( i ) - numericalUnitVectorDerivative( i ) ), 1.0E-18 );

             }
             BOOST_CHECK_SMALL( std::fabs( analyticalProjectedVelocityDerivative / physical_constants::SPEED_OF_LIGHT -
                                           numericalProjectedVelocityDerivative ), 1.0E-22 );
         }


         // Partials for fixed transmitter
         {
             // Compute numerical derivative of receiver state for acceleration)
             Eigen::Vector6d numericalStateDerivative = numerical_derivatives::computeCentralDifference(
                         receiverStateFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute unit vector derivative numerically
             boost::function< Eigen::Vector3d( const double ) > unitVectorFunction =
                     boost::bind( &computeUnitVectorToReceiverFromReceiverState,
                                  receiverStateFunction, nominalTransmitterState.segment( 0, 3 ), _1 );
             Eigen::Vector3d numericalUnitVectorDerivative = numerical_derivatives::computeCentralDifference(
                         unitVectorFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute projected velocoty vector derivative numerically
             boost::function< double( const double) > projectedVelocityFunction =
                     boost::bind( &calculateLineOfSightVelocityAsCFractionFromReceiverStateFunction< double, double >,
                                  receiverStateFunction, nominalTransmitterState.segment( 0, 3 ), _1 );
             double numericalProjectedVelocityDerivative =
                     numerical_derivatives::computeCentralDifference(
                         projectedVelocityFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

             // Compute analytical partial derivatives
             Eigen::Vector3d analyticalUnitVectorDerivative =
                     computePartialOfUnitVectorWrtLinkEndTime(
                         nominalVectorToReceiver, nominalVectorToReceiver.normalized( ),
                         nominalVectorToReceiver.norm( ), nominalReceiverState.segment( 3, 3 ) );
             double analyticalProjectedVelocityDerivative = computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                         nominalVectorToReceiver, nominalReceiverState.segment( 3, 3 ), \
                         numericalStateDerivative.segment( 3, 3 ), true );

             for( unsigned int i = 0; i < 3; i++ )
             {
                 BOOST_CHECK_SMALL( std::fabs( analyticalUnitVectorDerivative( i ) - numericalUnitVectorDerivative( i ) ), 1.0E-18 );

             }
             BOOST_CHECK_SMALL( std::fabs( analyticalProjectedVelocityDerivative / physical_constants::SPEED_OF_LIGHT -
                                           numericalProjectedVelocityDerivative ), 1.0E-22 );
         }

    }

    std::cout<<"*******************************************************"<<std::endl;
    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way doppler model
        boost::shared_ptr< ObservationModel< 1 > > oneWayDopplerModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    linkEnds, boost::make_shared< observation_models::ObservationSettings >(
                        observation_models::one_way_doppler ), bodyMap  );

        // Create parameter objects.
        boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodyMap, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayDopplerModel, bodyMap, fullEstimatableParameterSet, linkEnds, one_way_doppler, 1.0E-4, true, true );
    }

    std::cout<<"*******************************************************"<<std::endl;

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way doppler model
        boost::shared_ptr< ObservationModel< 1 > > oneWayDopplerModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    linkEnds, boost::make_shared< observation_models::ObservationSettings >(
                        observation_models::one_way_doppler ), bodyMap  );

        // Create parameter objects.
        boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodyMap, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayDopplerModel, bodyMap, fullEstimatableParameterSet, linkEnds, one_way_doppler, 1.0E-4, false, true );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




