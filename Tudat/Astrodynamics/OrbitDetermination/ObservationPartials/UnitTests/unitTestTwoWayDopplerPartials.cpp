/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

Eigen::VectorXd getProperTimeRateInVectorForm(
        boost::shared_ptr< DopplerProperTimeRateInterface > properTimeRateCalculator,
        const std::vector< double >& linkEndTimes,
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const LinkEndType linkEndAssociatedWithTime )
{
    return ( Eigen::Vector1d( ) << properTimeRateCalculator->getOberverProperTimeDeviation(
                 linkEndTimes, linkEndStates ) ).finished( );
}

//! Test partial derivatives of one-way doppler observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testTwoWayDopplerPartials )
{

    using namespace tudat::gravitation;
    using namespace tudat::gravitation;
    using namespace tudat::ephemerides;
    using namespace tudat::observation_models;
    using namespace tudat::simulation_setup;
    using namespace tudat::spice_interface;
    using namespace tudat::observation_partials;
    using namespace tudat::estimatable_parameters;

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ reflector1 ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];

        for( unsigned int estimationCase  = 0; estimationCase  < 2; estimationCase ++ )
        {
            std::cout << "Case " << estimationCase << std::endl;
            // Generate one-way doppler model
            boost::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            if( estimationCase  == 0 )
            {
                twoWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            linkEnds, boost::make_shared< observation_models::ObservationSettings >(
                                observation_models::two_way_doppler,
                                boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ), bodyMap  );
            }
            else
            {
                twoWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            linkEnds, boost::make_shared< TwoWayDopplerObservationSettings >
                            (  boost::make_shared< OneWayDopplerObservationSettings >(
                                   boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                       perturbingBodies ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) ),
                               boost::make_shared< OneWayDopplerObservationSettings >(
                                   boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                       perturbingBodies ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ) ) ), bodyMap );
            }

            // Create parameter objects.
            boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
            Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
            fullEstimatableParameterSet = createEstimatableParameters( bodyMap, 1.1E7 );

            printEstimatableParameterEntries( fullEstimatableParameterSet );

            testObservationPartials< 1 >(
                        twoWayDopplerModel, bodyMap, fullEstimatableParameterSet, linkEnds, two_way_doppler, 1.0E-5,
                        true, true, 10.0, parameterPerturbationMultipliers );

        }
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ reflector1 ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];

        for( unsigned int estimationCase = 0; estimationCase  < 2; estimationCase ++ )
        {
            std::cout << "Case (with motion): " << estimationCase << std::endl;
            // Generate two-way doppler model
            boost::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            if( estimationCase  == 0 )
            {
                twoWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            linkEnds, boost::make_shared< observation_models::ObservationSettings >(
                                observation_models::two_way_doppler,
                                boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ), bodyMap  );
            }
            else
            {
                twoWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            linkEnds, boost::make_shared< TwoWayDopplerObservationSettings >
                            (  boost::make_shared< OneWayDopplerObservationSettings >(
                                   boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                       perturbingBodies ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) ),
                               boost::make_shared< OneWayDopplerObservationSettings >(
                                   boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                       perturbingBodies ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                                   boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ) ) ), bodyMap );
            }
            // Create parameter objects.
            boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
            Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
            fullEstimatableParameterSet = createEstimatableParameters( bodyMap, 1.1E7 );

            printEstimatableParameterEntries( fullEstimatableParameterSet );

            testObservationPartials< 1 >(
                        twoWayDopplerModel, bodyMap, fullEstimatableParameterSet, linkEnds, two_way_doppler, 1.0E-4, false, true,
                        1.0, parameterPerturbationMultipliers );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




