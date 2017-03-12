/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;


BOOST_AUTO_TEST_SUITE( test_one_way_doppler_model )


BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
{
    // Load Spice kernels
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
    linkEnds[ receiver ] = std::make_pair( "Mars" , ""  );

    // Create observation settings
    boost::shared_ptr< ObservationSettings > observableSettings = boost::make_shared< ObservationSettings >
            ( one_way_doppler );

    // Create observation model.
    boost::shared_ptr< ObservationModel< 1, double, double> > observationModel =
           ObservationModelCreator< 1, double, double>::createObservationModel(
                linkEnds, observableSettings, bodyMap );

    boost::shared_ptr< OneWayDopplerObservationModel< double, double> > dopplerObservationModel =
            boost::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >( observationModel );

    // Test observable for both fixed link ends
    for( unsigned testCase = 0; testCase < 2; testCase++ )
    {

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;

        // Define link end
        LinkEndType referenceLinkEnd;
        if( testCase == 0 )
        {
            referenceLinkEnd = transmitter;
        }
        else
        {
            referenceLinkEnd = receiver;
        }

        // Compute observable
        double dopplerObservable = observationModel->computeObservationsWithLinkEndData(
                    observationTime, referenceLinkEnd, linkEndTimes, linkEndStates )( 0 );

        // Creare independent light time calculator object
        boost::shared_ptr< LightTimeCalculator< double, double > > lightTimeCalculator =
                createLightTimeCalculator( linkEnds[ transmitter ], linkEnds[ receiver ], bodyMap );
        Eigen::Vector6d transmitterState, receiverState;        
        // Compute light time
        double lightTime = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, observationTime, testCase );

        // Compare light time calculator link end conditions with observation model
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( receiverState, linkEndStates.at( 1 ), 1.0E-15 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transmitterState, linkEndStates.at( 0 ), 1.0E-15 );

            if( testCase == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime  - linkEndTimes.at( 0 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( observationTime + lightTime - linkEndTimes.at( 1 ) ), 1.0E-15 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime - linkEndTimes.at( 1 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( observationTime - lightTime - linkEndTimes.at( 0 ) ), 1.0E-15 );
            }
        }

        // Compute numerical partial derivative of light time.
        double timePerturbation = 100.0;
        double upPerturbedLightTime = lightTimeCalculator->calculateLightTime( linkEndTimes.at( 0 ) + timePerturbation, false );
        double downPerturbedLightTime = lightTimeCalculator->calculateLightTime( linkEndTimes.at( 0 ) - timePerturbation, false );

        double lightTimeSensitivity = ( upPerturbedLightTime - downPerturbedLightTime ) / ( 2.0 * timePerturbation );

        // Test numerical derivative against Doppler observable
        BOOST_CHECK_SMALL( std::fabs( lightTimeSensitivity  - dopplerObservable ), 1.0E-14 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


