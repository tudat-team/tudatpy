/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <map>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/keplerPropagatorTestData.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace spice_interface;
using namespace simulation_setup;
using namespace interpolators;

BOOST_AUTO_TEST_SUITE( test_MultiArcEphemeris )

void testEphemerisDifference(
        double testEpoch,
        const std::shared_ptr< Ephemeris > ephemeris,
        const std::shared_ptr< Ephemeris > defaultEphemeris,
        bool isEpochValid,
        const bool expectSmallError )
{
    bool exceptionCaught = false;
    Eigen::Vector6d stateDifference = Eigen::Vector6d::Constant( TUDAT_NAN );

    try
    {
        stateDifference = ephemeris->getCartesianState( testEpoch ) -
                defaultEphemeris->getCartesianState( testEpoch );
    }
    catch( ... )
    {
        exceptionCaught = true;
    }
    BOOST_CHECK_EQUAL( isEpochValid, !exceptionCaught );
    if( isEpochValid )
    {
        if( expectSmallError )
        {
            BOOST_CHECK_EQUAL( stateDifference.segment( 0, 3 ).norm( ) < 1.0E-4, true );
            BOOST_CHECK_EQUAL( stateDifference.segment( 3, 3 ).norm( ) < 1.0E-10, true );
        }
        else
        {
            BOOST_CHECK_EQUAL( stateDifference.segment( 0, 3 ).norm( ) > 1.0E2, true );
            BOOST_CHECK_EQUAL( stateDifference.segment( 3, 3 ).norm( ) > 1.0E-5, true );
        }
    }
    std::cout<<stateDifference.transpose( )<<std::endl;
    std::cout<<exceptionCaught<<std::endl;

//    if( isEpochValid )
//    {
//        std::cout << stateDifference.transpose( ) << std::endl;
//    }
//    else
//    {
//        std::cout << exceptionCaught << std::endl << std::endl;
//    }
}

BOOST_AUTO_TEST_CASE( testMultiArcEphemeris )
{
   spice_interface::loadStandardSpiceKernels( );

   int numberOfArcs = 2;
   double arcLength = 1.0E6;
   double arcGap = 4.0E6;
   double interpolationStep = 1000.0;

   std::vector< double > inRangeTestPoints = { 3.0 * interpolationStep, 3.7 * interpolationStep, ( 20.0 + mathematical_constants::PI ) * interpolationStep };
   std::vector< double > outOfLagrangeRangeValidTestPoints = { 0.0, interpolationStep };
   std::vector< double > outOfLagrangeRangeInvalidTestPoints = { 0.32 * interpolationStep, 1.5 * interpolationStep };
   std::vector< double > outOfRangeTestPoints = { -137343.3, -interpolationStep };

   std::map< double, std::shared_ptr< EphemerisSettings > > perArcEphemerisSettings;
   std::shared_ptr< EphemerisSettings > defaultEphemerisSettings = directSpiceEphemerisSettings( "SSB", "ECLIPJ2000", false, false, false );

   for( int hasDefaultEphemeris = 0; hasDefaultEphemeris < 1; hasDefaultEphemeris++ )
   {
       for( int interpolationBoundaries = 0; interpolationBoundaries < 4; interpolationBoundaries++ )
       {
           std::cout<<"RUN SETTINGS "<<hasDefaultEphemeris<<" "<<interpolationBoundaries<<std::endl;;
           bool lagrangeBoundaryIsValid, outOfRangeIsValid;
           LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling;
           BoundaryInterpolationType boundaryHandling;

           if( interpolationBoundaries % 2 == 0 )
           {
               boundaryHandling = use_nan_value;
               outOfRangeIsValid = false;
           }
           else
           {
               std::cout<<"EXTRAPOLATE ******************************** "<<std::endl;
               boundaryHandling = extrapolate_at_boundary;
               outOfRangeIsValid = true;
           }

           if( interpolationBoundaries < 2 )
           {
               lagrangeBoundaryHandling = lagrange_cubic_spline_boundary_interpolation;
               lagrangeBoundaryIsValid = true;
           }
           else
           {
               lagrangeBoundaryHandling = lagrange_boundary_nan_interpolation;
               lagrangeBoundaryIsValid = false;
               outOfRangeIsValid = false;
           }

           std::vector< double > arcStartTimes;
           std::vector< double > arcEndTimes;

           for( int i = 0; i < numberOfArcs; i++ )
           {
               double arcStart = static_cast< double >( i ) * ( arcLength + arcGap );
               double arcEnd = static_cast< double >( i ) * ( arcLength + arcGap ) + arcLength;

               arcStartTimes.push_back( arcStart );
               arcEndTimes.push_back( arcEnd );

               perArcEphemerisSettings[ arcStart ] = interpolatedSpiceEphemerisSettings(
                       arcStart,
                       arcEnd,
                       1000.0,
                       "SSB",
                       "ECLIPJ2000",
                       std::make_shared< LagrangeInterpolatorSettings >(
                               6, 0, huntingAlgorithm, lagrangeBoundaryHandling, boundaryHandling ) );
           }

           std::shared_ptr< EphemerisSettings > ehemerisSettings;
           if( hasDefaultEphemeris == 0 )
           {
               ehemerisSettings = multiArcEphemerisSettings( perArcEphemerisSettings );
           }
           else
           {
               ehemerisSettings = multiArcEphemerisSettings( perArcEphemerisSettings, "SSB", "ECLIPJ2000", defaultEphemerisSettings );
           }

           std::shared_ptr< Ephemeris > ephemeris = createBodyEphemeris( ehemerisSettings, "Earth" );
           std::shared_ptr< Ephemeris > defaultEphemeris = createBodyEphemeris( defaultEphemerisSettings, "Earth" );

           double testTime;
           Eigen::Vector6d stateDifference;
           bool exceptionCaught;
           for( int i = 0; i < numberOfArcs; i++ )
           {
               for( unsigned int j = 0; j < inRangeTestPoints.size( ); j++ )
               {
                   const double deltaT = inRangeTestPoints.at( j );

                   testEphemerisDifference( arcStartTimes.at( i ) + deltaT, ephemeris, defaultEphemeris, true, true );
                   testEphemerisDifference( arcEndTimes.at( i ) - deltaT, ephemeris, defaultEphemeris, true, true );
               }

               for( unsigned int j = 0; j < outOfRangeTestPoints.size( ); j++ )
               {
                   std::cout<<"Arc "<<i<<" Point "<<j<<std::endl;
                   const double deltaT = outOfRangeTestPoints.at( j );

                   std::cout<<"Arc start "<<std::endl;
                   testEphemerisDifference( arcStartTimes.at( i ) + deltaT, ephemeris, defaultEphemeris, outOfRangeIsValid,
                                            static_cast< bool >( hasDefaultEphemeris ) );
                   std::cout<<"Arc end "<<std::endl;
                   testEphemerisDifference( arcEndTimes.at( i ) - deltaT, ephemeris, defaultEphemeris, outOfRangeIsValid,
                                            static_cast< bool >( hasDefaultEphemeris ) );

               }

               for( unsigned int j = 0; j < outOfLagrangeRangeValidTestPoints.size( ); j++ )
               {
                   std::cout<<"Arc B "<<i<<" Point "<<j<<std::endl;
                   const double deltaT = outOfLagrangeRangeValidTestPoints.at( j );

                   std::cout<<"Arc start "<<std::endl;
                   testEphemerisDifference( arcStartTimes.at( i ) + deltaT, ephemeris, defaultEphemeris, lagrangeBoundaryIsValid, true );
                   std::cout<<"Arc end "<<std::endl;
                   testEphemerisDifference( arcEndTimes.at( i ) - deltaT, ephemeris, defaultEphemeris, lagrangeBoundaryIsValid, true );

               }

               for( unsigned int j = 0; j < outOfLagrangeRangeInvalidTestPoints.size( ); j++ )
               {
                   std::cout<<"Arc C "<<i<<" Point "<<j<<std::endl;
                   const double deltaT = outOfLagrangeRangeInvalidTestPoints.at( j );

                   std::cout<<"Arc start "<<std::endl;
                   testEphemerisDifference( arcStartTimes.at( i ) + deltaT, ephemeris, defaultEphemeris, lagrangeBoundaryIsValid,
                                            static_cast< bool >( hasDefaultEphemeris ) );
                   std::cout<<"Arc end "<<std::endl;
                   testEphemerisDifference( arcEndTimes.at( i ) - deltaT, ephemeris, defaultEphemeris, lagrangeBoundaryIsValid,
                                            static_cast< bool >( hasDefaultEphemeris )  );

               }
           }



       }
   }


}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
