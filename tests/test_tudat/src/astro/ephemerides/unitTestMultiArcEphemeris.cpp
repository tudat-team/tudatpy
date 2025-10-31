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

BOOST_AUTO_TEST_CASE( testMultiArcEphemeris )
{
   spice_interface::loadStandardSpiceKernels( );

   int numberOfArcs = 10;
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
                               6, 0, huntingAlgorithm, lagrange_cubic_spline_boundary_interpolation, extrapolate_at_boundary ) );
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
                   testTime = arcStartTimes.at( i ) + inRangeTestPoints.at( j );
                   stateDifference = ephemeris->getCartesianState( testTime ) -
                           defaultEphemeris->getCartesianState( testTime );
                   std::cout<<stateDifference.transpose( )<<std::endl;

                   testTime = arcEndTimes.at( i ) - inRangeTestPoints.at( j );
                   stateDifference = ephemeris->getCartesianState( testTime ) -
                           defaultEphemeris->getCartesianState( testTime );

                   std::cout<<stateDifference.transpose( )<<std::endl<<std::endl;
               }

               for( unsigned int j = 0; j < outOfLagrangeRangeValidTestPoints.size( ); j++ )
               {
                   testTime = arcStartTimes.at( i ) + outOfLagrangeRangeValidTestPoints.at( j );

                   exceptionCaught = false;
                   try
                   {
                       stateDifference = ephemeris->getCartesianState( testTime ) -
                               defaultEphemeris->getCartesianState( testTime );
                   }
                   catch( ... )
                   {
                       exceptionCaught = true;
                   }
                   if( lagrangeBoundaryIsValid )
                   {
                       std::cout<<stateDifference.transpose( )<<std::endl;
                   }
                   else
                   {
                       std::cout<<exceptionCaught<<std::endl<<std::endl;
                   }

                   testTime = arcEndTimes.at( i ) - outOfLagrangeRangeValidTestPoints.at( j );

                   exceptionCaught = false;
                   try
                   {
                       stateDifference = ephemeris->getCartesianState( testTime ) -
                               defaultEphemeris->getCartesianState( testTime );
                   }
                   catch( ... )
                   {
                       exceptionCaught = true;
                   }
                   if( lagrangeBoundaryIsValid )
                   {
                       std::cout<<stateDifference.transpose( )<<std::endl;
                   }
                   else
                   {
                       std::cout<<exceptionCaught<<std::endl<<std::endl;
                   }

               }

           }

       }
   }


}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
