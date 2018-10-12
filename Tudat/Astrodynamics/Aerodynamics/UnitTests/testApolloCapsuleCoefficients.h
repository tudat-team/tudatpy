/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TESTAPOLLOCAPSULECOEFFICIENTS_H
#define TUDAT_TESTAPOLLOCAPSULECOEFFICIENTS_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/GeometricShapes/capsule.h"

namespace tudat
{
namespace unit_tests
{

using Eigen::Vector6d;
using mathematical_constants::PI;
using namespace aerodynamics;

std::shared_ptr< HypersonicLocalInclinationAnalysis > getApolloCoefficientInterface( )
{

    // Create test capsule.
    std::shared_ptr< geometric_shapes::Capsule > capsule
            = std::make_shared< geometric_shapes::Capsule >(
                4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

    std::vector< int > numberOfLines;
    std::vector< int > numberOfPoints;
    std::vector< bool > invertOrders;
    numberOfLines.resize( 4 );
    numberOfPoints.resize( 4 );
    invertOrders.resize( 4 );

    // Set number of analysis points.
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    numberOfLines[ 1 ] = 31;
    numberOfPoints[ 1 ] = 31;
    numberOfLines[ 2 ] = 31;
    numberOfPoints[ 2 ] = 10;
    numberOfLines[ 3 ] = 11;
    numberOfPoints[ 3 ] = 11;
    invertOrders[ 0 ] = 0;
    invertOrders[ 1 ] = 0;
    invertOrders[ 2 ] = 0;
    invertOrders[ 3 ] = 0;

    Eigen::Vector3d momentReference;
    momentReference( 0 ) = -0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;

    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 15 );

    for ( int i = 0; i < 15; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }

    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );
    std::vector< std::vector< int > > selectedMethods;
    selectedMethods.resize( 2 );
    selectedMethods[ 0 ].resize( 4 );
    selectedMethods[ 1 ].resize( 4 );

    selectedMethods[ 0 ][ 0 ] = 1;
    selectedMethods[ 0 ][ 1 ] = 5;
    selectedMethods[ 0 ][ 2 ] = 5;
    selectedMethods[ 0 ][ 3 ] = 1;
    selectedMethods[ 1 ][ 0 ] = 6;
    selectedMethods[ 1 ][ 1 ] = 3;
    selectedMethods[ 1 ][ 2 ] = 3;
    selectedMethods[ 1 ][ 3 ] = 3;

    // Create analysis object and capsule database.
    return std::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * pow( capsule->getMiddleRadius( ), 2.0 ),
                3.9116, momentReference );
}

} // namespace unit_tests

} // namespace tudat
#endif // TUDAT_TESTAPOLLOCAPSULECOEFFICIENTS_H
