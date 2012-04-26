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
 *      110501    L. van der Ham    First creation of code.
 *      110607    L. van der Ham    Make code compatible with Tudat revision 114.
 *      110629    L. van der Ham    Modifications according to comments first code check.
 *      110710    K. Kumar          Restructured code; added subtests.
 *      111024    K. Kumar          Error spotted in L1/L2 tests; locations seem swapped.
 *                                  Tests commented out; needs to be fixed.
 *      111027    K. Kumar          Uncommented out tests as bugs fixed by L. van der Ham.
 *      120307    K. Kumar          Moved file.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References:
 */

// Temporary notes (move to class/function doxygen):
// Problem", http://www.math.utexas.edu/users/jjames/celestMech, 2006.
// 
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 
// Reference values for position Lagrange libration points are taken from
// (James, 2006).
// 

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace crtbp = tudat::astrodynamics::gravitation::circular_restricted_three_body_problem;
using std::cerr;
using std::endl;
using std::fabs;
using tudat::NewtonRaphson;
using tudat::bodies::Planet;

//! Test determination of L1 location.
bool testL1LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter )
{
    // Declare L1 Libration Point object.
    crtbp::LibrationPoint librationPointL1;

    // Set mass parameter.
    librationPointL1.setMassParameter( massParameter );

    // Set Newton-Raphson method.
    librationPointL1.setNewtonRaphsonMethod( boost::make_shared< NewtonRaphson >( ) );

    // Compute location of Lagrange libration point.
    librationPointL1.computeLocationOfLibrationPoint( crtbp::LibrationPoint::L1 );

    // Determine location of libration point in Earth-Moon system:
    Eigen::Vector3d positionOflibrationPointL1 = librationPointL1
            .getLocationOfLagrangeLibrationPoint( );

    // Test determination of location of L1 and output cerr statements if the test fails.
    if ( fabs( positionOflibrationPointL1.x( ) - 0.83629259089993 )
         / 0.83629259089993 > 1.0e-14
         || fabs( positionOflibrationPointL1.y( ) ) > 0.0
         || fabs( positionOflibrationPointL1.z( ) ) > 0.0 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The determination of the libration point L1 does not function well" << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLibrationPointComputationErroneous;
}


//! Test determination of L2 location.
bool testL2LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter )
{
    // Declare L2 Libration Point object.
    crtbp::LibrationPoint librationPointL2;

    // Set mass parameter.
    librationPointL2.setMassParameter( massParameter );

    // Set Newton-Raphson method.
    librationPointL2.setNewtonRaphsonMethod( boost::make_shared< NewtonRaphson >( ) );

    // Compute location of Lagrange libration point.
    librationPointL2.computeLocationOfLibrationPoint( crtbp::LibrationPoint::L2 );

    // Determine location of libration point in Earth-Moon system:
    Eigen::Vector3d positionOflibrationPointL2 = librationPointL2
            .getLocationOfLagrangeLibrationPoint( );

    // Test determination of location of L2 and output cerr statements if the test fails.
    if ( fabs( positionOflibrationPointL2.x( ) - 1.15616816590553 )
         / 1.15616816590553 > 1.0e-14
         || fabs( positionOflibrationPointL2.y( ) ) > 0.0
         || fabs( positionOflibrationPointL2.z( ) ) > 0.0 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The determination of the libration point L2 does not function well" << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLibrationPointComputationErroneous;
}

//! Test determination of L3 location.
bool testL3LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter )
{
    // Declare L3 Libration Point object.
    crtbp::LibrationPoint librationPointL3;

    // Set mass parameter.
    librationPointL3.setMassParameter( massParameter );

    // Set Newton-Raphson method.
    librationPointL3.setNewtonRaphsonMethod( boost::make_shared< NewtonRaphson >( ) );

    // Compute location of Lagrange libration point.
    librationPointL3.computeLocationOfLibrationPoint( crtbp::LibrationPoint::L3 );

    // Determine location of libration point in Earth-Moon system:
    Eigen::Vector3d positionOflibrationPointL3 = librationPointL3
            .getLocationOfLagrangeLibrationPoint( );

    // Test determination of location of L3 and output cerr statements if the test fails.
    if ( fabs( positionOflibrationPointL3.x( ) + 1.00511551160689 )
         / -1.00511551160689 > 1.0e-15
         || fabs( positionOflibrationPointL3.y( ) ) > 0.0
         || fabs( positionOflibrationPointL3.z( ) ) > 0.0 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The determination of the libration point L3 does not function well" << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLibrationPointComputationErroneous;
}

//! Test determination of L4 location.
bool testL4LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter )
{
    // Declare L4 Libration Point object.
    crtbp::LibrationPoint librationPointL4;

    // Set mass parameter.
    librationPointL4.setMassParameter( massParameter );

    // Compute location of Lagrange libration point.
    librationPointL4.computeLocationOfLibrationPoint( crtbp::LibrationPoint::L4 );

    // Determine location of libration point in Earth-Moon system:
    Eigen::Vector3d positionOflibrationPointL4 = librationPointL4
            .getLocationOfLagrangeLibrationPoint( );

    // Test determination of location of L4 and output cerr statements if the test fails.
    if ( fabs( positionOflibrationPointL4.x( ) - 0.487722529 )
         / 0.487722529 > std::numeric_limits< double >::epsilon( )
         || fabs( positionOflibrationPointL4.y( ) - 0.86602540378444 )
         / 0.86602540378444 > 1.0e-14
         || fabs( positionOflibrationPointL4.z( ) ) > 0.0 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The determination of the libration point L4 does not function well" << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLibrationPointComputationErroneous;
}

//! Test determination of L5 location.
bool testL5LibrationPointLocation( bool isLibrationPointComputationErroneous,
                                   double massParameter )
{
    // Declare L5 Libration Point object.
    crtbp::LibrationPoint librationPointL5;

    // Set mass parameter.
    librationPointL5.setMassParameter( massParameter );

    // Compute location of Lagrange libration point.
    librationPointL5.computeLocationOfLibrationPoint( crtbp::LibrationPoint::L5 );

    // Determine location of libration point in Earth-Moon system:
    Eigen::Vector3d positionOflibrationPointL5 = librationPointL5
            .getLocationOfLagrangeLibrationPoint( );

    // Test determination of location of L5 and output cerr statements if the test fails.
    if ( fabs( positionOflibrationPointL5.x( ) - 0.487722529 )
         / 0.487722529 > 1.0e-15
         || fabs( positionOflibrationPointL5.y( ) + 0.86602540378444 )
         / -0.86602540378444 > 1.0e-15
         || fabs( positionOflibrationPointL5.z( ) ) > 0.0 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The determination of the libration point L5 does not function well" << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLibrationPointComputationErroneous;
}

//! Test determination of libration point locations.
int main( )
{
    // Summary of tests.
    // Test 1: Test dimensionless mass parameter computation. Benchmark data is obtained through
    //         manually calculation.
    // Test 2: Compute location of L1 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    // Test 3: Compute location of L2 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    // Test 4: Compute location of L3 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    // Test 5: Compute location of L4 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    // Test 6: Compute location of L5 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).

    // Test result initialised to false.
    bool isLibrationPointComputationErroneous = false;

    // Test 1: Test dimensionless mass parameter computation.
    // Declare Libration Point object.
    crtbp::LibrationPoint librationPoint;

    // Create predefined Earth object.
    boost::shared_ptr< Planet > predefinedEarth = boost::make_shared< Planet >( );
    predefinedEarth->setPredefinedPlanetSettings( Planet::earth );

    // Create predefined Moon object.
    boost::shared_ptr< Planet > predefinedMoon = boost::make_shared< Planet >( );
    predefinedMoon->setPredefinedPlanetSettings( Planet::moon );

    // Set bodies.
    librationPoint.setPrimaryCelestialBody( predefinedEarth );
    librationPoint.setSecondaryCelestialBody( predefinedMoon );

    // Compute mass parameter.
    librationPoint.computeMassParameter( );

    // Check if computed Earth-Moon mass parameter is too large and output cerr statements.
    if ( fabs( ( librationPoint.getMassParameter( ) ) - 0.01215295290792761 )
         / 0.01215295290792761 > 1.0e-15 )
    {
        isLibrationPointComputationErroneous = true;

        cerr << "The computation of the mass parameter does not function well." << endl;
    }

    // Declare and initialize Earth-Moon mass parameter from (James, 2006).
    double earthMoonMassParameter = 0.012277471;

    // Test 2: Compute location of L1 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006). Set boolean to true if test fails.
    isLibrationPointComputationErroneous = testL1LibrationPointLocation(
                isLibrationPointComputationErroneous, earthMoonMassParameter );

    // Test 3: Compute location of L2 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    isLibrationPointComputationErroneous = testL2LibrationPointLocation(
                isLibrationPointComputationErroneous, earthMoonMassParameter );


    // Test 4: Compute location of L3 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    isLibrationPointComputationErroneous = testL3LibrationPointLocation(
                isLibrationPointComputationErroneous, earthMoonMassParameter );

    // Test 5: Compute location of L4 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    isLibrationPointComputationErroneous = testL4LibrationPointLocation(
                isLibrationPointComputationErroneous, earthMoonMassParameter );

    // Test 6: Compute location of L5 Lagrange libration point. Benchmark data is obtained from
    //         (James, 2006).
    isLibrationPointComputationErroneous = testL5LibrationPointLocation(
                isLibrationPointComputationErroneous, earthMoonMassParameter );

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isLibrationPointComputationErroneous )
    {
        cerr << "testLibrationPointLocations failed!" << endl;
    }

    return isLibrationPointComputationErroneous;
}
