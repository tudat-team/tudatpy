/*! \file unitTestBasicMathematicsFunctions.cpp
 *    Source file that defines the unitTestBasicMathematicsFunctions unit test,
 *    containing all basic mathematics functions contained in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Gondelach
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.gondelach@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : T. Secretin
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : T.A.LeitePintoSecretin@student.tudelft.nl
 *
 *    Date created      : 7 February, 2011
 *    Last modified     : 18 January, 2012
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical() function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean()
 *                                  and computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120118    D. Gondelach      Added unit tests for convertCylindricalToCartesianCoordinates,
 *                                  convertCylindricalToCartesianState,
 *                                  convertCartesianToCylindricalCoordinates and
 *                                  convertCartesianToCylindricalState functions.
 *                                  Removed unit test for convertCylindricalToCartesian function.
 */

// Include statements.
#include <cmath>
#include <Eigen/Core>
#include <iostream> 
#include <limits>
#include "Astrodynamics/States/cartesianPositionElements.h"
#include "Mathematics/basicMathematicsFunctions.h"

//! Test implementation of basic mathematics functions.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::fabs;
    using std::acos;
    using std::asin;
    using std::atan2;
    using std::cos;
    using std::sin;
    using tudat::mathematics::computeLinearInterpolation;
    using tudat::mathematics::computeModulo;
    using tudat::mathematics::convertSphericalToCartesian;
    using tudat::mathematics::convertCartesianToSpherical;
    using tudat::mathematics::convertCartesianToCylindricalCoordinates;
    using tudat::mathematics::convertCartesianToCylindricalState;
    using tudat::mathematics::convertCylindricalToCartesianCoordinates;
    using tudat::mathematics::convertCylindricalToCartesianState;
    using namespace tudat;

    // Declare and initialize test result to false.
    bool isBasicMathematicsFunctionsErroneous = false;

    // Test modulo function.
    // Test 10: Test 0.0 mod 0.0.
    // Test 11: Test 2.0 mod 0.0.
    // Test 12: Test 2.0 mod 2.0.
    // Test 13: Test 3.0 mod 2.5.
    // Test 14: Test 3.0 mod -2.5.
    double resultUsingModuloFunction = computeModulo( 0.0, 0.0 );

    if ( fabs( resultUsingModuloFunction - 0.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeModulo function does not "
             << "function correctly, as the computed value: "
             << resultUsingModuloFunction
             << " does not match the expected value: " << 0.0 << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    resultUsingModuloFunction = computeModulo( 2.0, 0.0 );

    if ( fabs( resultUsingModuloFunction - 2.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeModulo function does not "
             << "function correctly, as the computed value: "
             << resultUsingModuloFunction
             << " does not match the expected value: " << 2.0 << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    resultUsingModuloFunction = computeModulo( 2.0, 2.0 );

    if ( fabs( resultUsingModuloFunction - 0.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeModulo function does not "
             << "function correctly, as the computed value: "
             << resultUsingModuloFunction
             << " does not match the expected value: " << 0.0 << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    resultUsingModuloFunction = computeModulo( 3.0, 2.5 );

    if ( fabs( resultUsingModuloFunction - 0.5 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeModulo function does not "
             << "function correctly, as the computed value: "
             << resultUsingModuloFunction
             << " does not match the expected value: " << -0.5 << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    resultUsingModuloFunction = computeModulo( 3.0, -2.5 );

    if ( fabs( resultUsingModuloFunction + 2.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeModulo function does not "
             << "function correctly, as the computed value: "
             << resultUsingModuloFunction
             << " does not match the expected value: " << -2.0 << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from spherical to Cartesian coordinates.
    // Test 20: Test conversion of: ( 0.0, 0.0, 0.0 ).
    // Test 21: Test conversion of: ( 2.0, 225, 225 ).
    // Test 22: Test conversion of: ( 2.0, -225, -225 ).
    // Test 23: Test conversion of: ( 2.0, 180, 180 ).
    Eigen::VectorXd cartesianCoordinates3( 3 );

    convertSphericalToCartesian( 0.0, 0.0, 0.0, cartesianCoordinates3 );

    if ( fabs( cartesianCoordinates3( 0 ) + 0.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 1 ) - 0.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 2 ) - 0.0 ) >
         std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesian, function does not "
             << "function correctly. (test1)" << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    convertSphericalToCartesian( 2.0, 225.0 / 180.0 * M_PI,
                                 225.0 / 180.0 * M_PI, cartesianCoordinates3 );

    if ( fabs( cartesianCoordinates3( 0 ) - 1.0 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 1 ) - 1.0 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 2 ) + sqrt( 2.0 ) ) >
         std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesian function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinates3( 0 ) << ", " << cartesianCoordinates3( 1 )
             << " , " << cartesianCoordinates3( 2 ) << " ) do not match the "
             << "expected coordinates: ( " << 1.0 << ", " << 1.0 << ", "
             << -sqrt( 2.0 ) << " )" << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    convertSphericalToCartesian( 2.0, -225.0 / 180.0 * M_PI, -225.0 / 180.0 * M_PI,
                                 cartesianCoordinates3 );

    if ( fabs( cartesianCoordinates3( 0 ) + 1.0 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 1 ) - 1.0 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinates3( 2 ) + sqrt( 2.0 ) ) >
         std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesian function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinates3( 0 ) << ", " << cartesianCoordinates3( 1 )
             << " , " << cartesianCoordinates3( 2 ) << " ) do not match the "
             << "expected coordinates: ( " << -1.0 << ", " << 1.0 << ", "
             << -sqrt( 2.0 ) << " )" << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    convertSphericalToCartesian( 2.0, 180.0 / 180.0 * M_PI,
                                 180.0 / 180.0 * M_PI, cartesianCoordinates3 );

    if ( fabs( cartesianCoordinates3( 0 ) - 0.0 ) > 1.0e-15 ||
         fabs( cartesianCoordinates3( 1 ) - 0.0 ) > 1.0e-15 ||
         fabs( cartesianCoordinates3( 2 ) + 2.0 ) / 2.0 > 1.0e-15 )
    {
        cerr << "The convertCylindricalToCartesian function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinates3( 0 ) << ", " << cartesianCoordinates3( 1 )
             << " , " << cartesianCoordinates3( 2 ) << " ) do not match the "
             << "expected coordinates: ( " << 0.0 << ", " << 0.0 << ", "
             << -2.0 << " )" << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from Cartesian to spherical coordinates.
    // Test 24: Test conversion of: ( 0.0, 0.0, 0.0 ).
    // Test 25: Test conversion of: ( 2.0, 3.5, -4.1 ).
    // Test 26: Test conversion of: ( 5.2, -6.3, 0.0 ).
    // Test 27: Test conversion of: ( 0.0, 12.2, -0.9 ).

    // Cartesian coordinates.
    Eigen::VectorXd cartesianCoordinatesTest24_( 3 );
    cartesianCoordinatesTest24_.setZero( 3 );

    Eigen::VectorXd cartesianCoordinatesTest25_( 3 );
    cartesianCoordinatesTest25_( 0 ) = 2.0;
    cartesianCoordinatesTest25_( 1 ) = 3.5;
    cartesianCoordinatesTest25_( 2 ) = -4.1;

    Eigen::VectorXd cartesianCoordinatesTest26_( 3 );
    cartesianCoordinatesTest26_( 0 ) = 5.2;
    cartesianCoordinatesTest26_( 1 ) = -6.3;
    cartesianCoordinatesTest26_( 2 ) = 0.0;

    Eigen::VectorXd cartesianCoordinatesTest27_( 3 );
    cartesianCoordinatesTest27_( 0 ) = 0.0;
    cartesianCoordinatesTest27_( 1 ) = 12.2;
    cartesianCoordinatesTest27_( 2 ) = -0.9;

    // Expected vectors in spherical coordinates.
    Eigen::VectorXd expectedSphericalCoordinatesTest24_( 3 );
    expectedSphericalCoordinatesTest24_.setZero( 3 );

    Eigen::VectorXd expectedSphericalCoordinatesTest25_( 3 );
    expectedSphericalCoordinatesTest25_( 0 )
            = sqrt( pow( cartesianCoordinatesTest25_( 0 ), 2.0 )
                    + pow( cartesianCoordinatesTest25_( 1 ), 2.0 )
                    + pow( cartesianCoordinatesTest25_( 2 ), 2.0 ) );
    expectedSphericalCoordinatesTest25_( 1 )
            = atan2( cartesianCoordinatesTest25_( 1 ),
                     cartesianCoordinatesTest25_( 0 ) );
    expectedSphericalCoordinatesTest25_( 2 )
            = acos( cartesianCoordinatesTest25_( 2 ) /
                     expectedSphericalCoordinatesTest25_( 0 ) );

    Eigen::VectorXd expectedSphericalCoordinatesTest26_( 3 );
    expectedSphericalCoordinatesTest26_( 0 )
            = sqrt( pow( cartesianCoordinatesTest26_( 0 ), 2.0 )
                    + pow( cartesianCoordinatesTest26_( 1 ), 2.0 )
                    + pow( cartesianCoordinatesTest26_( 2 ), 2.0 ) );
    expectedSphericalCoordinatesTest26_( 1 )
            = atan2( cartesianCoordinatesTest26_( 1 ),
                     cartesianCoordinatesTest26_( 0 ) );
    expectedSphericalCoordinatesTest26_( 2 )
            = acos( cartesianCoordinatesTest26_( 2 ) /
                     expectedSphericalCoordinatesTest26_( 0 ) );

    Eigen::VectorXd expectedSphericalCoordinatesTest27_( 3 );
    expectedSphericalCoordinatesTest27_( 0 )
            = sqrt( pow( cartesianCoordinatesTest27_( 0 ), 2.0 )
                    + pow( cartesianCoordinatesTest27_( 1 ), 2.0 )
                    + pow( cartesianCoordinatesTest27_( 2 ), 2.0 ) );
    expectedSphericalCoordinatesTest27_( 1 )
            = atan2( cartesianCoordinatesTest27_( 1 ),
                     cartesianCoordinatesTest27_( 0 ) );
    expectedSphericalCoordinatesTest27_( 2 )
            = acos( cartesianCoordinatesTest27_( 2 ) /
                     expectedSphericalCoordinatesTest27_( 0 ) );

    // Result vector in spherical coordinates.
    Eigen::VectorXd sphericalCoordinates_( 3 );

    // Test 24: Test conversion of: ( 0.0, 0.0, 0.0 ).

    // Declare absolute and relative differences.
    double absoluteDifference_;
    double relativeDifference_;

    // Compute conversions.
    convertCartesianToSpherical( cartesianCoordinatesTest24_,
                                 sphericalCoordinates_ );

    // Compute absolute and relative differences.
    absoluteDifference_ = fabs( sphericalCoordinates_.norm( )
            - expectedSphericalCoordinatesTest24_.norm( ) );

    relativeDifference_ = absoluteDifference_
            / expectedSphericalCoordinatesTest24_.norm( );

    // Check if relative error is too large.
    if ( relativeDifference_  > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToSpherical, function does not "
             << "function correctly. ( Test 24 )." << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test 25: Test conversion of: ( 2.0, 3.5, -4.1 ).

    // Compute conversions.
    convertCartesianToSpherical( cartesianCoordinatesTest25_,
                                 sphericalCoordinates_ );

    // Compute absolute and relative differences.
    absoluteDifference_ = fabs( sphericalCoordinates_.norm( )
                               - expectedSphericalCoordinatesTest25_.norm( ) );

    relativeDifference_ = absoluteDifference_
            / expectedSphericalCoordinatesTest25_.norm( );

    if ( relativeDifference_ > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToSpherical, function does not "
             << "function correctly. ( Test 25 )." << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test 26: Test conversion of: ( 5.2, -6.3, 0.0 ).

    // Compute conversion.
    convertCartesianToSpherical( cartesianCoordinatesTest26_, sphericalCoordinates_ );

    // Compute absolute and relative differences.
    absoluteDifference_ = fabs( sphericalCoordinates_.norm( )
                               - expectedSphericalCoordinatesTest26_.norm( ) );

    relativeDifference_ = absoluteDifference_ / expectedSphericalCoordinatesTest26_.norm( );

    if ( relativeDifference_ > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToSpherical, function does not "
             << "function correctly. ( Test 26 )." << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test 27: Test conversion of: ( 0.0, 12.2, -0.9 ).

    // Compute conversion.
    convertCartesianToSpherical( cartesianCoordinatesTest27_, sphericalCoordinates_ );

    // Compute absolute and relative differences.
    absoluteDifference_ = fabs( sphericalCoordinates_.norm( )
                               - expectedSphericalCoordinatesTest27_.norm( ) );

    relativeDifference_ = absoluteDifference_
            / expectedSphericalCoordinatesTest27_.norm( );

    if ( relativeDifference_ > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToSpherical, function does not "
             << "function correctly. ( Test 27 )." << endl;
        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test linear interpolation from independent and dependent vectors.
    // Test 28: Test vector data with expected result of 0.0.
    // Test 29: Test vector data with expected result of -20.5

    // Declare variables.
    // Vectors of data.
    Eigen::VectorXd sortedIndependentVariables( 3 );
    Eigen::VectorXd associatedDependentVariables( 3 );

    // Target independent value in vector data.
    double targetIndependentVariableValue = 0.0;

    // Interpolated value.
    double interpolatedValue;

    // Initialize variables.
    sortedIndependentVariables( 0 ) = 0.0;
    sortedIndependentVariables( 1 ) = 1.0;
    sortedIndependentVariables( 2 ) = 3.0;
    associatedDependentVariables( 0 ) = -20.0;
    associatedDependentVariables( 1 ) = 20.0;
    associatedDependentVariables( 2 ) = 21.0;
    targetIndependentVariableValue = 0.5;

    // Compute interpolation.
    interpolatedValue = computeLinearInterpolation(
            sortedIndependentVariables, associatedDependentVariables,
            targetIndependentVariableValue );

    if ( fabs( interpolatedValue - 0.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
       cerr << "The computeLinearInterpolation function for vector data does "
            << "not function correctly, as the computed value: "
            << interpolatedValue << " does not match the expected value: "
            << 0.0 << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }

    // Reinitialize vectors and target with new data.
    sortedIndependentVariables( 0 ) = 0.0;
    sortedIndependentVariables( 1 ) = 1.0;
    sortedIndependentVariables( 2 ) = 3.0;
    associatedDependentVariables( 0 ) = -20.0;
    associatedDependentVariables( 1 ) = 20.0;
    associatedDependentVariables( 2 ) = 21.0;
    targetIndependentVariableValue = 2.0;

    // Compute interpolation.
    interpolatedValue = computeLinearInterpolation(
            sortedIndependentVariables, associatedDependentVariables,
            targetIndependentVariableValue );

    if ( fabs( interpolatedValue - 20.5 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeLinearInterpolation function for vector data does "
             << "not function correctly, as the computed value: "
             << interpolatedValue << " does not match the expected value: "
             << 20.5 << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }

    // Test linear interpolation with map of vectors with keys as
    // independent variable.
    // Test 30: Test map of data.

    // Declare map of data and vectors for map value.
    std::map < double, Eigen::VectorXd > sortedIndepedentAndDependentVariables;
    Eigen::VectorXd vectorOne( 3 );
    Eigen::VectorXd vectorTwo( 3 );
    Eigen::VectorXd vectorThree( 3 );
    Eigen::VectorXd interpolatedVector( 3 );

    // Initialize vectors for map value.
    // Initialize first vector.
    vectorOne( 0 ) = 10.0;
    vectorOne( 1 ) = -10.0;
    vectorOne( 2 ) = 70.0;

    // Initialize second vector.
    vectorTwo( 0 ) = 20.0;
    vectorTwo( 1 ) = -5.0;
    vectorTwo( 2 ) = 80.0;

    // Initialize third vector.
    vectorThree( 0 ) = 30.0;
    vectorThree( 1 ) = 60.0;
    vectorThree( 2 ) = 90.0;

    // Set map values in map using vector data.
    sortedIndepedentAndDependentVariables[ 0.0 ] = vectorOne;
    sortedIndepedentAndDependentVariables[ 1.0 ] = vectorTwo;
    sortedIndepedentAndDependentVariables[ 2.0 ] = vectorThree;

    // Set target independent variable value for interpolation.
    targetIndependentVariableValue = 1.5;

    // Compute interpolation.
    interpolatedVector = computeLinearInterpolation(
            sortedIndepedentAndDependentVariables,
            targetIndependentVariableValue );

    if ( fabs( interpolatedVector( 0 ) - 25.0 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedVector( 1 ) - 27.5 ) > std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedVector( 2 ) - 85.0 ) > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeLinearInterpolation function for a map of "
             << "vectors, does not function correctly, as the compute vector "
             << "data: " <<  interpolatedVector << " does not match the "
             << "expected vector data: ( 25.0, 27.5, 85.0 ). " << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }


    // Test linear interpolation with map of States.
    // Test 31: Test run using three states. Interpolation between the first
    // and second, and second and third tested. An extrapolation is also
    // tested.

    // Declare variables.
    std::map < double, State* > stateMap;
    State* interpolatedState;
    CartesianPositionElements* cartesianPostitionStateOne =
            new CartesianPositionElements();
    CartesianPositionElements* cartesianPostitionStateTwo =
            new CartesianPositionElements();
    CartesianPositionElements* cartesianPostitionStateThree =
            new CartesianPositionElements();

    // Set position elements for three states.
    cartesianPostitionStateOne->setCartesianElementX( -7.0 );
    cartesianPostitionStateOne->setCartesianElementY( -3.0 );
    cartesianPostitionStateOne->setCartesianElementZ( -1.0 );

    cartesianPostitionStateTwo->setCartesianElementX( 7.0 );
    cartesianPostitionStateTwo->setCartesianElementY( 3.0 );
    cartesianPostitionStateTwo->setCartesianElementZ( 1.0 );

    cartesianPostitionStateThree->setCartesianElementX( 0.0 );
    cartesianPostitionStateThree->setCartesianElementY( 0.0 );
    cartesianPostitionStateThree->setCartesianElementZ( 0.0 );

    // Add states to the state map with key the independent variable.
    stateMap[ -2.0 ] = cartesianPostitionStateOne;
    stateMap[ 2.0 ] = cartesianPostitionStateTwo;
    stateMap[ 4.0 ] = cartesianPostitionStateThree;

    // Set target indepedent variable value.
    targetIndependentVariableValue = 0.0;

    // Compute interpolation.
    interpolatedState = computeLinearInterpolation(
            stateMap, targetIndependentVariableValue );

    if ( fabs( interpolatedState->state( 0 ) - 0.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 1 ) - 0.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 2)  - 0.0 ) >
         std::numeric_limits< double >::epsilon( ) )
    {
       cerr << "The computeLinearInterpolation function for a map of State "
            << "objects does not function correctly, as the computed state: "
            << interpolatedState->state << " does not match the expected "
            << "state: ( 0.0, 0.0, 0.0 )" << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }

    // Set new target and compute again.
    targetIndependentVariableValue = 3.0;

    // Compute interpolation.
    interpolatedState = computeLinearInterpolation(
            stateMap, targetIndependentVariableValue );

    if ( fabs( interpolatedState->state( 0 ) - 3.5 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 1 ) - 1.5 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 2 ) - 0.5 ) >
         std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeLinearInterpolation function for a map of State "
             << "objects does not function correctly, as the computed state: "
             << interpolatedState->state << " does not match the expected "
             << "state: ( 3.5, 1.5, 0.5 )" << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }

    // Set new target and compute again.
    targetIndependentVariableValue = 6.0;

    // Compute interpolation.
    interpolatedState = computeLinearInterpolation(
            stateMap, targetIndependentVariableValue );

    if ( fabs( interpolatedState->state( 0 ) + 7.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 1 ) + 3.0 ) >
         std::numeric_limits< double >::epsilon( ) ||
         fabs( interpolatedState->state( 2 ) + 1.0 ) >
         std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The computeLinearInterpolation function for a map of State "
             << "objects does not function correctly, as the computed state: "
             << interpolatedState->state << " does not match the expected "
             << "state: ( -7.0, -3.0, -1.0 )" << endl;
       isBasicMathematicsFunctionsErroneous = true;
    }

    // Test computation of sample mean and sample variance.
    // Test 32: Test computation of sample mean and sample variance on finite
    // population using unbiased estimators. The expected values are computed
    // using the Microsoft Excel AVERAGE() and VAR() functions.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Declare expected sample mean.
    double expectedSampleMean = 9.1;

    // Declare expected sample variance.
    double expectedSampleVariance = 24.665;

    // Declare and compute sample mean.
    double computedSampleMean = mathematics::computeSampleMean( sampleData );

    // Declare and compute sample variance.
    double computedSampleVariance
            = mathematics::computeSampleVariance( sampleData );

    // Check if differences between computed and expected sample means and
    // variances are too large; if so print cerr statements.
    if ( fabs( computedSampleMean - expectedSampleMean) / expectedSampleMean
         > std::numeric_limits< double >::epsilon( )
         || fabs( computedSampleVariance - expectedSampleVariance )
         / expectedSampleVariance > std::numeric_limits< double >::epsilon( ) )
    {
        isBasicMathematicsFunctionsErroneous = true;

        cerr << "The computeSampleMean() and/or computeSampleVariance "
             << "functions are erroneous, as the computed sample mean "
             << "( " << computedSampleMean << " ) and/or computed sample "
             << "variance ( " << computedSampleVariance << " ) are/is not "
             << "equal to the expected sample mean ( " << expectedSampleMean
             << " ) and/or expected sample variance ( "
             << expectedSampleVariance << " )." << endl;
    }

    // Test conversion from cylindrical (r, theta, z) to Cartesian (x, y, z) coordinates.
    // Test 33: Test conversion of ( 0.0, pi, 1.2 ).
    // Test 34: Test conversion of ( 2.3, -pi/2, -3.5 ).

    // Cylindrical coordinates (r,theta,z).
    Eigen::Vector3d cylindricalCoordinatesTest33, cylindricalCoordinatesTest34;
    cylindricalCoordinatesTest33 << 0.0, M_PI, 1.2;
    cylindricalCoordinatesTest34 << 2.3, -M_PI / 2.0, -3.5;

    // Expected Cartesian coordinates (x, y, z).
    Eigen::Vector3d expectedCartesianCoordinatesTest33, expectedCartesianCoordinatesTest34;
    expectedCartesianCoordinatesTest33 << 0.0, 0.0, 1.2;
    expectedCartesianCoordinatesTest34 << 2.3 * cos( -M_PI / 2.0 ), 2.3*sin( -M_PI / 2.0 ), -3.5;

    // Test Cartesian coordinates (x, y, z).
    Eigen::Vector3d cartesianCoordinatesTest33 = convertCylindricalToCartesianCoordinates(
                                                        cylindricalCoordinatesTest33 );
    Eigen::Vector3d cartesianCoordinatesTest34 = convertCylindricalToCartesianCoordinates(
                                                        cylindricalCoordinatesTest34 );

    if ( fabs( cartesianCoordinatesTest33(0) - expectedCartesianCoordinatesTest33(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinatesTest33(1) - expectedCartesianCoordinatesTest33(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinatesTest33(2) - expectedCartesianCoordinatesTest33(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesianCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinatesTest33(0) << ", " << cartesianCoordinatesTest33(1)
             << ", " << cartesianCoordinatesTest33(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianCoordinatesTest33(0) << ", "
             << expectedCartesianCoordinatesTest33(1) << ", "
             << expectedCartesianCoordinatesTest33(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cartesianCoordinatesTest34(0) - expectedCartesianCoordinatesTest34(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinatesTest34(1) - expectedCartesianCoordinatesTest34(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianCoordinatesTest34(2) - expectedCartesianCoordinatesTest34(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesianCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinatesTest34(0) << ", " << cartesianCoordinatesTest34(1)
             << ", " << cartesianCoordinatesTest34(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianCoordinatesTest34(0) << ", "
             << expectedCartesianCoordinatesTest34(1) << ", "
             << expectedCartesianCoordinatesTest34(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from cylindrical (r, theta, z, Vr, Vtheta, Vz) to
    // Cartesian (x, y, z, xdot, ydot, zdot) state.
    // Test 35: Test conversion of (2.1, pi/2.0, 1.2, 5.4, 4.5, -3.9)
    // Test 36: Test conversion of (0.0, 8.2*pi/3.0, -2.5, -5.8, 0.0, 1.7)

    // Cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd cylindricalStateTest35( 6 ), cylindricalStateTest36( 6 );
    cylindricalStateTest35 << 2.1, M_PI / 2.0, 1.2, 5.4, 4.5, -3.9;
    cylindricalStateTest36 << 0.0, 8.2 * M_PI / 3.0, -2.5, -5.8, 0.0, 1.7;

    // Expected Cartesian states (x, y, z, xdot, ydot, zdot).
    Eigen::VectorXd expectedCartesianStateTest35( 6 ), expectedCartesianStateTest36( 6 );
    expectedCartesianStateTest35 << 2.1 * cos( M_PI / 2.0 ),
                                    2.1 * sin( M_PI / 2.0 ),
                                    1.2,
                                    5.4 * cos( M_PI / 2.0 ) - 4.5 * sin( M_PI / 2.0 ),
                                    5.4 * sin( M_PI / 2.0 ) + 4.5 * cos( M_PI / 2.0 ),
                                    -3.9;
    expectedCartesianStateTest36 << 0.0,
                                    0.0,
                                    -2.5,
                                    -5.8 * cos( 8.2 * M_PI / 3.0 ),
                                    -5.8 * sin( 8.2 * M_PI / 3.0 ),
                                    1.7;

    // Test Cartesian states (x, y, z, xdot, ydot, zdot).
    Eigen::VectorXd cartesianStateTest35 = convertCylindricalToCartesianState(
                                                        cylindricalStateTest35 );
    Eigen::VectorXd cartesianStateTest36 = convertCylindricalToCartesianState(
                                                        cylindricalStateTest36 );

    if ( fabs( cartesianStateTest35(0) - expectedCartesianStateTest35(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest35(1) - expectedCartesianStateTest35(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest35(2) - expectedCartesianStateTest35(2) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest35(3) - expectedCartesianStateTest35(3) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest35(4) - expectedCartesianStateTest35(4) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest35(5) - expectedCartesianStateTest35(5) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesianState function does not "
             << "function correctly, as the computed sate: ( "
             << cartesianStateTest35(0) << ", " << cartesianStateTest35(1) << ", "
             << cartesianStateTest35(2) << ", " << cartesianStateTest35(3) << ", "
             << cartesianStateTest35(4) << ", " << cartesianStateTest35(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianStateTest35(0) << ", " << expectedCartesianStateTest35(1) << ", "
             << expectedCartesianStateTest35(2) << ", " << expectedCartesianStateTest35(3) << ", "
             << expectedCartesianStateTest35(4) << ", " << expectedCartesianStateTest35(5)
             << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cartesianStateTest36(0) - expectedCartesianStateTest36(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest36(1) - expectedCartesianStateTest36(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest36(2) - expectedCartesianStateTest36(2) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest36(3) - expectedCartesianStateTest36(3) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest36(4) - expectedCartesianStateTest36(4) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cartesianStateTest36(5) - expectedCartesianStateTest36(5) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCylindricalToCartesianState function does not "
             << "function correctly, as the computed sate: ( "
             << cartesianStateTest36(0) << ", " << cartesianStateTest36(1) << ", "
             << cartesianStateTest36(2) << ", " << cartesianStateTest36(3) << ", "
             << cartesianStateTest36(4) << ", " << cartesianStateTest36(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianStateTest36(0) << ", " << expectedCartesianStateTest36(1) << ", "
             << expectedCartesianStateTest36(2) << ", " << expectedCartesianStateTest36(3) << ", "
             << expectedCartesianStateTest36(4) << ", " << expectedCartesianStateTest36(5)
             << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from Cartesian (x, y, z) to cylindrical (r, theta, z) coordinates.
    // Test 37: Test conversion of ( 0.0, 0.0, 1.0 ).
    // Test 38: Test conversion of ( 0.0, 2.0, 1.0 ).
    // Test 39: Test conversion of ( 0.0, -2.0, -1.0 ).
    // Test 40: Test conversion of ( -5.0, -8.0, 5.0 ).

    // Cartesian coordinates (x, y, z).
    Eigen::Vector3d cartesianCoordinatesTest37, cartesianCoordinatesTest38,
                    cartesianCoordinatesTest39, cartesianCoordinatesTest40;
    cartesianCoordinatesTest37 << 0.0, 0.0, 0.0;
    cartesianCoordinatesTest38 << 0.0, 2.0, 1.0;
    cartesianCoordinatesTest39 << 0.0, -2.0, -1.0;
    cartesianCoordinatesTest40 << -5.0, -8.0, 5.0;

    // Expected cylindrical coordinates (r, theta, z).
    Eigen::Vector3d expectedCylindricalCoordinatesTest37, expectedCylindricalCoordinatesTest38,
                    expectedCylindricalCoordinatesTest39, expectedCylindricalCoordinatesTest40;
    expectedCylindricalCoordinatesTest37 << 0.0, 0.0, 0.0;
    expectedCylindricalCoordinatesTest38 << 2.0, M_PI / 2.0, 1.0;
    expectedCylindricalCoordinatesTest39 << 2.0, 3.0 * M_PI / 2.0, -1.0;
    expectedCylindricalCoordinatesTest40 << sqrt( 25.0 + 64.0 ), atan2(-8.0,-5.0) + 2.0 * M_PI, 5.0;

    // Test cylindrical coordinates (r, theta, z).
    Eigen::Vector3d cylindricalCoordinatesTest37 = convertCartesianToCylindricalCoordinates(
                                                                    cartesianCoordinatesTest37 );
    Eigen::Vector3d cylindricalCoordinatesTest38 = convertCartesianToCylindricalCoordinates(
                                                                    cartesianCoordinatesTest38 );
    Eigen::Vector3d cylindricalCoordinatesTest39 = convertCartesianToCylindricalCoordinates(
                                                                    cartesianCoordinatesTest39 );
    Eigen::Vector3d cylindricalCoordinatesTest40 = convertCartesianToCylindricalCoordinates(
                                                                    cartesianCoordinatesTest40 );

    if ( fabs( cylindricalCoordinatesTest37(0) - expectedCylindricalCoordinatesTest37(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest37(1) - expectedCylindricalCoordinatesTest37(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest37(2) - expectedCylindricalCoordinatesTest37(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest37(0) << ", " << cylindricalCoordinatesTest37(1)
             << ", " << cylindricalCoordinatesTest37(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest37(0) << ", "
             << expectedCylindricalCoordinatesTest37(1) << ", "
             << expectedCylindricalCoordinatesTest37(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest38(0) - expectedCylindricalCoordinatesTest38(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest38(1) - expectedCylindricalCoordinatesTest38(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest38(2) - expectedCylindricalCoordinatesTest38(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest38(0) << ", " << cylindricalCoordinatesTest38(1)
             << ", " << cylindricalCoordinatesTest38(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest38(0) << ", "
             << expectedCylindricalCoordinatesTest38(1) << ", "
             << expectedCylindricalCoordinatesTest38(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest39(0) - expectedCylindricalCoordinatesTest39(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest39(1) - expectedCylindricalCoordinatesTest39(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest39(2) - expectedCylindricalCoordinatesTest39(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest39(0) << ", " << cylindricalCoordinatesTest39(1)
             << ", " << cylindricalCoordinatesTest39(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest39(0) << ", "
             << expectedCylindricalCoordinatesTest39(1) << ", "
             << expectedCylindricalCoordinatesTest39(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest40(0) - expectedCylindricalCoordinatesTest40(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest40(1) - expectedCylindricalCoordinatesTest40(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalCoordinatesTest40(2) - expectedCylindricalCoordinatesTest40(2) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalCoordinates function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest40(0) << ", " << cylindricalCoordinatesTest40(1)
             << ", " << cylindricalCoordinatesTest40(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest40(0) << ", "
             << expectedCylindricalCoordinatesTest40(1) << ", "
             << expectedCylindricalCoordinatesTest40(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from Cartesian (x, y, z, xdot, ydot, zdot) to
    // cylindrical (r, theta, z, Vr, Vtheta, Vz) state.
    // Test 41: Test conversion of ( 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).
    // Test 42: Test conversion of ( 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).
    // Test 43: Test conversion of ( -7.0, -4.0, 3.0, 5.0, -3.0, 7.0 ).

    // Cartesian states (x,y,z,xdot,ydot,zdot).
    Eigen::VectorXd cartesianStateTest41( 6 ), cartesianStateTest42( 6 ),
                    cartesianStateTest43( 6 );
    cartesianStateTest41 << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0;
    cartesianStateTest42 << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0;
    cartesianStateTest43 << -7.0, -4.0, 3.0, 5.0, -3.0, 7.0;

    // Expected cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd expectedCylindricalStateTest41( 6 ), expectedCylindricalStateTest42( 6 ),
                    expectedCylindricalStateTest43( 6 );
    expectedCylindricalStateTest41 << 0.0, 0.0, 1.0, sqrt( 25.0 + 36.0 ), 0.0, -9.0;
    expectedCylindricalStateTest42 << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0;
    expectedCylindricalStateTest43 << sqrt(49.0+16.0), atan2( -4.0, -7.0 ) + 2.0 * M_PI, 3.0,
            ( -7.0 * 5.0 + ( -4.0 ) * -3.0 ) / sqrt( 49.0 + 16.0),
            ( -7.0 * -3.0 - ( -4.0 ) * 5.0 ) / sqrt( 49.0 + 16.0), 7.0;

    // Test cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd cylindricalStateTest41 = convertCartesianToCylindricalState(
                                                            cartesianStateTest41 );
    Eigen::VectorXd cylindricalStateTest42 = convertCartesianToCylindricalState(
                                                            cartesianStateTest42 );
    Eigen::VectorXd cylindricalStateTest43 = convertCartesianToCylindricalState(
                                                            cartesianStateTest43 );

    if ( fabs( cylindricalStateTest41(0) - expectedCylindricalStateTest41(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest41(1) - expectedCylindricalStateTest41(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest41(2) - expectedCylindricalStateTest41(2) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest41(3) - expectedCylindricalStateTest41(3) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest41(4) - expectedCylindricalStateTest41(4) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest41(5) - expectedCylindricalStateTest41(5) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalState function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest41(0) << ", " << cylindricalStateTest41(1) << ", "
             << cylindricalStateTest41(2) << ", " << cylindricalStateTest41(3) << ", "
             << cylindricalStateTest41(4) << ", " << cylindricalStateTest41(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest41(0) << ", "
             << expectedCylindricalStateTest41(1) << ", "
             << expectedCylindricalStateTest41(2) << ", "
             << expectedCylindricalStateTest41(3) << ", "
             << expectedCylindricalStateTest41(4) << ", "
             << expectedCylindricalStateTest41(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalStateTest42(0) - expectedCylindricalStateTest42(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest42(1) - expectedCylindricalStateTest42(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest42(2) - expectedCylindricalStateTest42(2) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest42(3) - expectedCylindricalStateTest42(3) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest42(4) - expectedCylindricalStateTest42(4) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest42(5) - expectedCylindricalStateTest42(5) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalState function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest42(0) << ", " << cylindricalStateTest42(1) << ", "
             << cylindricalStateTest42(2) << ", " << cylindricalStateTest42(3) << ", "
             << cylindricalStateTest42(4) << ", " << cylindricalStateTest42(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest42(0) << ", "
             << expectedCylindricalStateTest42(1) << ", "
             << expectedCylindricalStateTest42(2) << ", "
             << expectedCylindricalStateTest42(3) << ", "
             << expectedCylindricalStateTest42(4) << ", "
             << expectedCylindricalStateTest42(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalStateTest43(0) - expectedCylindricalStateTest43(0) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest43(1) - expectedCylindricalStateTest43(1) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest43(2) - expectedCylindricalStateTest43(2) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest43(3) - expectedCylindricalStateTest43(3) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest43(4) - expectedCylindricalStateTest43(4) )
         > std::numeric_limits< double >::epsilon( ) ||
         fabs( cylindricalStateTest43(5) - expectedCylindricalStateTest43(5) )
         > std::numeric_limits< double >::epsilon( ) )
    {
        cerr << "The convertCartesianToCylindricalState function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest43(0) << ", " << cylindricalStateTest43(1) << ", "
             << cylindricalStateTest43(2) << ", " << cylindricalStateTest43(3) << ", "
             << cylindricalStateTest43(4) << ", " << cylindricalStateTest43(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest43(0) << ", "
             << expectedCylindricalStateTest43(1) << ", "
             << expectedCylindricalStateTest43(2) << ", "
             << expectedCylindricalStateTest43(3) << ", "
             << expectedCylindricalStateTest43(4) << ", "
             << expectedCylindricalStateTest43(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isBasicMathematicsFunctionsErroneous )
    {
        cerr << "testBasicMathematicsFunctions failed!" << endl;
    }

    return isBasicMathematicsFunctionsErroneous;
}

// End of file.
