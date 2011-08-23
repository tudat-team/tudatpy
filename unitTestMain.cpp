/*! \file unitTestMain.cpp
 *    Source file for the main unit test executable that runs all unit
 *    tests and generates a unit_test_report.txt report.
 *    A 0 corresponds to success in the report.
 *
 *    Path              : /
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : B. Römgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 January, 2011
 *    Last modified     : 22 August, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110125    B. Römgens        File created.
 *      110217    B. Römgens        Added new unit tests for v0.2.
 *      110217    K. Kumar          Minor changes.
 *      110329    K. Kumar          Added testTextFileReader unit test.
 *      110701    K. Kumar          Added testApproximatePlanetPositions,
 *                                  testDeepSpaceManeuver, testEulerIntegrator,
 *                                  testRungeKutta4thOrderFixedStepsizeIntegrator,
 *                                  unit tests; updated layout.
 *      110822    D.Dirkx           Added testAerodynamicForcesAndMoments
 */

// Include statements.
#include <fstream>
#include <ctime>

// Mathematics unit test includes.
#include "unitTestBasicMathematicsFunctions.h"
#include "unitTestEulerIntegrator.h"
#include "unitTestLawgsSurfaceGeometry.h"
#include "unitTestCubicSplineInterpolation.h"
#include "unitTestNewtonRaphson.h"
#include "unitTestRandomNumberGenerator.h"
#include "unitTestRungeKutta4thOrderFixedStepsizeIntegrator.h"
#include "unitTestSimpleLinearRegression.h"
#include "unitTestUnitConversions.h"

// Astrodynamics unit test includes.
#include "unitTestAerodynamicMomentAndAerodynamicForce.h"
#include "unitTestAerodynamicsNamespace.h"
#include "unitTestApproximatePlanetPositions.h"
#include "unitTestCartesianElements.h"
#include "unitTestCoefficientGenerator.h"
#include "unitTestDeepSpaceManeuver.h"
#include "unitTestEscapeAndCapture.h"
#include "unitTestExponentialAtmosphere.h"
#include "unitTestGravityAssist.h"
#include "unitTestKeplerianElements.h"
#include "unitTestKeplerPropagator.h"
#include "unitTestLambertTargeter.h"
#include "unitTestLibrationPoints.h"
#include "unitTestNumericalPropagator.h"
#include "unitTestOrbitalElementConversions.h"
#include "unitTestPhysicalConstants.h"
#include "unitTestReferenceFrameTransformations.h"
#include "unitTestSphericalHarmonicsGravityField.h"
#include "unitTestTabulatedAtmosphere.h"

// Input unit test includes.
#include "unitTestTextFileReader.h"
#include "unitTestTwoLineElementsTextFileReader.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Execute all unit tests.
/*!
 * Executes all unit tests.
 */
int main( )
{
    // Return value, 0 on success, 1 on failure of one or more unit tests.
    int isErroneous = 0;

    // Run all unit tests.
    // Unit test functions return false ( 0 ) for success.

    // Run Mathematics unit tests.

    // testBasicMathematicsFunctions: Tests the following functions,
    // linear interpolation, absolute value, raise to integer power
    // modulo, spherical to Cartesian, cylindrical to Cartesian.
    bool testBasicMathematicsFunctions =
            unit_tests::testBasicMathematicsFunctions( );

    // testCubicSpline: Test the cubic spline interpolation.
    bool testCubicSplineInterpolation = unit_tests::testCubicSplineInterpolation( );

    // testEulerIntegrator:: Tests the Euler integrator against benchmark
    // data from (Burden and Faires, 2001).
    bool testEulerIntegrator = unit_tests::testEulerIntegrator( );

    // testExponentialRandomNumberGenerator: Tests the exponential random
    // number generator.
    bool testExponentialRandomNumberGenerator
            = unit_tests::testExponentialRandomNumberGenerator( );

    // testLawgsSurfaceGeometry: Tests a Lawgs mesh of a sphere.
    bool testLawgsSurfaceGeometry =
            unit_tests::testLawgsSurfaceGeometry( );

    // testNewtonRaphson: Tests the Newton-Raphson root-finder.
    bool testNewtonRaphson = unit_tests::testNewtonRaphsonMethod( );

    // testNormalRandomNumberGenerator: Tests the normal random number
    // generator.
    bool testNormalRandomNumberGenerator
            = unit_tests::testNormalRandomNumberGenerator( );

    // testRungeKutta4thOrderFixedStepsizeIntegrator: Tests the 4th-order,
    // fixed stepsize, Runge-Kutta integrator against benchmark data from
    // (Burden and Faires, 2001).
    bool testRungeKutta4thOrderFixedStepsizeIntegrator
            = unit_tests::testRungeKutta4thOrderFixedStepsizeIntegrator( );

    // testSimpleLinearRegression: Tests simple linear regression method
    // against benchmark data from (Burden and Faires, 2001).
    bool testSimpleLinearRegression
            = unit_tests::testSimpleLinearRegression( );

    // testUniformRandomNumberGenerator: Tests the uniform random number
    // generator.
    bool testUniformRandomNumberGenerator
            = unit_tests::testUniformRandomNumberGenerator( );

    // testUnitConversions: Tests conversions that are defined
    // in unitConversions.h.
    bool testUnitConversions = unit_tests::testUnitConversions( );

    // Run Astrodynamics unit tests.

    // testAerodynamicForcesAndMoments: Tests calculation of aerodynamic forces and moments
    bool testAerodynamicForcesAndMoments =
            unit_tests::testAerodynamicMomentAndAerodynamicForce( );

    // testAerodynamicsNamespace: Tests the following functions:
    // static pressure ratio, stagnation pressure coefficient,
    // modified Newtonian pressure coefficient, Tangent Cone pressure coef.,
    // Mach base pressure coefficient, Tangent Wedge pressure coefficient,
    // freestream Prandtl-Meyer function, vacuum pressure coefficient,
    // shock pressure ratio, shock density ratio, shock temperature ratio
    // shock total pressure ratio, shock wave total pressure ratio.
    bool testAerodynamicsNamespace =
            unit_tests::testAerodynamicsNameSpace( );

    // testApproximatePlanetPositions: Tests the approximate planet positions
    // ephemeris class implemented in the class ApproximatePlanetPositions. The
    // ephemeris data comes from http://ssd.jpl.nasa.gov/horizons.cgi.
    bool testApproximatePlanetPositions
            = unit_tests::testApproximatePlanetPositions( );

    // testCartesianElements: Tests the different set and get functions
    // of CartesianElements.
    bool testCartesianElements =
            unit_tests::testCartesianElements( );

    // testCoefficientGenerator: Tests aerothermodynamic coefficients for
    // a sphere and the apollo capsule.
    bool testCoefficientGenerator =
            unit_tests::testCoefficientGenerator( );

    // testDeepSpaceManeuver:: Tests implementation of a Deep Space Maneuver
    // (DSM).
    bool testDeepSpaceManeuver = unit_tests::testDeepSpaceManeuver( );

    // testEscapeAndCapture: Tests the computed delta-V of escape and capture.
    bool testEscapeAndCapture =
            unit_tests::testEscapeAndCapture( );

    // testExponentialAtmosphere: Tests the exponential atmosphere model.
    bool testExponentialAtmosphere =
            unit_tests::testExponentialAtmosphere( );

    // testGravityAssist: Tests delta-V for a gravity assist for the
    // case of equal hyperbolic excess velocities.
    bool testGravityAssist =
            unit_tests::testGravityAssist( );

    // testKeplerianElements: Tests the different set and get functions
    // of KeplerianElements.
    bool testKeplerianElements =
            unit_tests::testKeplerianElements( );

    // testKeplerPropagator: Tests the Kepler propagator. For this test, the
    // ROOT_PATH variable must be set in basicFunctions.h to the root
    // directory for Tudat.
    bool testKeplerPropagator = unit_tests::testKeplerPropagator( );

    // testLambertTargeter: Tests the values of semi major axis,
    // radial speed at departure, radial speed at arrival,
    // transverse speed at departure, transverse speed at arrival for
    // elliptical and hypterbolic test case.
    bool testLambertTargeter =
            unit_tests::testLambertTargeter( );

    // testLibrationPointComputations: Tests the values of the locations of
    // the five Lagrange libration points in the Earth-Moon system, using
    // benchmark data from (James, 2006). Also tests the computation of the
    // dimensionless mass parameter.
    bool testLibrationPointLocations
            = unit_tests::testLibrationPointLocations( );

    // testNumericalPropagator: Tests the numerical propagator.
    bool testNumericalPropagator = unit_tests::testNumericalPropagator( );

    // testOrbitalElementConversions: Tests the code for elliptical, parabolic,
    // hyperbolic and circular orbits. It also tests the conversion from
    // Cartesian to Keplerian and Keplerian to Cartesian.
    bool testOrbitalElementConversions =
            unit_tests::testOrbitalElementConversions( );

    // testPhysicalConstants: Tests the physical constants that are defined
    // in physicalConstants.h.
    bool testPhysicalConstants = unit_tests::testPhysicalConstants( );

    // TestFrameTransformations: Tests the reference frame transformations:
    // Rotating planetocentric (R) to Inertial (I) frame transformation,
    // Inertial (I) to Rotating planetocentric (R) frame transformation,
    // Aerodynamic (airspeed based) (AA) to Body (B) frame transformation.
    bool testReferenceFrameTransformations
            = unit_tests::testReferenceFrameTransformations( );

    // testSphericalHarmonicsGravityField: Tests the implementation of the
    // spherical harmonics gravity field class.
    bool testSphericalHarmonicsGravityField =
            unit_tests::testSphericalHarmonicsGravityField( );

    // testTabulatedAtmosphere: Test the tabulated atmosphere.
    bool testTabulatedAtmosphere =
            unit_tests::testTabulatedAtmosphere( );

    // Run Input unit tests.

    // testTextFileReader: Tests reading in a text file using the readAndStore,
    // skipLines, and skipLinesStartingWithCharacter functions.
    bool testTextFileReader =
            unit_tests::testTextFileReader( );

    // testTwoLineElementsTextFileReader: Tests reading and storing of TLE
    // catalog data, and a function to check the data integrity of the catalog.
    bool testTwoLineElementsTextFileReader =
            unit_tests::testTwoLineElementsTextFileReader( );

    // Check if all tests were successful and print cerr message stating
    // which test failed, if any.
    if ( testBasicMathematicsFunctions )
    {
        cerr << "testBasicMathematicsFunctions failed!" << endl;
        isErroneous = 1;
    }

    if ( testCubicSplineInterpolation )
    {
        cerr << "cubicSplineInterpolation failed! " << endl;
        isErroneous = 1;
    }

    if ( testEulerIntegrator )
    {
        cerr << "testEulerIntegrator failed!" << endl;
        isErroneous = 1;
    }

    if ( testExponentialRandomNumberGenerator )
    {
        cerr << "testExponentialRandomNumberGenerator failed!" << endl;
        isErroneous = 1;
    }

    if ( testLawgsSurfaceGeometry )
    {
        cerr << "testLawgsSurfaceGeometry failed!" << endl;
        isErroneous = 1;
    }

    if ( testNewtonRaphson )
    {
        cerr << "testNewtonRaphson failed!" << endl;
        isErroneous = 1;
    }

    if ( testNormalRandomNumberGenerator )
    {
        cerr << "testNormalRandomNumberGenerator failed!" << endl;
        isErroneous = 1;
    }

    if ( testRungeKutta4thOrderFixedStepsizeIntegrator )
    {
        cerr << "testRungeKutta4thOrderFixedStepsizeIntegrator failed!"
             << endl;
        isErroneous = 1;
    }

    if ( testSimpleLinearRegression )
    {
        cerr << "testSimpleLinearRegression failed!" << endl;
        isErroneous = 1;
    }

    if ( testUniformRandomNumberGenerator )
    {
        cerr << "testUniformRandomNumberGenerator failed!" << endl;
        isErroneous = 1;
    }

    if ( testUnitConversions )
    {
        cerr << "testUnitConversions failed!" << endl;
        isErroneous = 1;
    }

    if ( testAerodynamicForcesAndMoments )
    {
        cerr << "testAerodynamicForcesAndMoment failed!" << endl;
        isErroneous = 1;
    }

    if ( testAerodynamicsNamespace )
    {
        cerr << "testAerodynamicsNamespace failed!" << endl;
        isErroneous = 1;
    }

    if ( testApproximatePlanetPositions )
    {
        cerr << "testApproximatePlanetPositions failed!" << endl;
        isErroneous = 1;
    }


    if ( testCartesianElements )
    {
        cerr << "testCartesianElements failed!" << endl;
        isErroneous = 1;
    }

    if ( testCoefficientGenerator )
    {
        cerr << "testCoefficientGenerator failed!" << endl;
        isErroneous = 1;
    }

    if ( testDeepSpaceManeuver )
    {
        cerr << "testDeepSpaceManeuver failed!" << endl;
        isErroneous = 1;
    }

    if ( testEscapeAndCapture )
    {
        cerr << "testEscapeAndCapture failed!" << endl;
        isErroneous = 1;
    }

    if ( testExponentialAtmosphere )
    {
        cerr << "testExponentialAtmosphere failed!" << endl;
        isErroneous = 1;
    }

    if ( testGravityAssist )
    {
        cerr << "testGravityAssist failed!" << endl;
        isErroneous = 1;
    }

    if ( testKeplerianElements )
    {
        cerr << "testKeplerianElements failed!" << endl;
        isErroneous = 1;
    }

    if ( testKeplerPropagator )
    {
        cerr << "testKeplerPropagator failed!" << endl;
        isErroneous = 1;
    }

    if ( testLambertTargeter )
    {
        cerr << "testLambertTargeter failed!" << endl;
        isErroneous = 1;
    }

    if ( testLibrationPointLocations )
    {
        cerr << "testLibrationPointLocations failed!" << endl;
        isErroneous = 1;
    }

    if ( testNumericalPropagator )
    {
        cerr << "testNumericalPropagator failed!" << endl;
        isErroneous = 1;
    }

    if ( testOrbitalElementConversions )
    {
        cerr << "testOrbitalElementConversions failed!" << endl;
        isErroneous = 1;
    }

    if ( testPhysicalConstants )
    {
        cerr << "testPhysicalConstants failed!" << endl;
        isErroneous = 1;
    }

    if ( testReferenceFrameTransformations )
    {
        cerr << "referenceFrameTransformations failed" << endl;
        isErroneous = 1;
    }

    if ( testSphericalHarmonicsGravityField )
    {
        cerr << "testSphericalHarmonicsGravityField failed!" << endl;
        isErroneous = 1;
    }

    if ( testTabulatedAtmosphere )
    {
        cerr << "testTabulatedAtmosphere failed!" << endl;
        isErroneous = 1;
    }

    if ( testTextFileReader )
    {
        cerr << "testTextFileReader failed!" << endl;
        isErroneous = 1;
    }

    if ( testTwoLineElementsTextFileReader )
    {
        cerr << "testTwoLineElementsTextFileReader failed!" << endl;
        isErroneous = 1;
    }   

    // Generate unit test report file ( unitTestReport.txt ).

    // Generate time and date stamp for report.
    char cptime[ 50 ];
    time_t now = time( NULL );
    strftime( cptime, 50, "%d %B %Y", localtime( &now ) );
    std::string date = cptime;
    strftime( cptime, 50, "%H:%M:%S", localtime( &now ) );
    std::string time = cptime;

    // Generate time and date stamp for report filename YYMMDD_HHMMSS
    strftime( cptime, 50, "%y%m%d", localtime( &now ) );
    std::string dateForFile = cptime;
    strftime( cptime, 50, "%H%M%S", localtime( &now ) );
    std::string timeForFile = cptime;

    // Open output file stream for ( YYMMDD_HHMMSS_unitTestReport.txt ) unit
    // test report.
    std::ofstream unitTestReportOutputFile;
    std::string reportFileName = dateForFile + "_" + timeForFile
            + "_unitTestReport.txt";
    unitTestReportOutputFile.open( reportFileName.c_str( ) );
    unitTestReportOutputFile << "Tudat Unit Test Report" << std::endl;
    unitTestReportOutputFile << date << ", " << time
                             << std::endl << std::endl;

    // Write unit test results for Mathematics.
    unitTestReportOutputFile << "Mathematics" << endl;
    unitTestReportOutputFile << testBasicMathematicsFunctions
                             << "\tBasic Mathematics Functions" << endl;
    unitTestReportOutputFile << testCubicSplineInterpolation
                             << "\tCubic Spline Interpolation" << endl;
    unitTestReportOutputFile << testEulerIntegrator
                             << "\tEuler Integrator" << endl;
    unitTestReportOutputFile << testExponentialRandomNumberGenerator
                             << "\tExponential Random Number Generator"
                             << endl;
    unitTestReportOutputFile << testLawgsSurfaceGeometry
                             << "\tLawgs Surface Geometry" << endl;
    unitTestReportOutputFile << testNewtonRaphson
                             << "\tNewton Raphson Root Finder" << endl;
    unitTestReportOutputFile << testNormalRandomNumberGenerator
                             << "\tNormal Random Number Generator" << endl;
    unitTestReportOutputFile << testRungeKutta4thOrderFixedStepsizeIntegrator
                             << "\tRunge-Kutta 4th-Order Fixed-Stepsize "
                             << "Integrator" << endl;
    unitTestReportOutputFile << testSimpleLinearRegression
                             << "\tSimple Linear Regression" << endl;
    unitTestReportOutputFile << testUniformRandomNumberGenerator
                             << "\tUniform Random Number Generator" << endl;
    unitTestReportOutputFile << testUnitConversions
                             << "\tUnit Conversions" << endl;

    // Empty line.
    unitTestReportOutputFile << endl;

    // Write unit test results for Astrodynamics.
    unitTestReportOutputFile << "Astrodynamics" << endl;
    unitTestReportOutputFile << testAerodynamicForcesAndMoments
                             << "\tAerodynamic Forces and Moments" << endl;
    unitTestReportOutputFile << testAerodynamicsNamespace
                             << "\tAerodynamics Namespace" << endl;
    unitTestReportOutputFile << testApproximatePlanetPositions
                             << "\tApproximate Planet Positions" << endl;
    unitTestReportOutputFile << testCartesianElements
                             << "\tCartesian Elements" << endl;
    unitTestReportOutputFile << testCoefficientGenerator
                             << "\tAerodynamics Coefficient Generator" << endl;
    unitTestReportOutputFile << testDeepSpaceManeuver
                             << "\tDeep Space Maneuver" << endl;
    unitTestReportOutputFile << testEscapeAndCapture
                             << "\tEscape and Capture" << endl;
    unitTestReportOutputFile << testExponentialAtmosphere
                             << "\tExponential Atmosphere" << endl;
    unitTestReportOutputFile << testGravityAssist
                             << "\tGravity Assist" << endl;
    unitTestReportOutputFile << testKeplerianElements
                             << "\tKeplerian Elements" << endl;
    unitTestReportOutputFile << testKeplerPropagator
                             << "\tKepler Propagator" << endl;
    unitTestReportOutputFile << testLambertTargeter
                             << "\tLambert Targeter" << endl;
    unitTestReportOutputFile << testLibrationPointLocations
                             << "\tLibration Point Locations" << endl;
    unitTestReportOutputFile << testNumericalPropagator
                             << "\tNumerical Propagator" << endl;
    unitTestReportOutputFile << testOrbitalElementConversions
                             << "\tOrbital Element Conversions" << endl;
    unitTestReportOutputFile << testPhysicalConstants
                             << "\tPhysical Constants" << endl;
    unitTestReportOutputFile << testReferenceFrameTransformations
                             << "\tReference Frame Transformations" << endl;
    unitTestReportOutputFile << testSphericalHarmonicsGravityField
                             << "\tSpherical Harmonics Gravity Field" << endl;
    unitTestReportOutputFile << testTabulatedAtmosphere
                             << "\tTabulated Atmosphere" << endl;

    // Empty line.
    unitTestReportOutputFile << endl;

    // Write unit test results for Input.
    unitTestReportOutputFile << "Input" << endl;
    unitTestReportOutputFile << testTextFileReader
                             << "\tText File Reader" << endl;
    unitTestReportOutputFile << testTwoLineElementsTextFileReader
                             << "\tTwo-Line Elements Text File Reader" << endl;

    // Return success variable.
    return isErroneous;
}

// End of file.
