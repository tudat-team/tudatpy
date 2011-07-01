/*! \file unitTestMain.cpp
 *    Source file for the main unit test executable that runs all unit
 *    tests and generates a unit_test_report.txt report.
 *    A 0 corresponds to success in the report.
 *
 *    Path              : /
 *    Version           : 4
 *    Check status      : Unchecked
 *
 *    Author            : B. Römgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 January, 2011
 *    Last modified     : 29 March, 2011
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
 */

// Include statements.
#include <fstream>
#include <ctime>

// Mathematics unit test includes
#include "unitTestBasicMathematicsFunctions.h"
#include "unitTestRandomNumberGenerator.h"
#include "unitTestUnitConversions.h"
#include "unitTestNewtonRaphson.h"
#include "unitTestLawgsSurfaceGeometry.h"

// Astrodynamics unit test includes
#include "unitTestPhysicalConstants.h"
#include "unitTestCartesianElements.h"
#include "unitTestKeplerianElements.h"
#include "unitTestOrbitalElementConversions.h"
#include "unitTestSphericalHarmonicsGravityField.h"
#include "unitTestNumericalPropagator.h"
#include "unitTestKeplerPropagator.h"
#include "unitTestAerodynamicsNamespace.h"
#include "unitTestCoefficientGenerator.h"
#include "unitTestEscapeAndCapture.h"
#include "unitTestGravityAssist.h"
#include "unitTestLambertTargeter.h"
#include "unitTestReferenceFrameTransformations.h"

// Input unit test includes
#include "unitTestTextFileReader.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Execute all unit tests.
/*!
 * Executes all unit tests.
 */
int main()
{
    // Return value, 0 on success, 1 on failure of one or more unit tests.
    int success = 0;

    // Run all unit tests.
    // Unit test functions return false ( 0 ) for success.

    // Run Mathematics unit tests.

    // testBasicMathematicsFunctions: Tests the following functions,
    // linear interpolation, absolute value, raise to integer power
    // modulo, spherical to Cartesian, cylindrical to Cartesian.
    bool testBasicMathematicsFunctions =
            unit_tests::testBasicMathematicsFunctions( );

    // testUnitConversions: Tests conversions that are defined
    // in unitConversions.h.
    bool testUnitConversions = unit_tests::testUnitConversions( );

    // testRandomNumberGenerator: Tests the random number generator
    // defined in randomNumberGenerator.h
    bool testRandomNumberGenerator = unit_tests::testRandomNumberGenerator( );

    // testNewtonRaphson: Tests the Newton-Raphson root-finder.
    bool testNewtonRaphson = unit_tests::testNewtonRaphsonMethod( );

    // testLawgsSurfaceGeometry: Tests a Lawgs mesh of a sphere.
    bool testLawgsSurfaceGeometry =
            unit_tests::testLawgsSurfaceGeometry( );

    // Run Astrodynamics unit tests.

    // testPhysicalConstants: Tests the physical constants that are defined
    // in physicalConstants.h.
    bool testPhysicalConstants = unit_tests::testPhysicalConstants( );

    // testCartesianElements: Tests the different set and get functions
    // of CartesianElements.
    bool testCartesianElements =
            unit_tests::testCartesianElements( );

    // testKeplerianElements: Tests the different set and get functions
    // of KeplerianElements.
    bool testKeplerianElements =
            unit_tests::testKeplerianElements( );

    // testOrbitalElementConversions: Tests the code for elliptical, parabolic,
    // hyperbolic and circular orbits. It also tests the conversion from
    // Cartesian to Keplerian and Keplerian to Cartesian.
    bool testOrbitalElementConversions =
            unit_tests::testOrbitalElementConversions( );

    // testSphericalHarmonicsGravityField: Tests the implementation of the
    // spherical harmonics gravity field class.
    bool testSphericalHarmonicsGravityField =
            unit_tests::testSphericalHarmonicsGravityField( );

    // testNumericalPropagator: Tests the numerical propagator.
    bool testNumericalPropagator = unit_tests::testNumericalPropagator( );

    // testKeplerPropagator: Tests the Kepler propagator. For this test, the
    // ROOT_PATH variable must be set in basicFunctions.h to the root
    // directory for Tudat.
    bool testKeplerPropagator = unit_tests::testKeplerPropagator( );

    // testAerodynamicsNamespace: Tests the following functions:
    // static pressure ratio, stagnation pressure coefficient,
    // modified Newtonian pressure coefficient, Tangent Cone pressure coef.,
    // Mach base pressure coefficient, Tangent Wedge pressure coefficient,
    // freestream Prandtl-Meyer function, vacuum pressure coefficient,
    // shock pressure ratio, shock density ratio, shock temperature ratio
    // shock total pressure ratio, shock wave total pressure ratio.
    bool testAerodynamicsNamespace =
            unit_tests::testAerodynamicsNameSpace( );

    // testCoefficientGenerator: Tests aerothermodynamic coefficients for
    // a sphere and the apollo capsule.
    bool testCoefficientGenerator =
            unit_tests::testCoefficientGenerator( );

    // testEscapeAndCapture: Tests the computed delta-V of escape and capture.
    bool testEscapeAndCapture =
            unit_tests::testEscapeAndCapture( );

    // testGravityAssist: Tests delta-V for a gravity assist for the
    // case of equal hyperbolic excess velocities.
    bool testGravityAssist =
            unit_tests::testGravityAssist( );

    // testLambertTargeter: Tests the values of semi major axis,
    // radial speed at departure, radial speed at arrival,
    // transverse speed at departure, transverse speed at arrival for
    // elliptical and hypterbolic test case.
    bool testLambertTargeter =
            unit_tests::testLambertTargeter( );

    // TestFrameTransformations: Tests the reference frame transformations:
    // Rotating planetocentric (R) to Inertial (I) frame transformation,
    // Inertial (I) to Rotating planetocentric (R) frame transformation,
    // Aerodynamic (airspeed based) (AA) to Body (B) frame transformation.
    bool testReferenceFrameTransformations
            = unit_tests::testReferenceFrameTransformations( );

    // Run Input unit tests.

    // testTextFileReader: Tests reading in a text file using the readAndStore,
    // skipLines, and skipLinesStartingWithCharacter functions.
    bool testTextFileReader =
            unit_tests::testTextFileReader( );

    // Check if all tests were successful and print cerr message stating
    // which test failed, if any.
    if ( testBasicMathematicsFunctions )
    {
        cerr << "testBasicMathematicsFunctions failed!" << endl;
        success = 1;
    }

    if ( testUnitConversions )
    {
        cerr << "testUnitConversions failed!" << endl;
        success = 1;
    }

    if ( testRandomNumberGenerator )
    {
        cerr << "testRandomNumberGenerator failed!" << endl;
        success = 1;
    }

    if ( testNewtonRaphson )
    {
        cerr << "testNewtonRaphson failed!" << endl;
        success = 1;
    }

    if ( testLawgsSurfaceGeometry )
    {
        cerr << "testLawgsSurfaceGeometry failed!" << endl;
        success = 1;
    }

    if ( testPhysicalConstants )
    {
        cerr << "testPhysicalConstants failed!" << endl;
        success = 1;
    }

    if ( testCartesianElements )
    {
        cerr << "testCartesianElements failed!" << endl;
        success = 1;
    }

    if ( testKeplerianElements )
    {
        cerr << "testKeplerianElements failed!" << endl;
        success = 1;
    }

    if ( testOrbitalElementConversions )
    {
        cerr << "testOrbitalElementConversions failed!" << endl;
        success = 1;
    }

    if ( testSphericalHarmonicsGravityField )
    {
        cerr << "testSphericalHarmonicsGravityField failed!" << endl;
        success = 1;
    }

    if ( testNumericalPropagator )
    {
        cerr << "testNumericalPropagator failed!" << endl;
        success = 1;
    }

    if ( testKeplerPropagator )
    {
        cerr << "testKeplerPropagator failed!" << endl;
        success = 1;
    }

    if ( testAerodynamicsNamespace )
    {
        cerr << "testAerodynamicsNamespace failed!" << endl;
        success = 1;
    }

    if ( testCoefficientGenerator )
    {
        cerr << "testCoefficientGenerator failed!" << endl;
        success = 1;
    }

    if ( testEscapeAndCapture )
    {
        cerr << "testEscapeAndCapture failed!" << endl;
        success = 1;
    }

    if ( testGravityAssist )
    {
        cerr << "testGravityAssist failed!" << endl;
        success = 1;
    }

    if ( testLambertTargeter )
    {
        cerr << "testLambertTargeter failed!" << endl;
        success = 1;
    }

    if ( testTextFileReader )
    {
        cerr << "testTextFileReader failed!" << endl;
        success = 1;
    }

    if ( testReferenceFrameTransformations )
    {
        cerr << "referenceFrameTransformations failed" << endl;
        success = 1;
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
    unitTestReportOutputFile << testUnitConversions
            << "\tUnit Conversions" << endl;
    unitTestReportOutputFile << testRandomNumberGenerator
            << "\tRandom Number Generator" << endl;
    unitTestReportOutputFile << testNewtonRaphson
            << "\tNewton Raphson Root Finder" << endl;
    unitTestReportOutputFile << testLawgsSurfaceGeometry
            << "\tLawgs Surface Geometry" << endl;




    // Empty line.
    unitTestReportOutputFile << endl;

    // Write unit test results for Astrodynamics.
    unitTestReportOutputFile << "Astrodynamics" << endl;
    unitTestReportOutputFile << testPhysicalConstants
            << "\tPhysical Constants" << endl;
    unitTestReportOutputFile << testCartesianElements
            << "\tCartesian Elements" << endl;
    unitTestReportOutputFile << testKeplerianElements
            << "\tKeplerian Elements" << endl;
    unitTestReportOutputFile << testOrbitalElementConversions
            << "\tOrbital Element Conversions" << endl;
    unitTestReportOutputFile << testSphericalHarmonicsGravityField
            << "\tSpherical Harmonics Gravity Field" << endl;
    unitTestReportOutputFile << testNumericalPropagator
            << "\tNumerical Propagator" << endl;
    unitTestReportOutputFile << testKeplerPropagator
            << "\tKepler Propagator" << endl;
    unitTestReportOutputFile << testAerodynamicsNamespace
            << "\tAerodynamics Namespace" << endl;
    unitTestReportOutputFile << testCoefficientGenerator
            << "\tAerodynamics Coefficient Generator" << endl;
    unitTestReportOutputFile << testEscapeAndCapture
            << "\tEscape and Capture" << endl;
    unitTestReportOutputFile << testGravityAssist
            << "\tGravity Assist" << endl;
    unitTestReportOutputFile << testLambertTargeter
            << "\tLambert Targeter" << endl;
    unitTestReportOutputFile << testReferenceFrameTransformations
            << "\tReference Frame Transformations" << endl;

    // Empty line.
    unitTestReportOutputFile << endl;

    // Write unit test results for Input.
    unitTestReportOutputFile << "Astrodynamics" << endl;
    unitTestReportOutputFile << testTextFileReader
            << "\tPhysical Constants" << endl;

    // Return success variable.
    return success;
}

// End of file.
