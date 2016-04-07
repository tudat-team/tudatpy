/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      102511    D. Dirkx          First version of file.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      120328    D. Dirkx          Updated code to use shared_ptrs instead of raw pointers.
 *      120605    J. Vandamme       Boostified unit test.
 *      120825    A. Ronse          Corrected capsule test to check case with angle of attack = 0.
 *                                  Extended spherical test to check for correct settings.
 *      121108    A. Ronse          Updated unit test to new generator architecture, corrected
 *                                  Apollo expected values and tolerances.
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *        Program, Volume II - Program Formulation, Douglas Aircraft Aircraft Company, 1973.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/array.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"

namespace tudat
{
namespace unit_tests
{

using basic_mathematics::Vector6d;
using mathematical_constants::PI;
using std::vector;

using namespace aerodynamics;

BOOST_AUTO_TEST_SUITE( test_aerodynamic_coefficient_generator )

//! Test coefficient generator.
BOOST_AUTO_TEST_CASE( testAerodynamicCoefficientGenerator )
{


    // Set units of coefficients
    const double expectedValueOfForceCoefficient = 1.0;

    // Tolerance in absolute units.
    const double toleranceForceCoefficient = 1.0e-2;
    const double toleranceAerodynamicCoefficients3 = 1.0e-4;
    const double toleranceAerodynamicCoefficients4 = 1.0e-2;
    const double toleranceAerodynamicCoefficients5 = 1.0e-4;

    // Create test sphere.
    boost::shared_ptr< geometric_shapes::SphereSegment > sphere
            = boost::make_shared< geometric_shapes::SphereSegment >( 1.0 );

    // Set vehicle in analysis with 10,000 panels.
    vector< int > numberOfLines( 1, 31 );
    vector< int > numberOfPoints( 1, 31 );
    vector< bool > invertOrder( 1, false );

    // Create analysis object.
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );

    independentVariableDataPoints[ 0 ] =
            getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    independentVariableDataPoints[ 1 ] =
            getDefaultHypersonicLocalInclinationAngleOfAttackPoints( );
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );


    // Set methods to use for aerodynamic analysis.
    std::vector< std::vector< int > > analysisMethod( 2, std::vector< int >( 1, 0 ));
    analysisMethod[ 1 ][ 0 ] = 1;

    // Generate sphere database of aerodynamic coefficients.
    boost::shared_ptr< HypersonicLocalInclinationAnalysis > coefficientInterface =
            boost::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, sphere,
                numberOfLines, numberOfPoints,
                invertOrder, analysisMethod, PI, 1.0,
                Eigen::Vector3d::Zero( ) );

    // Test basic properties of coefficient generator
    BOOST_CHECK_EQUAL(
                coefficientInterface->getIndependentVariableNames( ).size( ), 3 );
    BOOST_CHECK_EQUAL( coefficientInterface->getIndependentVariableName( 0 ),
                       mach_number_dependent );
    BOOST_CHECK_EQUAL( coefficientInterface->getIndependentVariableName( 1 ),
                       angle_of_attack_dependent );
    BOOST_CHECK_EQUAL( coefficientInterface->getIndependentVariableName( 2 ),
                       angle_of_sideslip_dependent );

    bool isVariableIndexTooHigh = 0;
    try
    {
        coefficientInterface->getIndependentVariableName( 3 );
    }
    catch ( std::runtime_error )
    {
        isVariableIndexTooHigh = 1;
    }
    BOOST_CHECK_EQUAL( isVariableIndexTooHigh, 1 );


    // Allocate memory for independent variables to pass to analysis for retrieval.
    boost::array< int, 3 > independentVariables;
    independentVariables[ 0 ] = 0;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;
    std::vector< double > independentVariablesVector( 3 );
    std::vector< double > interpolatingIndependentVariablesVector( 3 );

    // Declare local test variables.
    Vector6d aerodynamicCoefficients_ = Vector6d::Zero( );
    double forceCoefficient_;

    // Iterate over all angles of attack to verify sphere coefficients. Total force coefficient
    // should be one; all moment coefficients should be zero.
    // The functionality is tested directly from the generator, as well as from the
    // coefficient interface, both interpolated at the nodes, and halfway between the nodes.
    for ( int i = 0; i < coefficientInterface->getNumberOfValuesOfIndependentVariable( 0 ); i++ )
    {
        independentVariables[ 0 ] = i;
        independentVariablesVector[ 0 ] = coefficientInterface->getIndependentVariablePoint( 0, i );
        if( i < coefficientInterface->getNumberOfValuesOfIndependentVariable( 0 ) - 1 )
        {
            interpolatingIndependentVariablesVector[ 0 ] =
                    coefficientInterface->getIndependentVariablePoint( 0, i ) + 0.5 * (
                        coefficientInterface->getIndependentVariablePoint( 0, i + 1 ) -
                        coefficientInterface->getIndependentVariablePoint( 0, i ) );
        }


        for ( int j = 0; j <
              coefficientInterface->getNumberOfValuesOfIndependentVariable( 1 ); j++ )
        {
            independentVariables[ 1 ] = j;
            independentVariablesVector[ 1 ] =
                    coefficientInterface->getIndependentVariablePoint( 1, j );
            if( j < coefficientInterface->getNumberOfValuesOfIndependentVariable( 1 ) - 1 )
            {
                interpolatingIndependentVariablesVector[ 1 ] =
                        coefficientInterface->getIndependentVariablePoint( 1, j ) + 0.5 * (
                            coefficientInterface->getIndependentVariablePoint( 1, j + 1 ) -
                            coefficientInterface->getIndependentVariablePoint( 1, j ) );
            }

            for ( int k = 0; k <
                  coefficientInterface->getNumberOfValuesOfIndependentVariable( 2 ); k++ )
            {
                independentVariables[ 2 ] = k;
                independentVariablesVector[ 2 ] =
                        coefficientInterface->getIndependentVariablePoint( 2, k );
                if( k < coefficientInterface->getNumberOfValuesOfIndependentVariable( 2 ) - 1 )
                {
                    interpolatingIndependentVariablesVector[ 2 ] =
                            coefficientInterface->getIndependentVariablePoint( 2, k ) + 0.5 * (
                                coefficientInterface->getIndependentVariablePoint( 2, k + 1 ) -
                                coefficientInterface->getIndependentVariablePoint( 2, k ) );
                }

                // Retrieve aerodynamic coefficients.
                aerodynamicCoefficients_ =
                        coefficientInterface->getAerodynamicCoefficientsDataPoint(
                            independentVariables );
                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

                // Test if the computed force coefficient corresponds to the expected value
                // within the specified tolerance.
                BOOST_CHECK_CLOSE_FRACTION( forceCoefficient_,
                                            expectedValueOfForceCoefficient,
                                            toleranceForceCoefficient );

                // Test if the computed moment coefficients correspond to the expected value (0.0)
                // within the specified tolerance.
                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 3 ),
                                   toleranceAerodynamicCoefficients3 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 4 ),
                                   toleranceAerodynamicCoefficients4 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 5 ),
                                   toleranceAerodynamicCoefficients5 );

                // Retrieve aerodynamic coefficients from coefficient interface.
                coefficientInterface->updateCurrentCoefficients( independentVariablesVector );

                aerodynamicCoefficients_ =
                        coefficientInterface->getCurrentAerodynamicCoefficients( );
                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

                // Test if the computed force coefficient corresponds to the expected value
                // within the specified tolerance.
                BOOST_CHECK_CLOSE_FRACTION( forceCoefficient_,
                                            expectedValueOfForceCoefficient,
                                            toleranceForceCoefficient );

                // Test if the computed moment coefficients correspond to the expected value (0.0)
                // within the specified tolerance.
                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 3 ),
                                   toleranceAerodynamicCoefficients3 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 4 ),
                                   toleranceAerodynamicCoefficients4 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 5 ),
                                   toleranceAerodynamicCoefficients5 );

                // Retrieve aerodynamic coefficients from coefficient interface.
                coefficientInterface->updateCurrentCoefficients(
                            interpolatingIndependentVariablesVector );

                aerodynamicCoefficients_ =
                        coefficientInterface->getCurrentAerodynamicCoefficients( );
                forceCoefficient_ = ( aerodynamicCoefficients_.head( 3 ) ).norm( );

                // Test if the computed force coefficient corresponds to the expected value
                // within the specified tolerance.
                BOOST_CHECK_CLOSE_FRACTION( forceCoefficient_,
                                            expectedValueOfForceCoefficient,
                                            toleranceForceCoefficient );

                // Test if the computed moment coefficients correspond to the expected value (0.0)
                // within the specified tolerance.
                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 3 ),
                                   toleranceAerodynamicCoefficients3 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 4 ),
                                   toleranceAerodynamicCoefficients4 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_( 5 ),
                                   toleranceAerodynamicCoefficients5 );
            }
        }
    }
}

boost::shared_ptr< HypersonicLocalInclinationAnalysis > getApolloCoefficientInterface( )
{

    // Create test capsule.
    boost::shared_ptr< geometric_shapes::Capsule > capsule
            = boost::make_shared< geometric_shapes::Capsule >(
                4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

    vector< int > numberOfLines( 4 );
    vector< int > numberOfPoints( 4 );
    vector< bool > invertOrders( 4 );

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
    momentReference( 2 ) = -0.1369;

    std::vector< std::vector< double > > independentVariableDataPoints( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints( 7 );

    for ( int i = 0; i < 7; i++ )
    {
        angleOfAttackPoints[ i ] = static_cast< double >( i - 6 ) * 5.0 * PI / 180.0;
    }

    independentVariableDataPoints[ 1 ] = angleOfAttackPoints;
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );
    std::vector< std::vector< int > > selectedMethods( 2, std::vector< int >( 4 ));

    selectedMethods[ 0 ][ 0 ] = 1;
    selectedMethods[ 0 ][ 1 ] = 5;
    selectedMethods[ 0 ][ 2 ] = 5;
    selectedMethods[ 0 ][ 3 ] = 1;
    selectedMethods[ 1 ][ 0 ] = 6;
    selectedMethods[ 1 ][ 1 ] = 3;
    selectedMethods[ 1 ][ 2 ] = 3;
    selectedMethods[ 1 ][ 3 ] = 3;

    // Create analysis object and capsule database.
    return boost::make_shared< HypersonicLocalInclinationAnalysis >(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * pow( capsule->getMiddleRadius( ), 2.0 ),
                3.9116, momentReference );
}
//! Apollo capsule test case.
BOOST_AUTO_TEST_CASE( testApolloCapsule )
{
    // Set units of coefficients.
    const double expectedValueOfAerodynamicCoefficients0 = -1.51;
    const double expectedValueOfAerodynamicCoefficients4 = -0.052;

    // Tolerance in absolute units.
    const double toleranceAerodynamicCoefficients0 = 0.1;
    const double toleranceAerodynamicCoefficients1 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients2 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients3 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients4 = 0.05;
    const double toleranceAerodynamicCoefficients5 = std::numeric_limits< double >::epsilon( );

    // Create aerodynamic coefficients.
    boost::shared_ptr< HypersonicLocalInclinationAnalysis > coefficientInterface = getApolloCoefficientInterface( );

    // Retrieve coefficients at zero angle of attack for comparison.
    boost::array< int, 3 > independentVariables;

    independentVariables[ 0 ] =
            coefficientInterface->getNumberOfValuesOfIndependentVariable( 0 ) - 1;
    independentVariables[ 1 ] = 6;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::VectorXd aerodynamicCoefficients_;
    aerodynamicCoefficients_ = coefficientInterface->getAerodynamicCoefficientsDataPoint(
                independentVariables );

    // Compare values to database values.
    BOOST_CHECK_SMALL(
                aerodynamicCoefficients_( 0 ) - expectedValueOfAerodynamicCoefficients0,
                toleranceAerodynamicCoefficients0 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 1 ),
                       toleranceAerodynamicCoefficients1 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 2 ),
                       toleranceAerodynamicCoefficients2 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 3 ),
                       toleranceAerodynamicCoefficients3 );

    BOOST_CHECK_SMALL(
                aerodynamicCoefficients_( 4 ) - expectedValueOfAerodynamicCoefficients4,
                toleranceAerodynamicCoefficients4 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 5 ),
                       toleranceAerodynamicCoefficients5 );


    std::vector< double > independentVariablesVector( 3 );
    independentVariablesVector[ 0 ] = coefficientInterface->getIndependentVariablePoint(
                0, coefficientInterface->getNumberOfValuesOfIndependentVariable( 0 ) - 1 );
    independentVariablesVector[ 1 ] = coefficientInterface->getIndependentVariablePoint(
                1, 6 );
    independentVariablesVector[ 2 ] = coefficientInterface->getIndependentVariablePoint(
                2, 0 );

    coefficientInterface->updateCurrentCoefficients( independentVariablesVector );
    aerodynamicCoefficients_ = coefficientInterface->getCurrentAerodynamicCoefficients( );

    // Compare values to database values.
    BOOST_CHECK_SMALL(
                aerodynamicCoefficients_( 0 ) - expectedValueOfAerodynamicCoefficients0,
                toleranceAerodynamicCoefficients0 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 1 ),
                       toleranceAerodynamicCoefficients1 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 2 ),
                       toleranceAerodynamicCoefficients2 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 3 ),
                       toleranceAerodynamicCoefficients3 );

    BOOST_CHECK_SMALL(
                aerodynamicCoefficients_( 4 ) - expectedValueOfAerodynamicCoefficients4,
                toleranceAerodynamicCoefficients4 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_( 5 ),
                       toleranceAerodynamicCoefficients5 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
