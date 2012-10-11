/*    Copyright (c) 2010-2012, Delft University of Technology
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
 */

#define BOOST_TEST_MAIN

#include <boost/array.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_coefficient_generator )

//! Test coefficient generator.
BOOST_AUTO_TEST_CASE( testAerodynamicCoefficientGenerator )
{
    using tudat::mathematics::PI;
    using std::vector;
    using namespace aerodynamics;

    // Set units of coefficients
    const double expectedValueOfForceCoefficient = 1.0;

    // Tolerance in absolute units.
    const double toleranceForceCoefficient = 1.0e-2;
    const double toleranceAerodynamicCoefficients3 = 1.0e-4;
    const double toleranceAerodynamicCoefficients4 = 1.0e-2;
    const double toleranceAerodynamicCoefficients5 = 1.0e-4;

    // Create test sphere.
    boost::shared_ptr< mathematics::geometric_shapes::SphereSegment > sphere
            = boost::make_shared< mathematics::geometric_shapes::SphereSegment >( 1.0 );

    // Set vehicle in analysis with 10,000 panels.
    vector< int > numberOfLines;
    vector< int > numberOfPoints;
    vector< bool > invertOrder;
    numberOfLines.resize( 1 );
    numberOfPoints.resize( 1 );
    invertOrder.resize( 1 );
    numberOfLines[ 0 ] = 31;
    numberOfPoints[ 0 ] = 31;
    invertOrder[ 0 ] = 0;

    // Create analysis object.
    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );

    independentVariableDataPoints[ 0 ] =
            getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    independentVariableDataPoints[ 1 ] =
            getDefaultHypersonicLocalInclinationAngleOfAttackPoints( );
    independentVariableDataPoints[ 2 ] =
            getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

    std::vector< std::vector< int > > analysisMethod;
    analysisMethod.resize( 2 );
    analysisMethod[ 0 ].resize( 1 );
    analysisMethod[ 1 ].resize( 1 );
    analysisMethod[ 0 ][ 0 ] = 0;
    analysisMethod[ 1 ][ 0 ] = 1;

    HypersonicLocalInclinationAnalysis analysis( independentVariableDataPoints, sphere,
                                                 numberOfLines, numberOfPoints,
                                                 invertOrder, analysisMethod, PI, 1.0,
                                                 Eigen::Vector3d::Zero( ) );

    // Generate sphere database.
    analysis.generateCoefficients( );

    // Allocate memory for independent variables to pass to analysis for retrieval.
    boost::array< int, 3 > independentVariables;
    independentVariables[ 0 ] = 0;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::Matrix< double, 6, 1 > aerodynamicCoefficients_ =
            Eigen::Matrix< double, 6, 1 > ::Zero( );
    double forceCoefficient_;

    // Iterate over all angles of attack to verify sphere coefficients.
    for ( int i = 0; i < analysis.getNumberOfValuesOfIndependentVariable( 0 ); i++ )
    {
        independentVariables[ 0 ] = i;

        for ( int j = 0; j < analysis.getNumberOfValuesOfIndependentVariable( 1 ); j++ )
        {
            independentVariables[ 1 ] = j;

            for ( int k = 0; k < analysis.getNumberOfValuesOfIndependentVariable( 2 ); k++ )
            {
                independentVariables[ 2 ] = k;

                // Retrieve aerodynamic coefficients.
                aerodynamicCoefficients_ = analysis.getAerodynamicCoefficients(
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
            }
        }
    }
}

//! Apollo capsule test case.
BOOST_AUTO_TEST_CASE( testApolloCapsule )
{
    using tudat::mathematics::PI;
    using std::vector;
    using namespace aerodynamics;

    // Set units of coefficients.
    const double expectedValueOfAerodynamicCoefficients0 = -1.51;
    const double expectedValueOfAerodynamicCoefficients4 = -0.052;

    // Tolerance in absolute units.
    const double toleranceAerodynamicCoefficients0 = 0.05;
    const double toleranceAerodynamicCoefficients1 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients2 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients3 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients4 = 0.975;
    const double toleranceAerodynamicCoefficients5 = std::numeric_limits< double >::epsilon( );

    // Create test capsule.
    boost::shared_ptr< mathematics::geometric_shapes::Capsule > capsule
            = boost::make_shared< mathematics::geometric_shapes::Capsule >(
                    4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );

    vector< int > numberOfLines;
    vector< int > numberOfPoints;
    vector< bool > invertOrders;
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
    momentReference( 2 ) = -0.1369;

    std::vector< std::vector< double > > independentVariableDataPoints;
    independentVariableDataPoints.resize( 3 );
    independentVariableDataPoints[ 0 ] = getDefaultHypersonicLocalInclinationMachPoints( "Full" );
    std::vector< double > angleOfAttackPoints;
    angleOfAttackPoints.resize( 7 );

    for ( int i = 0; i < 7; i++ )
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

    // Create analysis object.
    HypersonicLocalInclinationAnalysis analysis = HypersonicLocalInclinationAnalysis(
                independentVariableDataPoints, capsule, numberOfLines, numberOfPoints,
                invertOrders, selectedMethods, PI * pow( capsule->getMiddleRadius( ), 2.0 ),
                3.9116, momentReference );

    // Generate capsule database.
    analysis.generateCoefficients( );

    // Retrieve coefficients at zero angle of attack for comparison.
    boost::array< int, 3 > independentVariables;

    independentVariables[ 0 ] = analysis.getNumberOfValuesOfIndependentVariable( 0 ) - 1;
    independentVariables[ 1 ] = 6;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::VectorXd aerodynamicCoefficients_;
    aerodynamicCoefficients_ = analysis.getAerodynamicCoefficients( independentVariables );

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
