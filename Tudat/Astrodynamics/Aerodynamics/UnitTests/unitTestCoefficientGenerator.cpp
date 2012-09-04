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
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *        Program, Volume II - Program Formulation, Douglas Aircraft Aircraft Company, 1973.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

using tudat::mathematics::PI;
using std::vector;

namespace tudat
{
namespace unit_tests
{

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
    boost::shared_ptr< mathematics::geometric_shapes::SphereSegment > sphere
            = boost::make_shared< mathematics::geometric_shapes::SphereSegment >( 1.0 );
    boost::shared_ptr< bodies::VehicleExternalModel > externalModel
            = boost::make_shared< bodies::VehicleExternalModel >( );
    externalModel->setVehicleGeometry( sphere );
    bodies::Vehicle vehicle;
    vehicle.setExternalModel( externalModel );

    // Create analysis object.
    aerodynamics::HypersonicLocalInclinationAnalysis analysis;

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
    analysis.setVehicle( vehicle, numberOfLines, numberOfPoints, invertOrder );

    // Set mach-points at which analysis is performed. For angle-of-attack and -sideslip, defaults
    // are used.
    analysis.setNumberOfMachPoints( 5 );
    for ( unsigned int machIterator = 0; machIterator < 5; machIterator++ )
    {
        analysis.setMachPoint( machIterator, static_cast< double >( 2 * machIterator + 10 ) );
    }

    // Set reference quantities.
    analysis.setReferenceArea( PI );
    analysis.setReferenceLength( 1.0 );
    analysis.setMomentReferencePoint( Eigen::Vector3d::Zero( ) );

    // Set pure Newtonian compression method for test purposes.
    analysis.setSelectedMethod( 0, 0, 0 );

    // Generate sphere database.
    analysis.generateDatabase( );

    // Allocate memory for independent variables to pass to analysis for retrieval.
    vector< int > independentVariables;
    independentVariables.resize( 3 );
    independentVariables[ 0 ] = 0;
    independentVariables[ 1 ] = 0;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::VectorXd aerodynamicCoefficients_;
    double forceCoefficient_;

    // Iterate over all angles of attack to verify sphere coefficients.
    for ( int i = 0; i < analysis.getNumberOfMachPoints( ); i++ )
    {
        independentVariables[ 0 ] = i;

        for ( int j = 0; j < analysis.getNumberOfAngleOfAttackPoints( ); j++ )
        {
            independentVariables[ 1 ] = j;

            for ( int k = 0; k < analysis.getNumberOfAngleOfSideslipPoints( ); k++ )
            {
                independentVariables[ 2 ] = k;

                // Test if the independent variables have correctly been set.
                BOOST_CHECK_EQUAL( analysis.getMachPoint( i ),
                                   static_cast< double >( 2 * i + 10 ) );
                BOOST_CHECK_EQUAL( analysis.getAngleOfAttackPoint( j ),
                                   static_cast< double >( j ) * 5.0 * PI / 180.0 );
                BOOST_CHECK_EQUAL( analysis.getAngleOfSideslipPoint( k ),
                                   static_cast< double >( k ) * 1.0 * PI / 180.0 );

                // Retrieve aerodynamic coefficients.
                aerodynamicCoefficients_ = analysis.getAerodynamicCoefficients(
                        independentVariables );
                forceCoefficient_ = sqrt( aerodynamicCoefficients_.x( )
                                        * aerodynamicCoefficients_.x( )
                                        + aerodynamicCoefficients_.y( )
                                        * aerodynamicCoefficients_.y( )
                                        + aerodynamicCoefficients_.z( )
                                        * aerodynamicCoefficients_.z( ) );

                // Test if the computed force coefficient corresponds to the expected value
                // within the specified tolerance.
                BOOST_CHECK_CLOSE_FRACTION( forceCoefficient_,
                                            expectedValueOfForceCoefficient,
                                            toleranceForceCoefficient );

                // Test if the computed moment coefficients correspond to the expected value (0.0)
                // within the specified tolerance.
                BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 3 ],
                                   toleranceAerodynamicCoefficients3 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 4 ],
                                   toleranceAerodynamicCoefficients4 );

                BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 5 ],
                                   toleranceAerodynamicCoefficients5 );
            }
        }
    }
}

//! Apollo capsule test case.
BOOST_AUTO_TEST_CASE( testApolloCapsule )
{
    // Set units of coefficients.
    const double expectedValueOfAerodynamicCoefficients0 = 1.51;
    const double expectedValueOfAerodynamicCoefficients4 = -0.052;

    // Tolerance in absolute units.
    const double toleranceAerodynamicCoefficients0 = 0.1 * 100.0
            / expectedValueOfAerodynamicCoefficients0;
    const double toleranceAerodynamicCoefficients1 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients2 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients3 = std::numeric_limits< double >::epsilon( );
    const double toleranceAerodynamicCoefficients4 = 0.01 * 100.0
            / expectedValueOfAerodynamicCoefficients4;
    const double toleranceAerodynamicCoefficients5 = std::numeric_limits< double >::epsilon( );

    // Create test capsule.
    boost::shared_ptr< mathematics::geometric_shapes::Capsule > capsule
            = boost::make_shared< mathematics::geometric_shapes::Capsule >(
    4.694, 1.956, 2.662, -1.0 * 33.0 * PI / 180.0, 0.196 );
    boost::shared_ptr< bodies::VehicleExternalModel > externalModel
            = boost::make_shared< bodies::VehicleExternalModel >( );
    externalModel->setVehicleGeometry( capsule );
    bodies::Vehicle vehicle;
    vehicle.setExternalModel( externalModel );

    // Create analysis object.
    aerodynamics::HypersonicLocalInclinationAnalysis analysis;

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
    invertOrders[ 0 ] = 1;
    invertOrders[ 1 ] = 1;
    invertOrders[ 2 ] = 1;
    invertOrders[ 3 ] = 1;

    analysis.setVehicle( vehicle, numberOfLines, numberOfPoints, invertOrders );

    // Set reference quantities.
    analysis.setReferenceArea( PI * pow( capsule->getMiddleRadius( ), 2.0 ) );
    analysis.setReferenceLength( 3.9116 );
    Eigen::VectorXd momentReference = Eigen::VectorXd( 3 );
    momentReference( 0 ) = 0.6624;
    momentReference( 1 ) = 0.0;
    momentReference( 2 ) = 0.1369;
    analysis.setMomentReferencePoint( momentReference );

    // Set pure Newtonian compression method for test purposes.
    analysis.setSelectedMethod( 0, 0, 0 );

    // Set angle of attack analysis points.
    analysis.setNumberOfAngleOfAttackPoints( 7 );
    int i;
    for ( i = 0; i < 7; i++ )
    {
        analysis.setAngleOfAttackPoint( i, static_cast< double >( i - 6 ) * 5.0 * PI / 180.0 );
    }

    // Generate capsule database.
    analysis.generateDatabase( );

    // Retrieve coefficients at zero angle of attack for comparison.
    vector< int > independentVariables;
    independentVariables.resize( 3 );

    independentVariables[ 0 ] = analysis.getNumberOfMachPoints( ) - 1;
    independentVariables[ 1 ] = 6;
    independentVariables[ 2 ] = 0;

    // Declare local test variables.
    Eigen::VectorXd aerodynamicCoefficients_;
    aerodynamicCoefficients_ = analysis.getAerodynamicCoefficients( independentVariables );

    // Compare values to database values.
    BOOST_CHECK_CLOSE_FRACTION( aerodynamicCoefficients_[ 0 ],
                                expectedValueOfAerodynamicCoefficients0,
                                toleranceAerodynamicCoefficients0 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 1 ],
                       toleranceAerodynamicCoefficients1 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 2 ],
                       toleranceAerodynamicCoefficients2 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 3 ],
                       toleranceAerodynamicCoefficients3 );

    BOOST_CHECK_CLOSE_FRACTION( aerodynamicCoefficients_[ 4 ],
                                expectedValueOfAerodynamicCoefficients4,
                                toleranceAerodynamicCoefficients4 );

    BOOST_CHECK_SMALL( aerodynamicCoefficients_[ 5 ],
                       toleranceAerodynamicCoefficients5 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
