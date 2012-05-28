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
 *      120203    K. Kumar          File created.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/unit_test.hpp>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_coefficients )

using mathematics::numerical_integrators::RungeKuttaCoefficients;

void checkValidityOfCoefficientSet( const RungeKuttaCoefficients::CoefficientSets& coefficientSet )
{
    // Declare coefficient set.
    RungeKuttaCoefficients coefficients;
    coefficients = coefficients.get( coefficientSet );

    // Set tolerance for unit tests.
    double tolerance = 1.0e-14;

    // Check that the sum of the b-coefficients for both the integrated order and the
    // error-checking order is one.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                Eigen::VectorXd::Constant( 2, 1.0 ),
                coefficients.bCoefficients.rowwise( ).sum( ), tolerance );

    // Check that the first c-coefficient is zero.
    BOOST_CHECK_SMALL( coefficients.cCoefficients( 0 ), tolerance );

    // Check that the c-coefficient/a-coefficient relation holds.
    for ( int i = 1; i < coefficients.cCoefficients.size( ); i++ )
    {
        if ( std::fabs( coefficients.cCoefficients( i ) ) < tolerance )
        {
            BOOST_CHECK_SMALL( coefficients.aCoefficients.row( i ).sum( ), tolerance );
        }

        else
        {
            BOOST_CHECK_CLOSE_FRACTION( coefficients.cCoefficients( i ),
                                        coefficients.aCoefficients.row( i ).sum( ), tolerance );
        }
    }
}

BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg45Coefficients )
{
    // Check validity of Runge-Kutta-Fehlberg 45 coefficients.
    checkValidityOfCoefficientSet( RungeKuttaCoefficients::rungeKuttaFehlberg45 );
}

BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg56Coefficients )
{
    // Check validity of Runge-Kutta-Fehlberg 56 coefficients.
    checkValidityOfCoefficientSet( RungeKuttaCoefficients::rungeKuttaFehlberg56 );
}

BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg78Coefficients )
{
    // Check validity of Runge-Kutta-Fehlberg 78 coefficients.
    // Note, for some reason, the rkf78 set fails the unit test when the tolerance in
    // checkValidityOfCoefficientSet() is set to lower than this value (row 8 and 9 of
    // aCoefficients matrix sum does not correspond to cCoefficient counterpart with tolerance less
    // than 1.0e-14).
    checkValidityOfCoefficientSet( RungeKuttaCoefficients::rungeKuttaFehlberg78 );
}

BOOST_AUTO_TEST_CASE( testRungeKutta87DormandAndPrinceCoefficients )
{
    // Check validity of Runge-Kutta 87 (Dormand and Prince) coefficients.
    checkValidityOfCoefficientSet( RungeKuttaCoefficients::rungeKutta87DormandPrince );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
