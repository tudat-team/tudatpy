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
 *      110208    D. Dirkx          First version of file.
 *      110210    L. Adulkadir      Code check.
 *      110211    K. Kumar          Corrected Doxygen errors; corrected layout errors; corrected
 *                                  double precision; updated function-naming.
 *      120605    J. Vandamme       Boostified unit test.
 *      120912    K. Kumar          Removed erroneous using-statements; added missing <cmath>
 *                                  include.
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill, 2001.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition,
 *          AIAA Education Series, 2006.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamics_namespace )

//! Test aerodynamic namespace pressure functions.
BOOST_AUTO_TEST_CASE( testAerodynamicNamespacePressureFunctions )
{
    using namespace tudat;
    using namespace aerodynamics;
    using mathematical_constants::PI;

    // Set default test conditions.
    const double machNumber_ = 12.0;
    const double ratioOfSpecificHeats_ = 1.4;

    // Test local to static pressure ratio.
    const double localToStaticPressureRatio_ = computeLocalToStaticPressureRatio(
           machNumber_, ratioOfSpecificHeats_ );

    const double expectedLocalToStaticPressureRatio_ = 1.0 / 0.1445e6;
    const double toleranceRatioCoefficient_ = 0.1445e6 * 1.0e-8 * 100.0
            / expectedLocalToStaticPressureRatio_;

    // Check if computed local-to-static pressure ratio matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( localToStaticPressureRatio_,
                                expectedLocalToStaticPressureRatio_,
                                toleranceRatioCoefficient_ );

    // Test stagnation pressure coefficient.
    const double stagnationPressureCoefficient_
            = computeStagnationPressure( machNumber_, ratioOfSpecificHeats_ );

    const double expectedStagnationPressureCoefficient_ = 1.83402;
    const double toleranceStagnationPressureCoefficient_ = 1.0e-5 * 100.0
            / expectedStagnationPressureCoefficient_;

    // Check if computed stagnation pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( stagnationPressureCoefficient_,
                                expectedStagnationPressureCoefficient_,
                                toleranceStagnationPressureCoefficient_ );

    // Test modified Newtonian pressure coefficient.
    const double newtonianPressureCoefficient_ = computeModifiedNewtonianPressureCoefficient(
                PI / 2.0, stagnationPressureCoefficient_ );

    const double expectedNewtonianPressureCoefficient_ = stagnationPressureCoefficient_;
    const double toleranceNewtonianPressureCoefficient_ = 1.0e-15 * 100.0
            / expectedNewtonianPressureCoefficient_;

    // Check if computed Newtonian pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( newtonianPressureCoefficient_,
                                expectedNewtonianPressureCoefficient_,
                                toleranceNewtonianPressureCoefficient_ );

    // Test empirical tangent cone pressure coefficient.
    const double empiricalTangentConePressureCoefficient_
            = computeEmpiricalTangentConePressureCoefficient( PI / 2.0, machNumber_ );

    const double expectedEmpiricalTangentConePressureCoefficient_ = 2.08961;
    const double toleranceEmpiricalTangentConePressureCoefficient_ = 1.0e-5 * 100.0
            / expectedEmpiricalTangentConePressureCoefficient_;

    // Check if computed empirical tangent cone pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( empiricalTangentConePressureCoefficient_,
                                expectedEmpiricalTangentConePressureCoefficient_,
                                toleranceEmpiricalTangentConePressureCoefficient_ );

    // Test high Mach base pressure coefficient.
    const double highMachBasePressure_ = computeHighMachBasePressure( machNumber_ );

    const double expectedHighMachBasePressure_ = -1.0 / ( machNumber_ * machNumber_ );
    const double toleranceHighMachBasePressure_ = std::fabs( 1.0e-15 * 100.0 / expectedHighMachBasePressure_ );

    // Check if computed high Mach base pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( highMachBasePressure_,
                                expectedHighMachBasePressure_,
                                toleranceHighMachBasePressure_ );

    // Test empirical tangent wedge pressure coefficient.
    const double empiricalTangentWedgePressureCoefficient_
            = computeEmpiricalTangentWedgePressureCoefficient( PI / 2.0, machNumber_ );

    const double expectedEmpiricalTangentWedgePressureCoefficient_ = 2.38867;
    const double toleranceEmpiricalTangentWedgePressureCoefficient_ = 1.0e-5 * 100.0
            / expectedEmpiricalTangentWedgePressureCoefficient_;

    // Check if computed tangent wedge pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( empiricalTangentWedgePressureCoefficient_,
                                expectedEmpiricalTangentWedgePressureCoefficient_,
                                toleranceEmpiricalTangentWedgePressureCoefficient_ );

    // Test free-stream Prandtl-Meyer function.
    const double freestreamPrandtlMeyerFunction_ = computePrandtlMeyerFunction(
           machNumber_, ratioOfSpecificHeats_ );

    const double expectedFreeStreamPrandtlMeyerFunction_ = 106.9 * PI / 180.0;
    const double toleranceFreeStreamPrandtlMeyerFunction_ = 1.0e-3 * 100.0
            / expectedFreeStreamPrandtlMeyerFunction_;

    // Check if computed free-stream Prandtl-Meyer function matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( freestreamPrandtlMeyerFunction_,
                                expectedFreeStreamPrandtlMeyerFunction_,
                                toleranceFreeStreamPrandtlMeyerFunction_ );

    // Test vacuum pressure coefficient function.
    const double vacuumPressureCoefficient_ = computeVacuumPressureCoefficient(
           machNumber_, ratioOfSpecificHeats_ );

    const double expectedVacuumPressureCoefficient_ = -2.0 / ( ratioOfSpecificHeats_
                                                               * machNumber_ * machNumber_ );
    const double toleranceVacuumPressureCoefficient_ = 1.0e-15 * 100.0
            / expectedVacuumPressureCoefficient_;

    // Check if computed vacuum pressure coefficient matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( vacuumPressureCoefficient_,
                                expectedVacuumPressureCoefficient_,
                                toleranceVacuumPressureCoefficient_ );
}

//! Test pressure coefficients.
BOOST_AUTO_TEST_CASE( testPressureCoefficients )
{
    using std::sin;
    using namespace tudat;
    using namespace aerodynamics;
    using mathematical_constants::PI;

    // Set default conditions.
    const double machNumber_ = 12.0;
    const double ratioOfSpecificHeats_ = 1.4;

    // Iterate over 10 angles and test pressure coefficients.
    for ( int i = 0; i < 10; i++ )
    {
         // Set angle.
         const double angle_ = static_cast< double > ( i ) * PI / 10.0;

         // Compute and compare Newtonian pressure coefficient.
         const double newtonianPressureCoefficient_
                 = computeNewtonianPressureCoefficient( angle_ );

         const double expectedNewtonianPressureCoefficient_ =  2.0 * sin( angle_ ) * sin( angle_ );
         const double toleranceNewtonianPressureCoefficient_ = 1.0e-15 * 100.0
                 / expectedNewtonianPressureCoefficient_;

         // Check if computed Newtonian pressure coefficient matches expected value.
         BOOST_CHECK_CLOSE_FRACTION( newtonianPressureCoefficient_,
                                     expectedNewtonianPressureCoefficient_,
                                     toleranceNewtonianPressureCoefficient_ );


         // Compute Prandtl-Meyer pressure coefficient and test if it is not
         // lower than vacuum pressure coefficient.
         const double freestreamPrandtlMeyerFunction_ = computePrandtlMeyerFunction(
                     machNumber_, ratioOfSpecificHeats_ );

         const double prandtlMeyerPressureCoefficient_ =
                 computePrandtlMeyerFreestreamPressureCoefficient(
                     -1.0 * angle_, machNumber_, ratioOfSpecificHeats_,
                     freestreamPrandtlMeyerFunction_ );

         const double vacuumPressureCoefficient_ = computeVacuumPressureCoefficient(
                     machNumber_, ratioOfSpecificHeats_ );

         bool isPrandtlMeyerPressureCoefficient = false;

         if ( vacuumPressureCoefficient_ < prandtlMeyerPressureCoefficient_ + 1.0e-15 )
         {
             isPrandtlMeyerPressureCoefficient = true;
         }

         // Check if Prandt-Meyer pressure coefficient is lower than vacuum pressure coefficient.
         BOOST_CHECK( isPrandtlMeyerPressureCoefficient );

         // Test shock pressure ratio.
         const double shockPressureRatio_ = computeShockPressureRatio( machNumber_,
                                                                       ratioOfSpecificHeats_ );

         const double expectedShockPressureRatio_ = 167.8;
         const double toleranceShockPressureRatio_ = 0.1 * 100.0 / expectedShockPressureRatio_;

         // Check if computed shock pressure ratio matches expected value.
         BOOST_CHECK_CLOSE_FRACTION( shockPressureRatio_,
                                     expectedShockPressureRatio_,
                                     toleranceShockPressureRatio_ );

         // Test shock density ratio.
         const double shockDensityRatio_ = computeShockDensityRatio( machNumber_,
                                                                     ratioOfSpecificHeats_ );

         const double expectedShockDensityRatio_ = 5.799;
         const double toleranceShockDensityRatio_ = 0.001 * 100 / expectedShockDensityRatio_;

         // Check if computed shock density ratio matches expected value.
         BOOST_CHECK_CLOSE_FRACTION( shockDensityRatio_,
                                     expectedShockDensityRatio_,
                                     toleranceShockDensityRatio_ );

         // Test shock temperature ratio.
         const double shockTemperatureRatio_ = computeShockTemperatureRatio(
                     machNumber_, ratioOfSpecificHeats_ );

         const double expectedShockTemperatureRatio_ = 28.94;
         const double toleranceShockTemperatureRatio_ = 0.01 * 100.0
                 / expectedShockTemperatureRatio_;

         // Check if shock temperature ratio matches expected value.
         BOOST_CHECK_CLOSE_FRACTION( shockTemperatureRatio_,
                                     expectedShockTemperatureRatio_,
                                     toleranceShockTemperatureRatio_ );

         // Test shock total pressure ratio.
         const double shockTotalPressureRatio_ = computeShockTotalPressureRatio(
                     machNumber_, ratioOfSpecificHeats_, 287.058 );

         const double expectedShockTotalPressureRatio_ = 0.001287;
         const double toleranceShockTotalPressureRatio_ = 1.0e-6 *100.0
                 / expectedShockTotalPressureRatio_;

         // Check if shock total pressure ratio matches expected value.
         BOOST_CHECK_CLOSE_FRACTION( shockTotalPressureRatio_,
                                     expectedShockTotalPressureRatio_,
                                     toleranceShockTotalPressureRatio_ );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
