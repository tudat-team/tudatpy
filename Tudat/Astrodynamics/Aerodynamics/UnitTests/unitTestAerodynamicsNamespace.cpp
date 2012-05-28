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
 *      110208    D. Dirkx          First version of file.
 *      110210    L. Adulkadir      Code check.
 *      110211    K. Kumar          Corrected Doxygen errors; corrected layout errors; corrected
 *                                  double precision; updated function-naming.
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill, 2001.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition,
 *          AIAA Education Series, 2006.
 *
 */

#include <cmath>
#include <iostream>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"

using tudat::mathematics::PI;

//! Test aerodynamics namespace.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::fabs;
    using std::pow;
    using std::sin;
    using namespace tudat;
    using namespace tudat::aerodynamics;

    bool isAerodynamicsNamespaceBad = false;

    // Set default test conditions.
    double machNumber_ = 12.0;
    double ratioOfSpecificHeats_ = 1.4;

    // Test local to static pressure ratio.
    double localToStaticPressureRatio_ = computeLocalToStaticPressureRatio(
                machNumber_, ratioOfSpecificHeats_ );

    if ( fabs( localToStaticPressureRatio_ - 1.0 / ( 0.1445e6 ) ) > 1.0e-8 )
    {
        cerr << "Error in local to static pressure ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test stagnation pressure coefficient.
    double stagnationPressureCoefficient_
            = computeStagnationPressure( machNumber_, ratioOfSpecificHeats_ );

    if ( fabs( stagnationPressureCoefficient_ - 1.83402 ) > 1.0e-5 )
    {
        cerr << "Error in stagnation pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test modified Newtonian pressure coefficient.
    if ( fabs( computeModifiedNewtonianPressureCoefficient(
                   PI / 2.0, stagnationPressureCoefficient_ ) -
               stagnationPressureCoefficient_ ) > 1.0e-15 )
    {
        cerr << "Error in modified Newtonian pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test empirical Tangent Cone pressure coefficient.
    if ( fabs( computeEmpiricalTangentConePressureCoefficient(
                   PI / 2.0, machNumber_) - 2.08961 ) > 1.0e-5 )
    {
        cerr << "Error in empirical Tangent Cone pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test high Mach base pressure coefficient.
    if ( fabs( computeHighMachBasePressure( machNumber_ )
               + 1.0 / pow( machNumber_, 2.0 ) ) > 1.0e-15 )
    {
        cerr << "Error in high Mach base pressure." << endl;
    }

    // Test empirical Tangent Wedge pressure coefficient.
    if ( fabs( computeEmpiricalTangentWedgePressureCoefficient(
                   PI / 2.0, machNumber_ ) - 2.38867 ) > 1.0e-5 )
    {
        cerr << "Error in empirical Tangent Wedge pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test freestream Prandtl-Meyer function.
    double freestreamPrandtlMeyerFunction_ = computePrandtlMeyerFunction(
                machNumber_, ratioOfSpecificHeats_ );

    if  ( fabs( freestreamPrandtlMeyerFunction_ - 106.9 * PI / 180.0 ) > 1.0e-3 )
    {
        cerr << "Error in freestream Prandtl-Meyer function." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test vacuum pressure coefficient.
    double vacuumPressureCoefficient_ = computeVacuumPressureCoefficient(
                machNumber_, ratioOfSpecificHeats_ );
    if ( fabs( vacuumPressureCoefficient_
               + 2.0 / ( ratioOfSpecificHeats_ * pow( machNumber_, 2.0 ) ) ) > 1.0e-15 )
    {
        cerr << "Error in vacuum pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Declare local variables.
    double angle_;
    double newtonianPressureCoefficient_;
    double prandtlMeterPressureCoefficient_;

    // Iterate over 10 angles and test pressure coefficients.
    for ( int i = 0; i < 10; i++ )
    {
        // Set angle.
        angle_ = static_cast< double > ( i ) * PI / 10.0;

        // Compute and compare Newtonian pressure coefficient.
        newtonianPressureCoefficient_
                = computeNewtonianPressureCoefficient( angle_ );

        if ( fabs( newtonianPressureCoefficient_ - 2.0 * pow( sin( angle_ ), 2.0 ) ) > 1.0e-15 )
        {
            cerr << "Error in Newtonian pressure coefficient." << endl;
            isAerodynamicsNamespaceBad = true;
        }

        // Compute Prandtl-Meyer pressure coefficient and test if it is not
        // lower than vacuum pressure coefficient.
        prandtlMeterPressureCoefficient_ =
            computePrandtlMeyerFreestreamPressureCoefficient( -1.0 * angle_,
            machNumber_, ratioOfSpecificHeats_, freestreamPrandtlMeyerFunction_ );
        if ( prandtlMeterPressureCoefficient_ - vacuumPressureCoefficient_ < -1.0e-15 )
        {
            cerr << "Error, Prandtl-Meyer pressure coefficient "
                 << "lower than vacuum pressure coefficient." << endl;
            isAerodynamicsNamespaceBad = true;
        }
    }

    // Test shock pressure ratio.
    if ( fabs( computeShockPressureRatio( machNumber_, ratioOfSpecificHeats_ )
            - 167.8 ) > 0.1 )
    {
        cerr << "Error in shock wave pressure ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock density ratio.
    if ( fabs( computeShockDensityRatio( machNumber_, ratioOfSpecificHeats_ )
            - 5.799 ) > 0.001 )
    {
        cerr << "Error in shock wave density ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock temperature ratio.
    if ( fabs( computeShockTemperatureRatio( machNumber_, ratioOfSpecificHeats_ )
            - 28.94 ) > 0.01 )
    {
        cerr << " Error in shock wave temperature ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock total pressure ratio.
    if ( fabs( computeShockTotalPressureRatio(
                   machNumber_, ratioOfSpecificHeats_, 287.058 ) - 0.001287 ) > 1.0e-6 )
    {
        cerr << "Error in shock wave total pressure ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isAerodynamicsNamespaceBad )
    {
        cerr << "testAerodynamicsNamespace failed!" << std::endl;
    }

    return isAerodynamicsNamespaceBad;
}
