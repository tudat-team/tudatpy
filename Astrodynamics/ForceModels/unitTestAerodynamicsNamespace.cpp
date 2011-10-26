/*! \file unitTestAerodynamicsNamespace.cpp
 *    This file contains the definition of the aerodynamics namespace unit test.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : L. Abdulkadir
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *    Date created      : 8 November, 2010
 *    Last modified     : 11 February, 2011
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition,
 *          McGraw Hill, 2001.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *          edition, AIAA Education Series, 2006.
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
 *      110208    D. Dirkx          First version of file.
 *      110210    L. Adulkadir      Code check.
 *      110211    K. Kumar          Corrected Doxygen errors; corrected layout errors; corrected
 *                                  double precision; updated function-naming.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/ForceModels/aerodynamics.h"

//! Test aerodynamics namespace.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using std::fabs;
    using std::pow;
    using std::sin;

    bool isAerodynamicsNamespaceBad = false;

    // Set default test conditions.
    double machNumber_ = 12.0;
    double ratioOfSpecificHeats_ = 1.4;

    // Test local to static pressure ratio.
    double localToStaticPressureRatio_ = aerodynamics::computeLocalToStaticPressureRatio(
                machNumber_, ratioOfSpecificHeats_ );

    if ( fabs( localToStaticPressureRatio_ - 1.0 / ( 0.1445e6 ) ) > 1.0e-8 )
    {
        cerr << "Error in local to static pressure ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test stagnation pressure coefficient.
    double stagnationPressureCoefficient_
            = aerodynamics::computeStagnationPressure( machNumber_, ratioOfSpecificHeats_ );

    if ( fabs( stagnationPressureCoefficient_ - 1.83402 ) > 1.0e-5 )
    {
        cerr << "Error in stagnation pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test modified Newtonian pressure coefficient.
    if ( fabs( aerodynamics::computeModifiedNewtonianPressureCoefficient(
                   M_PI / 2.0, stagnationPressureCoefficient_ ) -
               stagnationPressureCoefficient_ ) > 1.0e-15 )
    {
        cerr << "Error in modified Newtonian pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test empirical Tangent Cone pressure coefficient.
    if ( fabs( aerodynamics::computeEmpiricalTangentConePressureCoefficient(
                   M_PI / 2.0, machNumber_) - 2.08961 ) > 1.0e-5 )
    {
        cerr << "Error in empirical Tangent Cone pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test high Mach base pressure coefficient.
    if ( fabs( aerodynamics::computeHighMachBasePressure( machNumber_ )
               + 1.0 / pow( machNumber_, 2.0 ) ) > 1.0e-15 )
    {
        cerr << "Error in high Mach base pressure." << endl;
    }

    // Test empirical Tangent Wedge pressure coefficient.
    if ( fabs( aerodynamics::computeEmpiricalTangentWedgePressureCoefficient(
                   M_PI / 2.0, machNumber_ ) - 2.38867 ) > 1.0e-5 )
    {
        cerr << "Error in empirical Tangent Wedge pressure coefficient." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test freestream Prandtl-Meyer function.
    double freestreamPrandtlMeyerFunction_ = aerodynamics::computePrandtlMeyerFunction(
                machNumber_, ratioOfSpecificHeats_ );

    if  ( fabs( freestreamPrandtlMeyerFunction_ - 106.9 * M_PI / 180.0 ) > 1.0e-3 )
    {
        cerr << "Error in freestream Prandtl-Meyer function." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test vacuum pressure coefficient.
    double vacuumPressureCoefficient_ = aerodynamics::computeVacuumPressureCoefficient(
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
        angle_ = static_cast< double > ( i ) * M_PI / 10.0;

        // Compute and compare Newtonian pressure coefficient.
        newtonianPressureCoefficient_
                = aerodynamics::computeNewtonianPressureCoefficient( angle_ );

        if ( fabs( newtonianPressureCoefficient_ - 2.0 * pow( sin( angle_ ), 2.0 ) ) > 1.0e-15 )
        {
            cerr << "Error in Newtonian pressure coefficient." << endl;
            isAerodynamicsNamespaceBad = true;
        }

        // Compute Prandtl-Meyer pressure coefficient and test if it is not
        // lower than vacuum pressure coefficient.
        prandtlMeterPressureCoefficient_ = aerodynamics::
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
    if ( fabs( aerodynamics::computeShockPressureRatio( machNumber_, ratioOfSpecificHeats_ )
            - 167.8 ) > 0.1 )
    {
        cerr << "Error in shock wave pressure ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock density ratio.
    if ( fabs( aerodynamics::computeShockDensityRatio( machNumber_, ratioOfSpecificHeats_ )
            - 5.799 ) > 0.001 )
    {
        cerr << "Error in shock wave density ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock temperature ratio.
    if ( fabs( aerodynamics::computeShockTemperatureRatio( machNumber_, ratioOfSpecificHeats_ )
            - 28.94 ) > 0.01 )
    {
        cerr << " Error in shock wave temperature ratio." << endl;
        isAerodynamicsNamespaceBad = true;
    }

    // Test shock total pressure ratio.
    if ( fabs( aerodynamics::computeShockTotalPressureRatio(
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

// End of file.
