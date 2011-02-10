/*! \file unitTestAerodynamicsNamespace.cpp
 *  This file contains the definition of the aerodynamics namespace unit test.
 *
 *  Path              : Astrodynamics/EnvironmentModels/Aerodynamics
 *  Version           : 1
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : L. Abdulkadir
 *  Affiliation       : TU Delft
 *  E-mail address    : L.Abdulkadir@student.tudelft.nl
 *
 *  Date created      : 08 November, 2010
 *  Last modified     : 10 February, 2011
 *
 *  References
 *
 *  Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill
 *  2001.
 *  Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *  edition, AIAA Education Series, 2006
 *
 *  Notes
 *
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *      YYMMDD    author        comment
 *      110208    D. Dirkx      First version of file.
 *      110210    L. Adulkadir  Code check.
 */

#include "unitTestAerodynamicsNamespace.h"

bool unit_tests::testAerodynamicsNameSpace( )
{
    bool isAerodynamicsNamespaceBad = 0;

    // Set default test conditions.
    double machNumber_ = 12.0;
    double ratioOfSpecificHeats_ = 1.4;

    // Test local to static pressure ratio.
    double localToStaticPressureRatio_ = aerodynamics::
                                         calculateLocalToStaticPressureRatio(
                                         machNumber_, ratioOfSpecificHeats_ );
    if( fabs(localToStaticPressureRatio_ - 1 / ( 0.1445E6 ))>1E-8 )
    {
        std::cerr << "Error in local to static pressure ratio." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test stagnation pressure coefficient.
    double stagnationPressureCoefficient_ = aerodynamics::calculateStagnationPressure(
            machNumber_, ratioOfSpecificHeats_ );
    if( fabs( stagnationPressureCoefficient_ - 1.83402 ) > 1E-5 )
    {
        std::cerr << "Error in stagnation pressure coefficient." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test modified Newtonian pressure coefficient.
    if( fabs( aerodynamics::calculateModifiedNewtonianPressureCoefficient(
            M_PI/2.0, stagnationPressureCoefficient_ ) -
              stagnationPressureCoefficient_ ) > 1E-15 )
    {
        std::cerr << "Error in modified Newtonian pressure coefficient." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test empirical Tangent Cone pressure coefficient.
    if( fabs( aerodynamics::calculateEmpiricalTangentConePressureCoefficient(
            M_PI/2.0, machNumber_) - 2.08961 ) > 1E-5)
    {
        std::cerr << "Error in empirical Tangent Cone pressure coefficient." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test high Mach base pressure coefficient.
    if( fabs( aerodynamics::calculateHighMachBasePressure( machNumber_ ) +
              1.0/( machNumber_ * machNumber_ ) ) > 1.0E-15 )
    {
        std::cerr << "Error in high Mach base pressure." << std::endl;
    }

    // Test empirical Tangent Wedge pressure coefficient.
    if( fabs( aerodynamics::calculateEmpiricalTangentWedgePressureCoefficient(
            M_PI/2.0, machNumber_) - 2.38867 ) > 1.0E-5)
    {
        std::cerr << "Error in empirical Tangent Wedge pressure coefficient." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test freestream Prandtl-Meyer function.
    double freestreamPrandtlMeyerFunction = aerodynamics::
                   evaluatePrandtlMeyerFunction(  machNumber_, ratioOfSpecificHeats_ );

    if( fabs( freestreamPrandtlMeyerFunction - 106.9 * M_PI / 180.0 ) > 1E-3 )
    {
        std::cerr << "Error in freestream Prandtl-Meyer function." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test vacuum pressure coefficient.
    double vacuumPressureCoefficient_ = aerodynamics::
                calculateVacuumPressureCoefficient( machNumber_,
                                                    ratioOfSpecificHeats_ );
    if( fabs( vacuumPressureCoefficient_ + 2.0/(ratioOfSpecificHeats_ *
                                       machNumber_ * machNumber_  ) ) > 1.0E-15 )
    {
        std::cerr << "Error in vacuum pressure coefficient." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Declare local variables.
    double angle_, newtonianPressureCoefficient_, prandtlMeterPressureCoefficient_;

    // Iterate over 10 angles and test pressure coefficients.
    for( int i = 0; i < 10; i++ )
    {
        // Set angle.
        angle_ = static_cast< double >( i ) * M_PI / 10;

        // Calsulate and compare Newtonian pressure coefficient.
        newtonianPressureCoefficient_ = aerodynamics::
                             calculateNewtonianPressureCoefficient( angle_ );
        if( fabs( newtonianPressureCoefficient_ - 2 * pow( sin( angle_ ),
                                                           2.0 ) ) > 1.0E-15 )
        {
            std::cerr << "Error in Newtonian pressure coefficient." << std::endl;
            isAerodynamicsNamespaceBad = 1;
        }

        // Calculate Prandtl-Meyer pressure coefficient and test if it is not
        // lower than vacuum pressure coefficient.
        prandtlMeterPressureCoefficient_ = aerodynamics::
            calculatePrandtlMeyerFreestreamPressureCoefficient( -1.0 * angle_,
            machNumber_, ratioOfSpecificHeats_, freestreamPrandtlMeyerFunction );
        if( prandtlMeterPressureCoefficient_ - vacuumPressureCoefficient_ < -1.0E-15 )
        {
            std::cerr << "Error, Prandtl-Meyer pressure coefficient" <<
                    "lower than vacuum pressure coefficient." <<std::endl;
            isAerodynamicsNamespaceBad = 1;
        }
    }

    // Test shock pressure ratio.
    if( fabs( aerodynamics::calculateShockPressureRatio( machNumber_,
                             ratioOfSpecificHeats_ ) - 167.8 ) > 0.1 )
    {
        std::cerr << "Error in shock wave pressure ratio." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test shock density ratio.
    if( fabs( aerodynamics::calculateShockDensityRatio( machNumber_,
                             ratioOfSpecificHeats_ ) - 5.799 ) > 0.001 )
    {
        std::cerr << "Error in shock wave density ratio." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test shock temperature ratio.
    if( fabs( aerodynamics::calculateShockTemperatureRatio( machNumber_,
                             ratioOfSpecificHeats_ ) - 28.94 ) > 0.01 )
    {
        std::cerr << " Error in shock wave temperature ratio." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    // Test shock total pressure ratio.
    if( fabs( aerodynamics::calculateShockTotalPressureRatio( machNumber_,
                                        ratioOfSpecificHeats_, 287.058) -
                                        0.001287 ) > 1.0E-6 )
    {
        std::cerr << "Error in shock wave total pressure ratio." << std::endl;
        isAerodynamicsNamespaceBad = 1;
    }

    return isAerodynamicsNamespaceBad;
}

