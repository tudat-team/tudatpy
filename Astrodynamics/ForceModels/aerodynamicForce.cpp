/*! \file aerodynamicForce.cpp
*    This header file contains the aerodynamic forces model included in Tudat.
*
*    Path              : /Astrodynamics/ForceModels/
*    Version           : 2
*    Check status      : Checked
*
*    Checker           : F. M. Engelen
*    Affiliation       : Delft University of Technology
*    E-mail address    : F.M.Engelen@student.tudelft.nl
*
*    Checker           : D. Dirkx
*    Affiliation       : Delft University of Technology
*    E-mail address    : D..Dirkx@tudelft.nl
*
*    Date created      : 17 May, 2011
*    Last modified     : 22 August, 2011
*
*    References
*
*    Notes
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
*      110617    F.M. Engelen      First creation of code.
*      110822    D.Dirkx           Removed pointer to double member, minor changes.
*/

// Include statements.
#include "aerodynamicForce.h"

//! Default constructor.
AerodynamicForce::AerodynamicForce( ): pointerToAerodynamicCoefficientInterface_( NULL ),
                                       dynamicPressure_( -0.0 )
{
}

//! Default destructor.
AerodynamicForce::~AerodynamicForce( )
{
}

//! Set aerodynamic coefficient interface.
void AerodynamicForce::setAerodynamicCoefficientInterface(
        AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface )
{
    pointerToAerodynamicCoefficientInterface_ = pointerToAerodynamicCoefficientInterface;
}

//! Get aerodynamic coefficient interface.
AerodynamicCoefficientInterface*
        AerodynamicForce::getAerodynamicCoefficientInterface( )
{
    return pointerToAerodynamicCoefficientInterface_;
}

//! Set dynamic pressure.
void AerodynamicForce::setDynamicPressure( const double& dynamicPressure )
{
    dynamicPressure_ = dynamicPressure;
}

//! Get dynamic pressure.
double& AerodynamicForce::getDynamicPressure( )
{
    return dynamicPressure_;
}

//! Compute aerodynamic force.
void AerodynamicForce::computeForce( State* pointerToState )
{
    // Calculate the aerodynamic forces F_x = 1/2 * C_x * rho * V^2 * S etc...
    force_ = dynamicPressure_ *
                        pointerToAerodynamicCoefficientInterface_->getReferenceArea( ) *
                        pointerToAerodynamicCoefficientInterface_->getCurrentForceCoefficients( );
}

// End of file.
