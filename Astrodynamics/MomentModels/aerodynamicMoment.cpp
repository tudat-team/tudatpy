/*! \file aerodynamicMoment.cpp
*    This source file contains the aerodynamic moment model included in Tudat.
*
*    Path              : /Astrodynamics/MomentModels/
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
*      110617    F.M. Engelen      First creation of code.
*      110822    D.Dirkx           Removed pointer to double member, removed location of center
*                                  of mass, minor changes.
*/

// Include statements.
#include <cmath>
#include "Astrodynamics/MomentModels/aerodynamicMoment.h"

//! Compute aerodynamic moment.
void AerodynamicMoment::computeMoment( State* pointerToState )
{
    // Calculate moment.
    moment_ = dynamicPressure_ *
              pointerToAerodynamicCoefficientInterface_->getReferenceArea( ) *
              pointerToAerodynamicCoefficientInterface_->getReferenceLength( ) *
              pointerToAerodynamicCoefficientInterface_->getCurrentMomentCoefficients( );

    // Check if pointer to forcemodel is set.
    if ( pointerToForceModel_ != NULL )
    {
        // Update moment with additional moment due to force.
        Vector3d force_ = pointerToForceModel_->getForce( );
        moment_ += forceApplicationArm_.cross( force_ );
    }
}

// End of file.
