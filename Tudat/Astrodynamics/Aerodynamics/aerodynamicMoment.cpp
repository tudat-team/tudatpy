/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *
 *    References
 *
 */

#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"

namespace tudat
{

//! Compute aerodynamic moment.
void AerodynamicMoment::computeMoment( State* pointerToState, double time )
{

    // Calculate moment.
    moment_ = dynamicPressure_ * pointerToAerodynamicCoefficientInterface_->getReferenceArea( ) *
            pointerToAerodynamicCoefficientInterface_->getReferenceLength( ) *
            pointerToAerodynamicCoefficientInterface_->getCurrentMomentCoefficients( );
}

} // namespace tudat
