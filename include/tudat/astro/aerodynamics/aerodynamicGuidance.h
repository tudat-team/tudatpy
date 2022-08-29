/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_ENTRYGUIDANCE_H
#define TUDAT_ENTRYGUIDANCE_H

#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/basics/deprecationWarnings.h"

namespace tudat
{

namespace aerodynamics
{

/*!
 * Functionality is no longer supported, interfaces are retained to print error messages with referral to website for people using old code.
 */
class AerodynamicGuidance
{
public:

    //! Constructor.
    AerodynamicGuidance( )
    {
        utilities::printDeprecationError(
                    "tudatpy.numerical_simulation.propagation.AerodynamicGuidance",
                    "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#aerodynamic-guidance" );

    }

    virtual ~AerodynamicGuidance( ){ }

    virtual void updateGuidance( const double currentTime ) = 0;


protected:

    double currentAngleOfAttack_;

    double currentAngleOfSideslip_;

    double currentBankAngle_;

};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_ENTRYGUIDANCE_H
