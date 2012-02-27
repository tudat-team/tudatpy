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
 *      110822    D. Dirkx          Removed pointer to double member, minor changes.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_FORCE_H
#define TUDAT_AERODYNAMIC_FORCE_H

#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

#include <iostream>
#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{

//! Aerodynamic force model.
/*!
 * Calculates the aerodynamic force based on C_D C_S and C_L, a reference lenght or area, the
 * local velocity and the local density.
 */
class AerodynamicForce : public ForceModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicForce( AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface ): dynamicPressure_( -0.0 )
    { pointerToAerodynamicCoefficientInterface_ = pointerToAerodynamicCoefficientInterface; }

    //! Get aerodynamic coefficient interface.
    /*!
     * Returns the pointer to the AerodynamicCoefficientInterface object.
     * \return Aerodynamic coefficient interface used to retrieve aerodynamic coefficients.
     */
    AerodynamicCoefficientInterface* getAerodynamicCoefficientInterface( )
    { return pointerToAerodynamicCoefficientInterface_; }

    //! Set dynamic pressure.
    /*!
     * Sets the dynamic pressure.
     * \param dynamicPressure Dynamic pressure.
     */
    void setDynamicPressure( double dynamicPressure )
    { dynamicPressure_ = dynamicPressure; }

    //! Get dynamic pressure.
    /*!
     * Returns the dynamic pressure.
     * \return Dynamic pressure.
     */
    double getDynamicPressure( ) { return dynamicPressure_; }

    //! Compute aerodynamic force.
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    void computeForce( State* pointerToState = NULL, double time = 0.0 )
    {
        force_ = dynamicPressure_ * pointerToAerodynamicCoefficientInterface_->getReferenceArea( )
                * pointerToAerodynamicCoefficientInterface_->getCurrentForceCoefficients( );
    }

protected:

private:

    //! Pointer to aerodynamic coefficient interface.
    /*!
     * Pointer to an aerodynamic coefficient interface.
     */
    AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface_;

    //! The dynamic pressure.
    /*!
     * The dynamic pressure.
     */
    double dynamicPressure_;
};

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_FORCE_H
