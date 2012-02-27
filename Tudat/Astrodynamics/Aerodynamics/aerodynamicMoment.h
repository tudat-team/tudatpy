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
 *      110822    D. Dirkx          Removed pointer to double member, removed location of center
 *                                  of mass, minor changes.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_MOMENT_H
#define TUDAT_AERODYNAMIC_MOMENT_H

#include <iostream>
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/pureMomentModel.h"

namespace tudat
{

class AerodynamicMoment : public PureMomentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicMoment( AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface ) :
        pointerToAerodynamicCoefficientInterface_( pointerToAerodynamicCoefficientInterface ),
        dynamicPressure_( -0.0 ) { }

    //! Set aerodynamic coefficient interface.
    /*!
     * Returns the pointer to the AerodunamicCoefficientDatabase object.
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

    //! Compute aerodynamic moment
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    void computeMoment( State* pointerToState = NULL, double time = 0.0 );

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

#endif // TUDAT_AERODYNAMIC_MOMENT_H
