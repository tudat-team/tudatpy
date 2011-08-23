/*! \file aerodynamicMoment.h
*    This header file contains the aerodynamic moment model included in Tudat.
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
*      110822    D. Dirkx          Removed pointer to double member, removed location of center
*                                  of mass, minor changes.
*/

#ifndef AERODYNAMICMOMENT_H
#define AERODYNAMICMOMENT_H

// Include statements.
#include <cmath>
#include "aerodynamicCoefficientInterface.h"
#include "momentModel.h"

class AerodynamicMoment : public MomentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicMoment( );

    //! Default destructor.
    /*!
     * Default constructor.
     */
    ~AerodynamicMoment( );

    //! Set aerodynamic coefficient interface.
    /*!
     * Sets the aerodynamic coefficient interface.
     * \param pointerToAerodynamicCoefficientInterface Aerodynamic coefficient interface
     *          used to retrieve aerodynamic coefficients.
     */
    void setAerodynamicCoefficientInterface(
            AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface );

    //! Set aerodynamic coefficient interface.
    /*!
     * Returns the pointer to the AerodunamicCoefficientDatabase object.
     * \return Aerodynamic coefficient interface used to retrieve aerodynamic coefficients.
     */
    AerodynamicCoefficientInterface* getAerodynamicCoefficientInterface( );

    //! Set dynamic pressure.
    /*!
     * Sets the dynamic pressure.
     * \param dynamicPressure Dynamic pressure.
     */
    void setDynamicPressure( const double& dynamicPressure );

    //! Get dynamic pressure.
    /*!
     * Returns the dynamic pressure.
     * \return Dynamic pressure.
     */
    double& getDynamicPressure( );

    //! Compute aerodynamic moment
    /*!
     * Computes the aerodynamic moment.
     * \param pointerToState Pointer to current state.
     */
    void computeMoment( State* pointerToState );

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

#endif // AERODYNAMICMOMENT_H

// End of file.
