/*! \file trajectoryDesignMethods.h
 *    This header file contains a base class for all trajectory design methods
 *    classes in Tudat.
 *
 *    Path              : /Astrodynamics/Trajectory/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2010
 *    Last modified     : 11 November, 2010
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
 *      YYMMDD    author        comment
 *      101111    E. Iorfida    First creation of code.
 *
 */

#ifndef TRAJECTORYDESIGNMETHODS_H
#define TRAJECTORYDESIGNMETHODS_H

// Include statements.
#include "linearAlgebra.h"

//! Trajectory design methods base class.
/*!
 * Trajectory design methods class.
 */
class TrajectoryDesignMethods
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    TrajectoryDesignMethods( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~TrajectoryDesignMethods( );

    virtual void execute( ) = 0;

protected:

private:
};

#endif // TRAJECTORYDESIGNMETHODS_H

// End of file.
