/*! \file forceModel.h
 *    Header file that defines the base class for all force models included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/ForceModel/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 14 September, 2010
 *    Last modified     : 29 September, 2010
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
 *      YYMMDD    author              comment
 *      100914    K. Kumar            File created.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor corrections to include statements
 *                                     and comments.
 */

#ifndef FORCEMODEL_H
#define FORCEMODEL_H

// Include statements.
#include "celestialBody.h"
#include "linearAlgebra.h"

//! Force model class.
/*!
 * Base class for all force models.
 */
class ForceModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ForceModel( ){};

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~ForceModel( ){};

    //! Compute state derivatives.
    /*!
     * This function computes the state derivatives.
     * \param stateVector State vector.
     * \param stateDerivativeVector Vector containing derivative of
     *          stateVector.
     */
    virtual void computeStateDerivatives( VectorXd& stateVector,
                                          VectorXd& stateDerivativeVector ) =0;

    //! Object of Celestial Body class.
    /*!
     * Object of Celestial Body class for gravity field expansion.
     */
    CelestialBody celestialBody_;

protected:

private:
};

#endif // FORCEMODEL_H

// End of file.
