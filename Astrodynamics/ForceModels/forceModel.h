/*! \file forceModel.h
 *    Header file that defines the base class for all force models included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 6
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
 *    Last modified     : 2 February, 2011
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
 *                                    and comments.
 *      110113    K. Kumar            Modified CelestialBody object to pointer;
 *                                    minor comment changes.
 *      110119    K. Kumar            Changed computeStateDerivatives() to
 *                                    computeForce().
 *      110202    K. Kumar            Updated code to make use of State class.
 */

#ifndef FORCEMODEL_H
#define FORCEMODEL_H

// Include statements.
#include "celestialBody.h"
#include "linearAlgebra.h"
#include "cartesianPositionElements.h"

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

    //! Compute force per unit mass.
    /*!
     * Computes the forces per unit mass.
     * \param pointerToState Pointer to a State object.
     */
    virtual VectorXd& computeForce( State* pointerToState ) = 0;

protected:

    //! Pointer to object of Celestial Body class.
    /*!
     * Pointer to object of Celestial Body class for gravity field expansion.
     */
    CelestialBody* celestialBody_;

    //! Force per unit mass.
    /*!
     * Force per unit mass.
     */
    VectorXd forcePerUnitMass_;

private:
};

#endif // FORCEMODEL_H

// End of file.
