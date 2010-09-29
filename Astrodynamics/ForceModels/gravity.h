/*! \file gravity.h
 *    Header file that defines the gravity force model included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 16 September, 2010
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
 *      100916    K. Kumar            File created.
 *      100916    K. Kumar            Filename modified.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor corrections to include statements
                                      and comments.
 */

#ifndef GRAVITY_H
#define GRAVITY_H

// Include statements.
#include <cmath>
#include "forceModel.h"
#include "celestialBody.h"

//! Gravity force class.
/*!
 * Class containing the gravity force model.
 */
class Gravity : public ForceModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    Gravity();

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~Gravity();

    //! Set body for gravity field expansion.
    /*!
     * This function sets the body for gravity field expansion.
     * \param celestialBody Celestial body which is set.
     */
    void setBody( CelestialBody& celestialBody );

    //! Compute state derivatives for gravity field expansion.
    /*!
     * This function computes the state derivatives for gravity field
     * expansion.
     * \param stateVector State vector of size 6; first three cartesian
     * position coordinates, followed by three cartesian velocity coordinates.
     * \param stateDerivativeVector Derivative of stateVector from gravity
     * field; first three entries are equal to velocity from state, second
     * three are computed from gravity field and position.
     */
    void computeStateDerivatives( Vector& stateVector,
                                  Vector& stateDerivativeVector );

protected:

private:
};

#endif // GRAVITY_H

// End of file.
