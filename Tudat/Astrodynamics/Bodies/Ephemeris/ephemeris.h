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
 *      110128    K. Kumar          First creation of code.
 *      110221    K. Kumar          Updated code to work as base class for derived ephemeris
 *                                  classes.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *
 *    References
 *
 */

#ifndef TUDAT_EPHEMERIS_H
#define TUDAT_EPHEMERIS_H

#include <sstream>

#include <TudatCore/Mathematics/BasicMathematics/linearAlgebra.h>

namespace tudat
{
namespace ephemerides
{

//! Ephemeris base class.
/*!
 * Ephemeris base class.
 */
class Ephemeris
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Ephemeris( ) { }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param julianDate Julian date given in Julian days.
     * \return State from ephemeris.
     */
    virtual Eigen::VectorXd getCartesianStateFromEphemeris( const double julianDate ) = 0;

protected:

private:
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_EPHEMERIS_H
