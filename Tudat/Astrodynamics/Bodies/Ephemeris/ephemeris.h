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
 *
 *    References
 *
 */

#ifndef TUDAT_EPHEMERIS_H
#define TUDAT_EPHEMERIS_H

#include <sstream>
#include "Tudat/Astrodynamics/States/cartesianElements.h"

namespace tudat
{

//! Ephemeris base class.
/*!
 * Ephemeris base class.
 */
class Ephemeris
{
public:

    //! Bodies with ephemeris data.
    /*!
     * Bodies with ephemeris data.
     */
    enum BodiesWithEphemerisData
    { mercury, venus, earthMoonBarycenter, mars, jupiter, saturn, uranus, neptune, pluto };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    Ephemeris( ) : ephemerisLineData_( ) { }

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
    virtual CartesianElements* getStateFromEphemeris( double julianDate ) = 0;

protected:

    //! String stream for ephemeris line data.
    /*!
     * String stream for ephemeris line data.
     */
    std::stringstream ephemerisLineData_;

private:
};

} // namespace tudat

#endif // TUDAT_EPHEMERIS_H
