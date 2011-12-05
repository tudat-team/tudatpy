/*! \file ephemeris.h
 *    This header file contains the definition of a base class for all ephemeris classes in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 28 January, 2011
 *    Last modified     : 21 February, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 */

#ifndef EPHEMERIS_H
#define EPHEMERIS_H

// Include statements.
#include <sstream>
#include "Astrodynamics/States/cartesianElements.h"
#include "Input/textFileReader.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
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

    //! Bodies with ephemeris data.
    /*!
     * Bodies with ephemeris data.
     */
    BodiesWithEphemerisData bodyWithEphemerisData_;

    //! Ephemeris text file reader.
    /*!
     * Ephemeris text file reader.
     */
    TextFileReader ephemerisTextFileReader_;

    //! String stream for ephemeris line data.
    /*!
     * String stream for ephemeris line data.
     */
    std::stringstream ephemerisLineData_;

private:
};

}

#endif // EPHEMERIS_H

// End of file.
