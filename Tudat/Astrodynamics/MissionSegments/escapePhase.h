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
 *      110201    E. Iorfida        First creation of code.
 *
 *
 *    References
 *
 */

#ifndef TUDAT_ESCAPE_PHASE_H
#define TUDAT_ESCAPE_PHASE_H

#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

namespace tudat
{

//! Escape phase class.
/*!
 * Escape phase class.
 */
class EscapePhase : public EscapeAndCapture
{
public:

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param escapePhase Escape phase.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, EscapePhase& escapePhase )
    {
        stream << "The computed delta-V is: " << escapePhase.computeDeltaV( ) << std::endl;
        return stream;
    }

protected:

private:
};

} // namespace tudat

#endif // TUDAT_ESCAPE_PHASE_H
