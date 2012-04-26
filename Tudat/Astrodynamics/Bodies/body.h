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
 *      100920    J. Melman         First creation of code.
 *      100929    K. Kumar          Minor comment changes.
 *      110115    J. Melman         Added set and get shape model functions.
 *      110815    K. Kumar          Added setMass( ) and getMass( ) functions.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120322    D. Dirkx          Moved Ephemeris member to here from CelestialBody
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#ifndef TUDAT_BODY_H
#define TUDAT_BODY_H

#include <map>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/States/state.h"
#include "Tudat/Astrodynamics/Bodies/Ephemeris/ephemeris.h"

namespace tudat
{
namespace bodies
{

//! Body base class.
/*!
 * Body base class.
 */
class Body
{
public:

    //! Create a body with a given mass.
    /*!
     * Creates a body with a given mass, 0 by default.
     * \param mass Mass of the body, default = 0.0.
     */
    Body( const double mass = 0.0 ) : mass_( mass ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Body( ) { }

    //! Set mass of body.
    /*!
     * Sets the mass of the body.
     * \param mass Mass.
     */
    void setMass( const double mass ) { mass_ = mass; }

    //! Get mass of body.
    /*!
     * Returns the mass of the body.
     * \return Mass.
     */
    double getMass( ) { return mass_; }

    //! Get ephemeris.
    /*!
     * Returns the body ephemeris.
     * \return Pointer to ephemeris.
     */
    boost::shared_ptr< ephemerides::Ephemeris > getEphemeris( ) { return ephemeris_; }

    //! Set ephemeris.
    /*!
     * Sets the ephemeris.
     */
    void setEphemeris( boost::shared_ptr< ephemerides::Ephemeris > ephemeris )
    {
        ephemeris_ = ephemeris;
    }

protected:

    //! Mass.
    /*!
     * Mass of body.
     */
    double mass_;

    //! Pointer to ephemeris.
    /*!
     * Pointer to ephemeris.
     */
    boost::shared_ptr< ephemerides::Ephemeris > ephemeris_;

private:
};

} // namespace bodies
} // namespace tudat

#endif // TUDAT_BODY_H
