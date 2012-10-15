/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100920    J. Melman         Creation of code.
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
