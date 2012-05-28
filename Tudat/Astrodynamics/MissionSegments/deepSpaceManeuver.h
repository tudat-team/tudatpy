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
 *      110224    E. Iorfida        Creation of code.
 *      110406    K. Kumar          Minor modifications.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120417    T. Secretin       Moved set functions to constructor.
 *
 *    References
 *
 */

#ifndef TUDAT_DEEP_SPACE_MANEUVER_H
#define TUDAT_DEEP_SPACE_MANEUVER_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/States/state.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

//! Deep space maneuver base class.
/*!
 * Deep space maneuver class.
 */
class DeepSpaceManeuver
{
public:

    //! Typedef of shared pointer to state.
    /*!
     * Typedef of shared pointer to state.
     */
    typedef boost::shared_ptr< states::State > StatePointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    DeepSpaceManeuver( const double deltaV,
                       const double timeOfDeepSpaceManeuver,
                       StatePointer state )
        : deltaV_( deltaV ),
          timeOfDeepSpaceManeuver_( timeOfDeepSpaceManeuver ),
          state_( state )
    { }

    //! Get time of deep space maneuver event.
    /*!
     * Returns the the time of deep space maneuver event.
     * \return Time of deep space maneuver event.
     */
    double getTime( ) { return timeOfDeepSpaceManeuver_; }

    //! Get state at deep space maneuver event.
    /*!
     * Returns a pointer to state at deep space maneuver event.
     * \return Pointer to state at deep space maneuver event.
     */
    StatePointer getState( ) { return state_; }

    //! Get delta-V of deep space maneuver event.
    /*!
     * Returns delta-V of deep space maneuver event.
     * \return Delta-V of deep space maneuver event.
     */
    double getDeltaV( ) { return deltaV_; }

protected:

private:

    //! Delta-V of deep space maneuver event.
    /*!
     * Delta-V of deep space maneuver event.
     */
    double deltaV_;

    //! Time of deep space maneuver event.
    /*!
     * Time of deep space maneuver event.
     */
    double timeOfDeepSpaceManeuver_;

    //! Shared pointer to state at deep space maneuver event.
    /*!
     * Shared pointer to state at deep space maneuver event.
     */
    StatePointer state_;
};

} // astrodynamics
} // mission_segments
} // namespace tudat

#endif // TUDAT_DEEP_SPACE_MANEUVER_H
