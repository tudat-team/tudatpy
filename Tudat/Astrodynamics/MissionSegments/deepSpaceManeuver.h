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
 *      110224    E. Iorfida        First creation of code.
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
