#ifndef CONSTANTEPHEMERIS_H
#define CONSTANTEPHEMERIS_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Ephemeris class that gives a constant (i.e. time independent state).
/*!
 *  Ephemeris class that gives a constant (i.e. time independent state).
 */
class ConstantEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianStateFromEphemeris;

    //! Constructor
    /*!
     *  \param constantStateFunction Function returning the constant state.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantEphemeris( const boost::function< basic_mathematics::Vector6d( ) > constantStateFunction,
                       const std::string referenceFrameOrigin = "SSB",
                       const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ), constantStateFunction_( constantStateFunction ) { }

    //! Constructor
    /*!
     *  \param constantState Constant state value.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantEphemeris( const basic_mathematics::Vector6d constantState,
                       const std::string referenceFrameOrigin = "SSB",
                       const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
        { constantStateFunction_ = boost::lambda::constant( constantState ); }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given time.
     * \param ephemerisTime Seconds since epoch at which ephemeris is to be evaluated (not used in this derived class)
     * \return Constant state given by constantStateFunction_
     */
    basic_mathematics::Vector6d getCartesianStateFromEphemeris(
            const double ephemerisTime = 0.0, const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000)
    {
        return constantStateFunction_( );
    }

    //! Modifies the constant state
    /*!
     * Changes the constant state value to a new value.
     * \param newState New value for constant state.
     */
    void updateConstantState( const basic_mathematics::Vector6d newState )
    {
        constantStateFunction_ = boost::lambda::constant( newState );
    }
private:

    //! Time-independent state function.
    /*!
     *  Function that returns a (constant) cartesian state.
     */
    boost::function< basic_mathematics::Vector6d( ) > constantStateFunction_;

};

}

}

#endif // CONSTANTEPHEMERIS_H
