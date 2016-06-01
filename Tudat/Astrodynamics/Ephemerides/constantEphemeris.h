#ifndef TUDAT_CONSTANTEPHEMERIS_H
#define TUDAT_CONSTANTEPHEMERIS_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Ephemeris class that gives a constant (i.e. time independent) state.
class ConstantEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianStateFromEphemeris;

    //! Constructor of a constant Ephemeris object
    /*!
     *  Constructor.
     *  \param constantStateFunction Function returning the constant state.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantEphemeris( const boost::function< basic_mathematics::Vector6d( ) > constantStateFunction,
                       const std::string& referenceFrameOrigin = "SSB",
                       const std::string& referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
                constantStateFunction_( constantStateFunction ) { }

    //! Constructor of a constant Ephemeris object
    /*!
     *  Constructor
     *  \param constantState Constant state value.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantEphemeris( const basic_mathematics::Vector6d constantState,
                       const std::string& referenceFrameOrigin = "SSB",
                       const std::string& referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
        { constantStateFunction_ = boost::lambda::constant( constantState ); }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given time.
     * \param seconsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated
              (not used in this derived class)
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return Constant state given by constantStateFunction_
     */
    basic_mathematics::Vector6d getCartesianStateFromEphemeris(
            const double seconsSinceEpoch = 0.0,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000)
    {
        return constantStateFunction_( );
    }

    //! Modifies the constant state
    /*!
     * Changes the constant state value to a new value.
     * \param newState New value for constant state.
     */
    void updateConstantState( const basic_mathematics::Vector6d& newState )
    {
        constantStateFunction_ = boost::lambda::constant( newState );
    }

private:

    //! Time-independent state function.
    /*!
     *  Function that returns a constant cartesian state.
     */
    boost::function< basic_mathematics::Vector6d( ) > constantStateFunction_;

};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CONSTANTEPHEMERIS_H
