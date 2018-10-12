/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CONSTANTEPHEMERIS_H
#define TUDAT_CONSTANTEPHEMERIS_H

#include <functional>
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

    using Ephemeris::getCartesianState;

    //! Constructor of a constant Ephemeris object
    /*!
     *  Constructor.
     *  \param constantStateFunction Function returning the constant state.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    ConstantEphemeris( const std::function< Eigen::Vector6d( ) > constantStateFunction,
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
    ConstantEphemeris( const Eigen::Vector6d constantState,
                       const std::string& referenceFrameOrigin = "SSB",
                       const std::string& referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
        { constantStateFunction_ = [ = ]( ){ return constantState; }; }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given time.
     * \param seconsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated
              (not used in this derived class)
     * \return Constant state given by constantStateFunction_
     */
    Eigen::Vector6d getCartesianState(
            const double seconsSinceEpoch = 0.0 )
    {
        return constantStateFunction_( );
    }

    //! Modifies the constant state
    /*!
     * Changes the constant state value to a new value.
     * \param newState New value for constant state.
     */
    void updateConstantState( const Eigen::Vector6d& newState )
    {
        constantStateFunction_ = [ = ]( ){ return newState; };
    }

private:

    //! Time-independent state function.
    /*!
     *  Function that returns a constant cartesian state.
     */
    std::function< Eigen::Vector6d( ) > constantStateFunction_;

};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CONSTANTEPHEMERIS_H
