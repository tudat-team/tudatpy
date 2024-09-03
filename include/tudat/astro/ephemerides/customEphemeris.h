/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CUSTOMEPHEMERIS_H
#define TUDAT_CUSTOMEPHEMERIS_H

#include <functional>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/ephemerides/ephemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Ephemeris class that gives a custom (i.e. arbitrarily defined as a function of time) state.
template< typename TimeType = double, typename StateScalarType = double >
class CustomEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianState;
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    //! Constructor.
    /*!
     *  Constructor of an Ephemeris object from a custom function
     *  \param stateFunction Function returning the state as a function of time
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    CustomEphemeris( const std::function< StateType( const TimeType ) > stateFunction,
                     const std::string& referenceFrameOrigin = "SSB",
                     const std::string& referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
        stateFunction_( stateFunction ) { }


    Eigen::Vector6d getCartesianState(
        const double secondsSinceEpoch )
    {
        return stateFunction_( static_cast< TimeType >( secondsSinceEpoch ) ).template cast< double >( );
    }

    Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
        const double secondsSinceEpoch )
    {
        return stateFunction_( static_cast< TimeType >( secondsSinceEpoch ) ).template cast< long double >( );
    }

    Eigen::Vector6d getCartesianStateFromExtendedTime(
        const Time& time )
    {
        return stateFunction_( static_cast< TimeType >( time ) ).template cast< double >( );
    }

    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
        const Time& time )
    {
        return stateFunction_( static_cast< TimeType >( time ) ).template cast< long double >( );
    }

private:

    //! Time-independent state function.
    /*!
     *  Function that returns a constant cartesian state.
     */
    std::function< StateType( const TimeType ) > stateFunction_;

};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CUSTOMEPHEMERIS_H
