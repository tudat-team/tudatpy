/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_TABULATEDEPHEMERIS_H
#define TUDAT_TABULATEDEPHEMERIS_H

#include <map>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace ephemerides
{

//! Class that determines an ephemeris from tabulated data.
/*!
 *  Class that determines an ephemeris from tabulated data, by using numerical interpolation of
 *  this data. Required input to this class is a OneDimensionalInterpolator, which may be reset.
 *  This class may for instance be used for setting the numerically integrated state of a body
 *  as its 'new' ephemeris
 */
template< typename StateScalarType = double, typename TimeType = double >
class TabulatedCartesianEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianState;
    using Ephemeris::getCartesianLongState;
    using Ephemeris::getCartesianStateFromExtendedTime;
    using Ephemeris::getCartesianLongStateFromExtendedTime;

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VariableStateType;

    //! Typedef for state interpolator
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator
    < TimeType, StateType  > > StateInterpolatorPointer;
    typedef std::shared_ptr< interpolators::OneDimensionalInterpolator
    < TimeType, VariableStateType  > > VariableStateInterpolatorPointer;

    //! Constructor, sets data interpolator and frame data.
    /*!
     *  Constructor, sets data interpolator and frame data.
     *  \param interpolator Interpolator that returns the interpolated state as a function of time.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    TabulatedCartesianEphemeris(
            const StateInterpolatorPointer interpolator,
            const std::string referenceFrameOrigin = "SSB",
            const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ), interpolator_( interpolator )
    {  }

    TabulatedCartesianEphemeris(
            const VariableStateInterpolatorPointer interpolator,
            const std::string referenceFrameOrigin = "SSB",
            const std::string referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation )
    {
       interpolator_ = interpolators::convertBetweenStaticDynamicEigenTypeInterpolators<
               TimeType, StateScalarType, Eigen::Dynamic, 1, 6, 1 >(
                   interpolator );
    }


    //! Destructor
    /*!
     *  Destructor
     */
    ~TabulatedCartesianEphemeris( ){ }

    //! Function to reset the state interpolator.
    /*!
     *  Function to reset the state interpolator, for instance following an update of the states of
     *  the body after a new numerical integration.
     *  \param interpolator New interpolator that returns the interpolated state as a function of
     *  time.
     */
    void resetInterpolator( const StateInterpolatorPointer interpolator )
    {
        interpolator_ = interpolator;
    }

    void resetInterpolator( const VariableStateInterpolatorPointer interpolator )
    {
        interpolator_ = interpolators::convertBetweenStaticDynamicEigenTypeInterpolators<
                TimeType, StateScalarType, Eigen::Dynamic, 1, 6, 1 >(
                    interpolator );
    }

    //! Get cartesian state from ephemeris.
    /*!
     * Returns cartesian state from ephemeris, as calculated from interpolator_.
     * \param secondsSinceEpoch Seconds since epoch.
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch );

    //! Get cartesian state from ephemeris (in long double precision).
    /*!
     * Returns cartesian state from ephemeris  (in long double precision), as calculated from interpolator_. For
     * double StateScalarType class template argument, this function returns the double precision interpolated values,
     * cast to long double. Only for long double StateScalarType argument is this function used to its fullest.
     * \param secondsSinceEpoch Seconds since epoch.
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
            const double secondsSinceEpoch );

    //! Get cartesian state from ephemeris (in double precision from Time input).
    /*!
     * Returns cartesian state from ephemeris  (in double precision from Time input), as calculated from interpolator_.
     * \param time Time at which ephemeris is to be evaluated
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Vector6d getCartesianStateFromExtendedTime(
            const Time& time );

    //! Get cartesian state from ephemeris (in long double precision from Time input).
    /*!
     * Returns cartesian state from ephemeris  (in long double precision from Time input), as calculated from interpolator_.
     * For double StateScalarType class template argument, this function returns the double precision interpolated values,
     * cast to long double. Only for long double StateScalarType argument is this function used to its fullest.
     * \param time Time at which ephemeris is to be evaluated
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
            const Time& time );


    //! Function to return the interpolator
    /*!
     *  Function to return the interpolator that is to be used to calculate the state.
     *  \return Interpolator that is to be used to calculate the state.
     */
    StateInterpolatorPointer getInterpolator( )
    {
        return interpolator_;
    }

    VariableStateInterpolatorPointer getDynamicVectorSizeInterpolator( )
    {
        return interpolators::convertBetweenStaticDynamicEigenTypeInterpolators<
                TimeType, StateScalarType, 6, 1, Eigen::Dynamic, 1 >(
                    interpolator_ );
    }

    //! Function that retrieves the time interval at which this ephemeris can be safely interrogated
    /*!
     * Function that retrieves the time interval at which this ephemeris can be safely interrogated. The interval
     * on which the interpolator inside this object is valid is checked and returned
     * \return The time interval at which the tabulated ephemeris can be safely interrogated
     */
    std::pair< double, double > getSafeInterpolationInterval( )
    {
        std::pair< double, double > safeInterpolationInterval;

        // Check interpolator type. If interpolator is not a Lagrange interpolator, return full domain
        if( std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, double > >(
                    interpolator_ ) == nullptr &&
                std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, long double > >(
                    interpolator_ ) == nullptr )
        {
            safeInterpolationInterval.first = interpolator_->getIndependentValues( ).at( 0 );
            safeInterpolationInterval.second = interpolator_->getIndependentValues( ).at(
                        interpolator_->getIndependentValues( ).size( ) - 1 );
        }
        // If interpolator is a Lagrange interpolator, return full domain minus edges where interpolator has reduced accuracy
        else if( std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, double > >(
                     interpolator_ ) != nullptr )
        {
            int numberOfNodes =
                    std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, double > >(
                        interpolator_ )->getNumberOfStages( );

            safeInterpolationInterval.first = interpolator_->getIndependentValues( ).at( 0 + numberOfNodes / 2 + 1 );
            safeInterpolationInterval.second = interpolator_->getIndependentValues( ).at(
                        interpolator_->getIndependentValues( ).size( ) - 1 - ( + numberOfNodes / 2 + 1 ) );
        }
        else if( std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, long double > >(
                     interpolator_ ) != nullptr )
        {
            int numberOfNodes =
                    std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< TimeType, StateType, long double > >(
                        interpolator_ )->getNumberOfStages( );
            safeInterpolationInterval.first = interpolator_->getIndependentValues( ).at( 0 + numberOfNodes / 2 + 1 );
            safeInterpolationInterval.second = interpolator_->getIndependentValues( ).at(
                        interpolator_->getIndependentValues( ).size( ) - 1 - ( + numberOfNodes / 2 + 1 ) );
        }
        return safeInterpolationInterval;
    }

private:

    //! Interpolator that returns body state as a function of time.
    /*!
     *  Interpolator that returns body state as a function of time by calling the interpolate
     *  function (i.e. time as independent variable and states as dependent variables ).
     */
    StateInterpolatorPointer interpolator_;
};


//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
/*!
 *  Function to check whether an ephemeris is a (type of) tabulated ephemeris, it checks all typical combinations of
 *  class template arguments are returns true if a dynamic cast is succesful
 *  \param ephemeris Ephemeris pointer for which it is to be checked whether it is a tabulated ephemeris
 *  \return True if ephemeris is a tabulated ephemeris
 */
bool isTabulatedEphemeris( const std::shared_ptr< Ephemeris > ephemeris );

//! Function that retrieves the time interval at which a tabulated ephemeris can be safely interrogated
/*!
 * Function that retrieves the time interval at which a tabulated ephemeris can be safely interrogated. The interval
 * on which the interpolator inside this object is valid is checked and returned
 * \param ephemeris Ephemeris model for which the interval is to be determined. AN exception is thrown if this is not
 * a tabulated ephemeris
 * \return The time interval at which the tabulated ephemeris can be safely interrogated
 */
std::pair< double, double > getTabulatedEphemerisSafeInterval( const std::shared_ptr< Ephemeris > ephemeris );

//! Function to create an empty (dummy) tabulated ephemeris
/*!
 *  Function to create an empty (dummy) tabulated ephemeris. This is used when for instance propagating a body for which
 *  the propagated result shopuld be saved in the epehemris, but no a priori ephemeris is available
 *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
 *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
 *  \return Empty tabulated ephemeris with given reference frame settings
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< Ephemeris > createEmptyTabulatedEphemeris(
        const std::string referenceFrameOrigin = "SSB",
        const std::string referenceFrameOrientation = "ECLIPJ2000"  )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    return std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateType > >( ),
                referenceFrameOrigin, referenceFrameOrientation );
}

//! Create a tabulated ephemeris from a given ephemeris model and interpolation settings
/*!
 * Create a tabulated ephemeris from a given ephemeris model and interpolation settings
 * \param ephemerisToInterrogate Ephemeris model from which the tabulated model is tpo be synthesized
 * \param startTime Start time for tabulated ephemeris
 * \param endTime End time for tabulated ephemeris
 * \param timeStep Constant time step for tabulated ephemeris
 * \param interpolatorSettings Interpolation settings for tabulated ephemeris
 * \return Tabulated ephemeris, as synthesized from a given ephemeris model and interpolation settings
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< Ephemeris > getTabulatedEphemeris(
        const std::shared_ptr< Ephemeris > ephemerisToInterrogate,
        const TimeType startTime,
        const TimeType endTime,
        const TimeType timeStep,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) )
{

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    // Create state map that is to be interpolated
    std::map< TimeType, StateType >  stateMap;
    TimeType currentTime = startTime;
    while( currentTime <= endTime )
    {
        stateMap[ currentTime ] = ephemerisToInterrogate->getTemplatedStateFromEphemeris<
                StateScalarType, TimeType >( currentTime );
        currentTime += timeStep;
    }

    // Create tabulated ephemeris model
    auto ephemeris = std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                interpolators::createOneDimensionalInterpolator( stateMap, interpolatorSettings ),
                ephemerisToInterrogate->getReferenceFrameOrigin( ),
                ephemerisToInterrogate->getReferenceFrameOrientation( ) );

    return ephemeris;

}

} // namespace ephemerides

} // namespace tudat
#endif // TUDAT_TABULATEDEPHEMERIS_H
