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

#ifndef TUDAT_TIME_EPHEMERIS_H
#define TUDAT_TIME_EPHEMERIS_H

#include <functional>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/timeType.h"

namespace tudat
{

//! Combine chained time-difference functions with configurable evaluation epochs.
/*!
 *  This utility evaluates a sequence of conversion functions and accumulates the resulting
 *  time differences. Each conversion function can be evaluated at either the original epoch
 *  or at one of the intermediate corrected epochs.
 *  \param timeDifferenceFunctions Ordered list of individual conversion functions.
 *  \param evaluationStepIndices Index list defining at which intermediate epoch each
 *  conversion function is evaluated.
 *  \param time Base input epoch.
 *  \return Total time difference obtained by summing all conversion steps.
 */
template< typename TimeType >
TimeType combineTimeDifferenceFunction( const std::vector< std::function< TimeType( const TimeType ) > >& timeDifferenceFunctions,
                                        const std::vector< int >& evaluationStepIndices,
                                        const TimeType time )
{
    TimeType timeDifference = TimeType( 0.0 );
    std::vector< TimeType > intermediateTimes;
    intermediateTimes.push_back( time );

    for( unsigned int i = 0; i < timeDifferenceFunctions.size( ); i++ )
    {
        timeDifference += timeDifferenceFunctions.at( i )( intermediateTimes.at( evaluationStepIndices.at( i ) ) );
        intermediateTimes.push_back( time + timeDifference );
    }
    return timeDifference;
}

namespace detail
{

template< typename TimeType >
inline TimeType convertTimeDifferenceFromExtendedTime( const Time& timeDifference )
{
    if constexpr ( std::is_same_v< TimeType, Time > )
    {
        return timeDifference;
    }
    else
    {
        return static_cast< TimeType >( timeDifference.getSeconds< long double >( ) );
    }
}

template< typename TimeType >
inline TimeType convertInterpolatorTimeFromExtendedTime( const Time& currentTime )
{
    if constexpr ( std::is_same_v< TimeType, Time > )
    {
        return currentTime;
    }
    else
    {
        return static_cast< TimeType >( currentTime.getSeconds< long double >( ) );
    }
}

}  // namespace detail

class TimeEphemeris
{
public:

    explicit TimeEphemeris( const std::string& centralBodyName ) : centralBodyName_( centralBodyName )
    { }

    virtual ~TimeEphemeris( ) { }

    //! Function to retrieve the time difference at a given time between two scales.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTime Epoch at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    template< typename TimeType >
    TimeType getTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                                const basic_astrodynamics::TimeScales outputScale,
                                const TimeType inputTime,
                                const std::string& pointIdentifier = "" )
    {
        const std::function< TimeType( const TimeType ) > timeDifferenceFunction =
                getTimeDifferenceFunction< TimeType >( inputScale, outputScale, pointIdentifier );
        return timeDifferenceFunction( inputTime );
    }

    //! Double-precision convenience overload.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTime Epoch at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$ in seconds.
     */
    double getTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                              const basic_astrodynamics::TimeScales outputScale,
                              const double inputTime,
                              const std::string& pointIdentifier = "" )
    {
        return getTimeDifference< double >( inputScale, outputScale, inputTime, pointIdentifier );
    }

    //! Function to retrieve multiple time differences between two scales.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTimes Epochs at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Map of input epoch to time difference \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    template< typename TimeType >
    std::map< TimeType, TimeType > getTimeDifferences( const basic_astrodynamics::TimeScales inputScale,
                                                        const basic_astrodynamics::TimeScales outputScale,
                                                        const std::vector< TimeType >& inputTimes,
                                                        const std::string& pointIdentifier = "" )
    {
        const std::function< TimeType( const TimeType ) > timeDifferenceFunction =
                getTimeDifferenceFunction< TimeType >( inputScale, outputScale, pointIdentifier );

        std::map< TimeType, TimeType > timeDifferences;
        for( const TimeType& inputTime : inputTimes )
        {
            timeDifferences[ inputTime ] = timeDifferenceFunction( inputTime );
        }
        return timeDifferences;
    }

    //! Retrieve the callable that computes the requested conversion at a given epoch.
    /*!
     *  This overload returns a function with templated time input/output type.
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Function object evaluating \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    template< typename TimeType >
    std::function< TimeType( const TimeType ) > getTimeDifferenceFunction(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" )
    {
        const std::function< Time( const Time ) > timeDifferenceFunction =
                getTimeDifferenceFunctionFromExtendedTime( inputScale, outputScale, pointIdentifier );
        return [=]( const TimeType inputTime )
        {
            return detail::convertTimeDifferenceFromExtendedTime< TimeType >( timeDifferenceFunction( Time( inputTime ) ) );
        };
    }

    //! Double-precision convenience overload.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param inputTimes Epochs at which the conversion is evaluated.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Map of input epoch to time difference in seconds.
     */
    std::map< double, double > getTimeDifferences( const basic_astrodynamics::TimeScales inputScale,
                                                   const basic_astrodynamics::TimeScales outputScale,
                                                   const std::vector< double >& inputTimes,
                                                   const std::string& pointIdentifier = "" )
    {
        return getTimeDifferences< double >( inputScale, outputScale, inputTimes, pointIdentifier );
    }

    //! Retrieve the callable that computes the requested conversion at a given epoch.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Function object evaluating \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    std::function< double( const double ) > getTimeDifferenceFunction(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" )
    {
        const std::function< Time( const Time ) > timeDifferenceFunction =
                getTimeDifferenceFunctionFromExtendedTime( inputScale, outputScale, pointIdentifier );
        return [=]( const double inputTime )
        {
            return timeDifferenceFunction( Time( inputTime ) ).getSeconds< double >( );
        };
    }

    //! Retrieve the callable that computes the requested conversion at a given epoch.
    /*!
     *  \param inputScale Input time scale.
     *  \param outputScale Output time scale.
     *  \param pointIdentifier Optional point identifier for topocentric/proper-time conversions.
     *  \return Function object evaluating \f$t_{\mathrm{output}}-t_{\mathrm{input}}\f$.
     */
    virtual std::function< Time( const Time ) > getTimeDifferenceFunctionFromExtendedTime(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) = 0;

protected:
    std::string centralBodyName_;
};

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_H
