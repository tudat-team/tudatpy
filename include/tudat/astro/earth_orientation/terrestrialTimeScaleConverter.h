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

#ifndef TUDAT_TERRESTRIALTIMESCALECONVERTER_H
#define TUDAT_TERRESTRIALTIMESCALECONVERTER_H

#include <functional>
#include <optional>
#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/earth_orientation/eopReader.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace earth_orientation
{

//! UT1 and its derivative with respect to TAI.
struct Ut1ValueAndRate {
    //! UT1 (seconds since J2000).
    Time ut1;

    //! Derivative of UT1 with respect to TAI; NaN if it was not computed.
    double dUt1DTai = TUDAT_NAN;
};

//! Class used to convert between the terrestrial time scales TAI, TT, TDB, UTC and UT1.
/*!
 *  Every conversion is done in two steps: from the input time scale to TAI, and from TAI to the output time scale. TAI is
 *  used as the intermediate time scale because, unlike UTC, it has no jumps at leap seconds. UT1 (and TDB at a point away
 *  from the geocentre) is computed from the Earth orientation data given to the constructor. Without these data, conversions
 *  to and from UT1 are only possible before 1962, using historical values of Delta T. Before 1962, UTC is taken to be equal
 *  to UT1.
 */
class TerrestrialTimeScaleConverter
{
public:
    //! Constructor of time scale conversion object.
    /*!
     *  Constructor of time scale conversion object.
     *  \param eopData Interpolated daily Earth orientation data, used to compute UT1. If nullptr, UT1 is only available
     *  before 1962.
     *  \param shortPeriodUt1CorrectionCalculator Object to calculate (sub-)diurnal variations in UT1 not captured by the daily
     *  values. If nullptr, these variations are not included.
     *  \param tdbToTtInterpolatorSettings Settings to interpolate TDB - TT at the geocentre, instead of computing it at each
     *  epoch. If nullptr, TDB - TT is computed at each epoch.
     */
    TerrestrialTimeScaleConverter(
            const std::shared_ptr< EarthOrientationDataEvaluator > eopData = nullptr,
            const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator =
                    getDefaultUT1CorrectionCalculator( ),
            const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > tdbToTtInterpolatorSettings = nullptr ):
        eopData_( eopData ), shortPeriodUt1CorrectionCalculator_( shortPeriodUt1CorrectionCalculator )
    {
        if( tdbToTtInterpolatorSettings != nullptr )
        {
            std::function< double( const double ) > correctionFunction =
                    std::bind< double( const double, const double, const double, const double ) >(
                            sofa_interface::getTDBminusTT, std::placeholders::_1, 0.0, 0.0, 0.0 );

            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > correctionInterpolator =
                    interpolators::createOneDimensionalInterpolator< double, double >( correctionFunction, tdbToTtInterpolatorSettings );
            tdbToTtInterpolators_[ std::make_tuple( 0.0, 0.0, 0.0 ) ] = correctionInterpolator;
        }

        // Definition of epoch at which UTC was introduced
        utcIntroductionEpochInTai_ = ( basic_astrodynamics::JULIAN_DAY_OF_EOP_INTRODUCTION - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY;
    }

    //! Function to convert a time value from the input to the output scale.
    /*!
     *  This function converts a time value from the input to the output scale.
     *  The available time scales are defined in the TimeScales enum.
     *  \param inputScale Time scale of inputTimeValue.
     *  \param outputScale Desired time scale for output value.
     *  \param inputTimeValue Time value that is to be converted.
     *  \param earthFixedPosition Earth-fixed position at which time conversions are to be evaluated
     *  \return Converted time value.
     */
    template< typename TimeType >
    TimeType getCurrentTime( const basic_astrodynamics::TimeScales inputScale,
                             const basic_astrodynamics::TimeScales outputScale,
                             const TimeType& inputTimeValue,
                             const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
    {
        if( inputScale == outputScale )
        {
            checkTimeScale( inputScale );
            return inputTimeValue;
        }
        return static_cast< TimeType >(
                getTimeFromTai( outputScale, getTaiTime( inputScale, Time( inputTimeValue ), earthFixedPosition ), earthFixedPosition ) );
    }

    //! Function to compute the difference between a time value in the output scale and in the input scale.
    /*!
     *  Function to compute the difference between a time value in the output scale and in the input scale, i.e. the converted
     *  time value minus the input time value.
     *  \param inputScale Time scale of inputTimeValue.
     *  \param outputScale Time scale to which inputTimeValue is converted.
     *  \param inputTimeValue Time value that is to be converted.
     *  \param earthFixedPosition Earth-fixed position at which time conversions are to be evaluated
     *  \return Converted time value minus input time value.
     */
    template< typename TimeType >
    TimeType getCurrentTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                                       const basic_astrodynamics::TimeScales outputScale,
                                       const TimeType& inputTimeValue,
                                       const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
    {
        if( inputScale == outputScale )
        {
            checkTimeScale( inputScale );
            return TimeType( 0.0 );
        }
        const Time input( inputTimeValue );
        return static_cast< TimeType >(
                getTimeFromTai( outputScale, getTaiTime( inputScale, input, earthFixedPosition ), earthFixedPosition ) - input );
    }

    //! Function to convert a vector of time values from the input to the output scale.
    /*!
     *  This function converts a vector of time values from the input to the output scale.
     *  The available time scales are defined in the TimeScales enum.
     *  \param inputScale Time scale of inputTimeValues.
     *  \param outputScale Desired time scale for output values.
     *  \param inputTimeValues Time values that are to be converted.
     *  \param earthFixedPosition Earth-fixed position at which time conversions are to be evaluated
     *  \return Converted time values.
     */
    template< typename TimeType >
    std::vector< TimeType > getCurrentTimes( const basic_astrodynamics::TimeScales inputScale,
                                             const basic_astrodynamics::TimeScales outputScale,
                                             const std::vector< TimeType >& inputTimeValues,
                                             const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
    {
        std::vector< TimeType > convertedTimes;

        for( unsigned int i = 0; i < inputTimeValues.size( ); ++i )
        {
            convertedTimes.push_back( getCurrentTime( inputScale, outputScale, inputTimeValues.at( i ), earthFixedPosition ) );
        }

        return convertedTimes;
    }

    //! Function to convert a vector of time values from the input to the output scale.
    /*!
     *  This function converts a vector of time values from the input to the output scale.
     *  The available time scales are defined in the TimeScales enum.
     *  \param inputScale Time scale of inputTimeValues.
     *  \param outputScale Desired time scale for output values.
     *  \param inputTimeValues Time values that are to be converted.
     *  \param earthFixedPositions Earth-fixed positions at which time conversions are to be evaluated
     *  \return Converted time values.
     */
    template< typename TimeType >
    std::vector< TimeType > getCurrentTimes( const basic_astrodynamics::TimeScales inputScale,
                                             const basic_astrodynamics::TimeScales outputScale,
                                             const std::vector< TimeType >& inputTimeValues,
                                             const std::vector< Eigen::Vector3d >& earthFixedPositions )
    {
        if( inputTimeValues.size( ) != earthFixedPositions.size( ) )
        {
            throw std::runtime_error(
                    "Error time values between scales: number of inputted time values and number of Earth-fixed positions are not "
                    "consistent." );
        }

        std::vector< TimeType > convertedTimes;

        for( unsigned int i = 0; i < inputTimeValues.size( ); ++i )
        {
            convertedTimes.push_back( getCurrentTime( inputScale, outputScale, inputTimeValues.at( i ), earthFixedPositions.at( i ) ) );
        }
        return convertedTimes;
    }

    template< typename TimeType >
    std::vector< TimeType > getCurrentTimesFromSinglePosition( const basic_astrodynamics::TimeScales inputScale,
                                                               const basic_astrodynamics::TimeScales outputScale,
                                                               const std::vector< TimeType >& inputTimeValues,
                                                               const Eigen::Vector3d& earthFixedPosition )
    {
        std::vector< TimeType > convertedTimes;
        convertedTimes.resize( inputTimeValues.size( ) );

        for( unsigned int i = 0; i < inputTimeValues.size( ); ++i )
        {
            convertedTimes[ i ] = getCurrentTime( inputScale, outputScale, inputTimeValues.at( i ), earthFixedPosition );
        }
        return convertedTimes;
    }

    //! Function to compute UT1 - UTC (in seconds) at a given epoch, including short-period variations.
    template< typename TimeType >
    double getUt1Correction( const basic_astrodynamics::TimeScales inputScale,
                             const TimeType& inputTimeValue,
                             const Eigen::Vector3d currentPosition = Eigen::Vector3d::Zero( ) )
    {
        const Time tai = getTaiTime( inputScale, Time( inputTimeValue ), currentPosition );
        return static_cast< double >( getUt1FromTai( tai ) - getTimeFromTai( basic_astrodynamics::utc_scale, tai, currentPosition ) );
    }

    //! Function to retrieve the interpolated daily Earth orientation data (nullptr if not provided).
    std::shared_ptr< EarthOrientationDataEvaluator > getEarthOrientationDataEvaluator( ) const
    {
        return eopData_;
    }

    //! Function to retrieve the object computing the short-period variations in UT1 (nullptr if not used).
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > getShortPeriodUt1CorrectionCalculator( ) const
    {
        return shortPeriodUt1CorrectionCalculator_;
    }

    //! UT1 (seconds since J2000) at a given TAI, from the daily data plus short-period variations.
    Time getUt1FromTai( const Time& tai );

    //! UT1 (seconds since J2000) at a given TAI, and its derivative with respect to TAI. The derivative is that of the
    //! interpolated daily data plus that of the short-period variations; it is not available before the start of the data.
    Ut1ValueAndRate getUt1FromTaiWithRate( const Time& tai );

    //! Derivative of TT with respect to TDB at the geocentre, at a given TDB (seconds since J2000). When TDB - TT is
    //! interpolated, this is the derivative of the interpolated value.
    double getGeocentricTtRateWithRespectToTdb( const Time& tdb );

    //! Convert an epoch in the given time scale to TAI (seconds since J2000).
    Time getTaiTime( const basic_astrodynamics::TimeScales inputScale,
                     const Time& input,
                     const Eigen::Vector3d& position = Eigen::Vector3d::Zero( ) );

private:
    //! Throw an exception if the time scale is not one of TAI, TT, TDB, UTC and UT1.
    void checkTimeScale( const basic_astrodynamics::TimeScales timeScale );

    //! Convert TAI (seconds since J2000) to the given time scale.
    Time getTimeFromTai( const basic_astrodynamics::TimeScales outputScale, const Time& tai, const Eigen::Vector3d& position );

    //! Load the table of historical values of Delta T (TT - UT1), when it is first needed.
    void loadHistoricalDeltaTInterpolator( );

    //! UT1 from historical values of Delta T, at a given TAI (before 1962).
    Time getHistoricalUt1FromTai( const Time& tai );

    //! TAI from historical values of Delta T, at a given UT1 (before 1962).
    Time getTaiFromHistoricalUt1( const Time& ut1 );

    //! Evaluate UT1 and, only when requested, its derivative with respect to TAI.
    Ut1ValueAndRate evaluateUt1FromTai( const Time& tai, bool requireRate );

    //! TAI at which UT1 (or TDB at a point away from the geocentre) has the requested value, within the Earth orientation data.
    Time invertTimeFromTai( basic_astrodynamics::TimeScales inputScale, const Time& requested, const Eigen::Vector3d& position );

    //! TDB - TT (seconds) at a given TT and Earth-fixed position.
    double getTdbMinusTt( const Time& tt, const Eigen::Vector3d& position );

    //! TT - TDB (seconds) at the geocentre, at a given TDB (seconds since J2000).
    double getGeocentricTtMinusTdb( const Time& tdb );

    //! Interpolated daily Earth orientation data; nullptr if not provided.
    std::shared_ptr< EarthOrientationDataEvaluator > eopData_;

    //! Object to compute the short-period variations in UT1
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator_;

    //! Interpolators of TDB - TT, per Earth-fixed position (only the geocentre is used).
    std::map< std::tuple< double, double, double >, std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > >
            tdbToTtInterpolators_;

    //! Epoch (TAI seconds since J2000) at which UTC was introduced; before this epoch, UTC is taken to be equal to UT1.
    double utcIntroductionEpochInTai_;

    //! Interpolator of historical values of Delta T (TT - UT1), as a function of the year.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > historicalDeltaTInterpolator_;
};

//! Function to create the default Earth time scales conversion object
/*!
 * Function to create the default Earth time scales conversion object. All (sub-)diurnal corrections to UTC-UT1
 * according to IERS 2010, and UT1 daily corrections published by IERS (read from input EOPReader).
 * \param eopReader Object that reads an Earth Orientation Parameters file.
 * \param reuseExistingConverter Reuse a previously created default converter. This option is only used when eopReader is nullptr.
 * \return Default Earth time scales conversion object
 */
std::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter(
        const std::shared_ptr< EOPReader > eopReader = nullptr,
        const bool reuseExistingConverter = true,
        const std::optional< EarthOrientationEvaluationSettings >& settings = std::nullopt );

// static const std::shared_ptr< TerrestrialTimeScaleConverter > defaultTimeConverter = createDefaultTimeConverter( );

template< typename TimeType >
TimeType convertTimeScale( const TimeType& inputTime,
                           const basic_astrodynamics::TimeScales inputScale,
                           const basic_astrodynamics::TimeScales outputScale,
                           const std::shared_ptr< TerrestrialTimeScaleConverter > timeConverter = createDefaultTimeConverter( ),
                           const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
{
    Time time = Time( inputTime );
    return static_cast< TimeType >( timeConverter->getCurrentTime( inputScale, outputScale, time, earthFixedPosition ) );
}

template< typename TimeType >
TimeType convertToTimeScale( const basic_astrodynamics::DateTime inputDateTime,
                             const basic_astrodynamics::TimeScales inputScale,
                             const basic_astrodynamics::TimeScales outputScale,
                             const std::shared_ptr< TerrestrialTimeScaleConverter > timeConverter = createDefaultTimeConverter( ),
                             const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
{
    Time time = inputDateTime.epoch< Time >( );
    return static_cast< TimeType >( timeConverter->getCurrentTime( inputScale, outputScale, time, earthFixedPosition ) );
}

}  // namespace earth_orientation

}  // namespace tudat

#endif  // TUDAT_TERRESTRIALTIMESCALECONVERTER_H
