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

//! Data structure to save the current time in several time scales (TAI, TT, TDB, UTC, UT1)
template< typename TimeType >
struct CurrentTimes {
    //! Default constructor
    CurrentTimes( ): tai( TUDAT_NAN ), tt( TUDAT_NAN ), tdb( TUDAT_NAN ), utc( TUDAT_NAN ), ut1( TUDAT_NAN ) { }

    //! Function to retrieve the current time in requested scale
    /*!
     * Function to retrieve the current time in requested scale
     * \param requestedScale Time scale for which time is to be returned
     * \return Current time in requested scale.
     */
    TimeType getTimeValue( basic_astrodynamics::TimeScales requestedScale )
    {
        TimeType valueToReturn = -0.0;
        switch( requestedScale )
        {
            case basic_astrodynamics::tai_scale:
                valueToReturn = tai;
                break;
            case basic_astrodynamics::tt_scale:
                valueToReturn = tt;
                break;
            case basic_astrodynamics::tdb_scale:
                valueToReturn = tdb;
                break;
            case basic_astrodynamics::utc_scale:
                valueToReturn = utc;
                break;
            case basic_astrodynamics::ut1_scale:
                valueToReturn = ut1;
                break;
            default:
                std::cerr << "Error when getting time value in CurrentTimes, found time scale " << requestedScale << std::endl;
        }
        return valueToReturn;
    }

    //! Current time in TAI.
    TimeType tai;

    //! Current time in TT.
    TimeType tt;

    //! Current time in TDB.
    TimeType tdb;

    //! Current time in UTC.
    TimeType utc;

    //! Current time in UT1.
    TimeType ut1;
};

//! Class used to convert between terrestrial time scales TAI, TT, TDB, UTC ans UT1
class TerrestrialTimeScaleConverter
{
public:
    //! Constructor of time scale conversion object.
    /*!
     *  Constructor of time scale conversion object. Input values are required for to/from UT1.
     *  \param dailyUtcUt1CorrectionInterpolator Object to interpolate daily corrections between UT1 and UTC,
     *  typically from IERS products.
     *  \param shortPeriodUt1CorrectionCalculator Object to calculate (sub-)diurnal variations in UT1 not captured by
     *  (typically daily) values given by interpolator.
     */
    TerrestrialTimeScaleConverter(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > dailyUtcUt1CorrectionInterpolator =
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > >( ),
            const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator =
                    getDefaultUT1CorrectionCalculator( ),
            const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > tdbToTtInterpolatorSettings = nullptr ):
        dailyUtcUt1CorrectionInterpolator_( dailyUtcUt1CorrectionInterpolator ),
        shortPeriodUt1CorrectionCalculator_( shortPeriodUt1CorrectionCalculator ), previousEarthFixedPosition_( Eigen::Vector3d::Zero( ) )
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

        // Definition of epoch at which UTC was introduced (with 2 hour buffer time)
        utcIntroductionEpochInTai_ = ( basic_astrodynamics::JULIAN_DAY_OF_EOP_INTRODUCTION - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY;

        // Read historical Delta T values
        std::map< double, double > historicalDeltaTMap = input_output::readFloatingPointMapFromFile< double, double >(
                paths::getEarthOrientationDataFilesPath( ) + "/historicalDeltaT.txt", ",", "#" );

        historicalDeltaTInterpolator_ =
                interpolators::createOneDimensionalInterpolator( historicalDeltaTMap, interpolators::cubicSplineInterpolation( ) );
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
        TimeType convertedTime;

        // Check whether desired conversion is possible (conversion to/from UT1 requires interpolator)
        if( inputScale == basic_astrodynamics::ut1_scale || outputScale == basic_astrodynamics::ut1_scale )
        {
            if( dailyUtcUt1CorrectionInterpolator_ == nullptr )
            {
                throw std::runtime_error(
                        "Error when converting to/from UT1 time scale: UTC to UT1 interpolator "
                        "was not provided." );
            }
        }

        // Check if any conversion should take place.
        if( inputScale == outputScale )
        {
            convertedTime = inputTimeValue;
        }
        else
        {
            // Check if update is required
            if( !( static_cast< TimeType >( getCurrentTimeList< TimeType >( ).getTimeValue( inputScale ) ) ==
                   static_cast< TimeType >( inputTimeValue ) ) ||
                !( getPreviousGroundStationPosition< TimeType >( ) == earthFixedPosition ) )
            {
                updateTimes< TimeType >( inputScale, inputTimeValue, earthFixedPosition );
            }
            convertedTime = getCurrentTimeList< TimeType >( ).getTimeValue( outputScale );
        }
        return convertedTime;
    }

    template< typename TimeType >
    TimeType getCurrentTimeDifference( const basic_astrodynamics::TimeScales inputScale,
                                       const basic_astrodynamics::TimeScales outputScale,
                                       const TimeType& inputTimeValue,
                                       const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
    {
        Time convertedTime = getCurrentTime< Time >( inputScale, outputScale, Time( inputTimeValue ), earthFixedPosition );
        return static_cast< TimeType >( convertedTime - inputTimeValue );
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

    //! Function to reset all current times at given precision to NaN.
    template< typename TimeType >
    void resetTimes( )
    {
        CurrentTimes< TimeType >& timesToUpdate = getCurrentTimeList< TimeType >( );
        timesToUpdate.tai = TUDAT_NAN;
        timesToUpdate.tt = TUDAT_NAN;
        timesToUpdate.tdb = TUDAT_NAN;
        timesToUpdate.ut1 = TUDAT_NAN;
        timesToUpdate.utc = TUDAT_NAN;
    }

    //! Function to recalculate time-values at all time scales from given unput values.
    /*!
     * Function to recalculate time-values at all time scales from given unput values.
     *  \param inputScale Time scale of inputTimeValue.
     *  \param inputTimeValue Time value from which there is to be converted.
     *  \param earthFixedPosition Earth-fixed position at which time conversions are to be evaluated
     */
    template< typename TimeType >
    void updateTimes( const basic_astrodynamics::TimeScales inputScale,
                      const TimeType& inputTimeValue,
                      const Eigen::Vector3d& earthFixedPosition )
    {
        // Retrieve CurrentTimes object that is to be updated
        CurrentTimes< TimeType >& timesToUpdate = getCurrentTimeList< TimeType >( );

        // Convert position to SOFA input valies
        setCurrentGroundStation< TimeType >( earthFixedPosition );

        TimeType tdbMinusTt;

        // Check input type, and call conversion functions accordingly
        switch( inputScale )
        {
            case basic_astrodynamics::tdb_scale:
                timesToUpdate.tdb = inputTimeValue;
                tdbMinusTt = static_cast< TimeType >( this->getTDBminusTT( inputTimeValue, earthFixedPosition ) );
                timesToUpdate.tt = timesToUpdate.tdb - tdbMinusTt;
                timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );
                calculateUniversalTimes< TimeType >( );
                break;

            case basic_astrodynamics::tt_scale:
                timesToUpdate.tt = inputTimeValue;
                tdbMinusTt = static_cast< TimeType >( this->getTDBminusTT( inputTimeValue, earthFixedPosition ) );
                timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;
                timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );

                calculateUniversalTimes< TimeType >( );
                break;

            case basic_astrodynamics::tai_scale:
                timesToUpdate.tai = inputTimeValue;
                timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
                tdbMinusTt = static_cast< TimeType >( this->getTDBminusTT( timesToUpdate.tt, earthFixedPosition ) );
                timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

                calculateUniversalTimes< TimeType >( );

                break;
            case basic_astrodynamics::utc_scale:
                // --- 1. Use function to calculate TAI and TT from UTC ---
                calculateAtomicTimesFromUtc< TimeType >( inputTimeValue );

                // --- 2. Calculate TDB from the newly computed TT ---
                tdbMinusTt = static_cast< TimeType >( this->getTDBminusTT( timesToUpdate.tt, earthFixedPosition ) );
                timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

                // --- 3. Calculate UT1 if possible (for modern era) ---
                if( dailyUtcUt1CorrectionInterpolator_ != nullptr && timesToUpdate.tai >= utcIntroductionEpochInTai_ )
                {
                    timesToUpdate.ut1 = static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.utc ) ) +
                            timesToUpdate.utc;
                    timesToUpdate.ut1 +=
                            static_cast< TimeType >( shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tdb ) );
                }
                break;
            case basic_astrodynamics::ut1_scale:

                timesToUpdate.ut1 = inputTimeValue;
                try
                {
                    timesToUpdate.utc = timesToUpdate.ut1 - static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.ut1 ) );
                }
                catch( std::runtime_error& caughtException )
                {
                    throw std::runtime_error( "Error in UTC-UT1 correction.\nOriginal error: " + std::string( caughtException.what( ) ) );
                }

                timesToUpdate.utc -= static_cast< TimeType >( shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.utc ) );
                timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
                timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
                tdbMinusTt = static_cast< TimeType >( this->getTDBminusTT( timesToUpdate.tt, earthFixedPosition ) );
                timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

                // Iterate conversion.
                try
                {
                    timesToUpdate.utc =
                        timesToUpdate.ut1 - static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.utc ) );
                }
                catch( std::runtime_error& caughtException )
                {
                    throw std::runtime_error( "Error in UTC-UT1 correction.\nOriginal error: " + std::string( caughtException.what( ) ) );
                }

                timesToUpdate.utc -= static_cast< TimeType >( shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tt ) );
                timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
                timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
                tdbMinusTt = this->getTDBminusTT( timesToUpdate.tt, earthFixedPosition );
                timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

                break;

            default:
                throw std::runtime_error( "Error when performing Earth time scales, input time not recognized" );
                break;
        }
    }

    template< typename TimeType >
    double getUt1Correction( const basic_astrodynamics::TimeScales inputScale,
                             const TimeType& inputTimeValue,
                             const Eigen::Vector3d currentPosition = Eigen::Vector3d::Zero( ) )
    {
        TimeType currentUtc = getCurrentTime( inputScale, basic_astrodynamics::utc_scale, inputTimeValue, currentPosition );
        TimeType currentTt = getCurrentTime( inputScale, basic_astrodynamics::tt_scale, inputTimeValue, currentPosition );

        if( dailyUtcUt1CorrectionInterpolator_ == nullptr )
        {
            throw std::runtime_error( "Error when getting UT1 correction: UTC to UT1 interpolator was not provided." );
        }

        try
        {
            return dailyUtcUt1CorrectionInterpolator_->interpolate( currentUtc ) +
                shortPeriodUt1CorrectionCalculator_->getCorrections( currentTt );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error in UT1 correction.\nOriginal error: " + std::string( caughtException.what( ) ) );
        }


    }

    //! Interpolator for UT1 corrections, values published daily by IERS
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > getDailyUtcUt1CorrectionInterpolator( )
    {
        return dailyUtcUt1CorrectionInterpolator_;
    }

private:
    double getTDBminusTT( const double ttOrTdbSinceJ2000, const Eigen::Vector3d earthFixedPosition )
    {
        auto tupleToCheck = std::make_tuple( earthFixedPosition[ 0 ], earthFixedPosition[ 1 ], earthFixedPosition[ 2 ] );
        if( tdbToTtInterpolators_.count( tupleToCheck ) != 0 )
        {
            try
            {
                return tdbToTtInterpolators_.at( tupleToCheck )->interpolate( ttOrTdbSinceJ2000 );
            }
            catch( std::runtime_error& caughtException )
            {
                throw std::runtime_error( "Error in TDB-TT interpolation.\nOriginal error: " + std::string( caughtException.what( ) ) );
            }


        }
        else
        {
            return sofa_interface::getTDBminusTT( ttOrTdbSinceJ2000, earthFixedPosition );
        }
    }

    //! Function to get current time list at requested numerical precision
    /*!
     *  Function to get current time list at requested numerical precision
     *  \return Current time list at requested numerical precision
     */
    template< typename TimeType >
    CurrentTimes< TimeType >& getCurrentTimeList( );

    //! Function to get ground station position used on last call to updateTimes function at requested numerical precision
    /*!
     * Function to get ground station position used on last call to updateTimes function at requested numerical precision
     *  \return Ground station position used on last call to updateTimes function at requested numerical precision
     */
    template< typename TimeType >
    Eigen::Vector3d& getPreviousGroundStationPosition( );

    //! Function to reset ground station position used on last call to updateTimes function at requested numerical precision
    /*!
     * Function to reset ground station position used on last call to updateTimes function at requested numerical precision
     *  \return Ground station position used on last call to updateTimes function at requested numerical precision
     */
    template< typename TimeType >
    void setCurrentGroundStation( const Eigen::Vector3d& currentGroundStation );

    //! Function to update the universal times (UT1 and UTC) in CurrentTimes member at requested precision
    template< typename TimeType >
    void calculateUniversalTimes( )
    {
        if( getCurrentTimeList< TimeType >( ).tai < utcIntroductionEpochInTai_ )
        {
            double approximateYear = getCurrentTimeList< TimeType >( ).tt / physical_constants::JULIAN_YEAR + 2000.0;
            try
            {
                getCurrentTimeList< TimeType >( ).ut1 =
                    getCurrentTimeList< TimeType >( ).tt - historicalDeltaTInterpolator_->interpolate( approximateYear );
            }
            catch( std::runtime_error& caughtException )
            {
                throw std::runtime_error( "Error in historical UT calculation.\nOriginal error: " + std::string( caughtException.what( ) ) );
            }


            getCurrentTimeList< TimeType >( ).utc = getCurrentTimeList< TimeType >( ).ut1;
        }
        else
        {
            getCurrentTimeList< TimeType >( ).utc = sofa_interface::convertTAItoUTC< TimeType >( getCurrentTimeList< TimeType >( ).tai );

            if( dailyUtcUt1CorrectionInterpolator_ != nullptr )
            {
                try
                {
                    getCurrentTimeList< TimeType >( ).ut1 = static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate(
                                                                    getCurrentTimeList< TimeType >( ).utc ) ) + getCurrentTimeList< TimeType >( ).utc;
                }
                catch( std::runtime_error& caughtException )
                {
                    throw std::runtime_error( "Error in UT1 calculation.\nOriginal error: " + std::string( caughtException.what( ) ) );
                }


                getCurrentTimeList< TimeType >( ).ut1 += static_cast< TimeType >(
                        shortPeriodUt1CorrectionCalculator_->getCorrections( getCurrentTimeList< TimeType >( ).tt ) );
            }
        }
    }


    template< typename TimeType >
    void calculateAtomicTimesFromUtc( const TimeType& inputUtcTime )
    {
        // --- NEW LOGIC: Check for historical date FIRST ---
        // The UTC introduction epoch is in TAI, but we can use it as an approximate
        // threshold for the input time (which is UTC-like) to decide which path to take.
        // This avoids calling SOFA functions with dates they cannot handle.
        if ( inputUtcTime < utcIntroductionEpochInTai_)
        {
            // --- 1. Historical Time (Pre-UTC Era) ---
            // This path is now correctly taken for the 1893 date.
            // In this period, the input "UTC" is effectively UT1.
            // The conversion path is UT1 -> TT -> TAI.
            getCurrentTimeList< TimeType >( ).ut1 = inputUtcTime;
            getCurrentTimeList< TimeType >( ).utc = inputUtcTime;

            // Approximate the year to look up the historical ΔT value.
            double approximateYear = ( inputUtcTime / physical_constants::JULIAN_YEAR ) + 2000.0;

            // Calculate TT from UT1 using the historical ΔT. (TT = UT1 + ΔT)
            getCurrentTimeList< TimeType >( ).tt = getCurrentTimeList< TimeType >( ).ut1 +
                historicalDeltaTInterpolator_->interpolate( approximateYear );

            // Convert TT to TAI using the standard constant offset.
            getCurrentTimeList< TimeType >( ).tai = basic_astrodynamics::convertTTtoTAI< TimeType >(
                getCurrentTimeList< TimeType >( ).tt );
        }
        else
        {
            // --- 2. Modern Time (Post-UTC Era) ---
            // The conversion is handled directly by the SOFA function, which correctly
            // accounts for all leap seconds.
            getCurrentTimeList< TimeType >( ).utc = inputUtcTime;
            getCurrentTimeList< TimeType >( ).tai = sofa_interface::convertUTCtoTAI< TimeType >(
                getCurrentTimeList< TimeType >( ).utc );

            // >>>>> THE FIX IS HERE <<<<<
            // We must also convert the newly calculated TAI to TT.
            getCurrentTimeList< TimeType >( ).tt = basic_astrodynamics::convertTAItoTT< TimeType >(
                getCurrentTimeList< TimeType >( ).tai );
        }
    }

    //! Interpolator for UT1 corrections, values published daily by IERS
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > dailyUtcUt1CorrectionInterpolator_;

    //! Object to compute the short-period variations in UT1
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator_;

    //! Object containing current times, as set by last updateTimes< double > function.
    CurrentTimes< double > currentTimes_;

    //! Object containing current times, as set by last updateTimes< Time > function.
    CurrentTimes< Time > currentTimesSplit_;

    //! Value of ground station position used on last call to updateTimes< double > function
    Eigen::Vector3d previousEarthFixedPosition_;

    //! Value of ground station position used on last call to updateTimes< Time > function
    Eigen::Vector3d previousEarthFixedPositionSplit_;

    std::map< std::tuple< double, double, double >, std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > >
            tdbToTtInterpolators_;

    double utcIntroductionEpochInTai_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > historicalDeltaTInterpolator_;
};

//! Function to create the default Earth time scales conversion object
/*!
 * Function to create the default Earth time scales conversion object. All (sub-)diurnal corrections to UTC-UT1
 * according to IERS 2010, and UT1 daily corrections published by IERS (read from input EOPReader).
 * \param eopReader Object that reads an Earth Orientation Parameters file.
 * \return Default Earth time scales conversion object
 */
std::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter(
        const std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >( ) );

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
