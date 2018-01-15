/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/function.hpp>

#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/eopReader.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace earth_orientation
{

//! Data structure to save the current time in several time scales (TAI, TT, TDB, UTC, UT1)
template< typename TimeType >
struct CurrentTimes
{
    //! Default constructor
    CurrentTimes( ){ }

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
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > >
            dailyUtcUt1CorrectionInterpolator =
            boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > >( ),
            const boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > >
            shortPeriodUt1CorrectionCalculator = getDefaultUT1CorrectionCalculator( ) ):
        dailyUtcUt1CorrectionInterpolator_( dailyUtcUt1CorrectionInterpolator ),
        shortPeriodUt1CorrectionCalculator_( shortPeriodUt1CorrectionCalculator ),
        previousEarthFixedPosition_( Eigen::Vector3d::Zero( ) )
    { }

    //! Function to convert a time value from the input to the output scale.
    /*!
     *  This function converts a time value from the input to the output scale.
     *  The availabel time scales are defined in the TimeScales enum.
     *  \param inputScale Time scale of inputTimeValue.
     *  \param outputScale Desired time scale for output value.
     *  \param inputTimeValue Time value that is to be converted.
     *  \param earthFixedPosition Earth-fixed position at which time conversions are to be evaluated
     *  \return Converted time value.
     */
    template< typename TimeType >
    TimeType getCurrentTime(
            const basic_astrodynamics::TimeScales inputScale, const basic_astrodynamics::TimeScales outputScale,
            const TimeType& inputTimeValue, const Eigen::Vector3d& earthFixedPosition = Eigen::Vector3d::Zero( ) )
    {
        TimeType convertedTime;

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
    void updateTimes( const basic_astrodynamics::TimeScales inputScale, const TimeType& inputTimeValue,
                      const Eigen::Vector3d& earthFixedPosition )
    {
        // Retrieve CurrentTimes object that is to be updated
        CurrentTimes< TimeType >& timesToUpdate = getCurrentTimeList< TimeType >( );

        // Convert position to SOFA input valies
        setCurrentGroundStation< TimeType >( earthFixedPosition );

        double siteLongitude = std::atan2( earthFixedPosition.y( ), earthFixedPosition.x( ) );
        double distanceFromSpinAxis = std::sqrt( earthFixedPosition.x( ) * earthFixedPosition.x( ) +
                                                 earthFixedPosition.y( ) * earthFixedPosition.y( ) );
        double distanceFromEquatorialPlane = earthFixedPosition.z( );

        TimeType tdbMinusTt;

        // Check input type, and call conversion functions accordingly
        switch( inputScale )
        {
        case basic_astrodynamics::tdb_scale:
            timesToUpdate.tdb = inputTimeValue;
            tdbMinusTt = static_cast< TimeType >( sofa_interface::getTDBminusTT(
                        inputTimeValue, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane ) );
            timesToUpdate.tt = timesToUpdate.tdb - tdbMinusTt;
            timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );

            calculateUniversalTimes< TimeType >( );
            break;

        case basic_astrodynamics::tt_scale:
            timesToUpdate.tt = inputTimeValue;
            tdbMinusTt = static_cast< TimeType >( sofa_interface::getTDBminusTT(
                        inputTimeValue, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane ) );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;
            timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );

            calculateUniversalTimes< TimeType >( );
            break;

        case basic_astrodynamics::tai_scale:
            timesToUpdate.tai = inputTimeValue;
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = static_cast< TimeType >( sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane ) );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            calculateUniversalTimes< TimeType >( );

            break;
        case basic_astrodynamics::utc_scale:
            timesToUpdate.utc = inputTimeValue;
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );

            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = static_cast< TimeType >( sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane ) );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            timesToUpdate.ut1 = static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.utc ) )
                    + timesToUpdate.utc;

            timesToUpdate.ut1 += static_cast< TimeType >(
                        shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tdb ) );

            break;
        case basic_astrodynamics::ut1_scale:

            timesToUpdate.ut1 = inputTimeValue;
            timesToUpdate.utc = timesToUpdate.ut1 -
                    static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.ut1 ) );
            timesToUpdate.utc -=
                    static_cast< TimeType >( shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.utc ) );
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = static_cast< TimeType >( sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane ) );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            // Iterate conversion.
            timesToUpdate.utc = timesToUpdate.ut1 -
                    static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.utc ) );
            timesToUpdate.utc -= static_cast< TimeType >(
                        shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tt ) );
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            break;

        default:
            throw std::runtime_error( "Error when performing Earth time scales, input time not recognized" );
            break;
        }
    }

    template< typename TimeType >
    double getUt1Correction(
            const basic_astrodynamics::TimeScales inputScale, const TimeType& inputTimeValue,
            const Eigen::Vector3d currentPosition = Eigen::Vector3d::Zero( ) )
    {
        TimeType currentUtc = getCurrentTime( inputScale, basic_astrodynamics::utc_scale, inputTimeValue, currentPosition );
        TimeType currentTt = getCurrentTime( inputScale, basic_astrodynamics::tt_scale, inputTimeValue, currentPosition );

        return dailyUtcUt1CorrectionInterpolator_->interpolate( currentUtc ) +
                shortPeriodUt1CorrectionCalculator_->getCorrections( currentTt );

    }

    //! Interpolator for UT1 corrections, values published daily by IERS
    boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > > getDailyUtcUt1CorrectionInterpolator( )
    {
        return dailyUtcUt1CorrectionInterpolator_;
    }

private:

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
        getCurrentTimeList< TimeType >( ).utc = sofa_interface::convertTAItoUTC< TimeType >(
                    getCurrentTimeList< TimeType >( ).tai );

        getCurrentTimeList< TimeType >( ).ut1 = static_cast< TimeType >( dailyUtcUt1CorrectionInterpolator_->interpolate(
                    getCurrentTimeList< TimeType >( ).utc ) ) + getCurrentTimeList< TimeType >( ).utc;

        getCurrentTimeList< TimeType >( ).ut1 +=
                static_cast< TimeType >( shortPeriodUt1CorrectionCalculator_->getCorrections(
                    getCurrentTimeList< TimeType >( ).tt ) );
    }

    //! Interpolator for UT1 corrections, values published daily by IERS
    boost::shared_ptr< interpolators::OneDimensionalInterpolator
    < double, double > > dailyUtcUt1CorrectionInterpolator_;

    //! Object to compute the short-period variations in UT1
    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator_;

    //! Object containing current times, as set by last updateTimes< double > function.
    CurrentTimes< double > currentTimes_;

    //! Object containing current times, as set by last updateTimes< Time > function.
    CurrentTimes< Time > currentTimesSplit_;

    //! Value of ground station position used on last call to updateTimes< double > function
    Eigen::Vector3d previousEarthFixedPosition_;

    //! Value of ground station position used on last call to updateTimes< Time > function
    Eigen::Vector3d previousEarthFixedPositionSplit_;
};

//! Function to create the default Earth time scales conversion object
/*!
 * Function to create the default Earth time scales conversion object. All (sub-)diurnal corrections to UTC-UT1
 * according to IERS 2010, and UT1 daily corrections published by IERS (read from input EOPReader).
 * \param eopReader Object that reads an Earth Orientation Parameters file.
 * \return Default Earth time scales conversion object
 */
boost::shared_ptr< TerrestrialTimeScaleConverter >  createDefaultTimeConverter( const boost::shared_ptr< EOPReader > eopReader =
        boost::make_shared< EOPReader >( ) );

}

}

#endif // TUDAT_TERRESTRIALTIMESCALECONVERTER_H
