#ifndef EARTHSIDEREALTIMECALCULATOR_H
#define EARTHSIDEREALTIMECALCULATOR_H

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

template< typename TimeType >
struct CurrentTimes
{
    CurrentTimes( ){ }

    TimeType tai;
    TimeType tt;
    TimeType tdb;
    TimeType utc;
    TimeType ut1;

    TimeType getTimeValue( basic_astrodynamics::TimeScales requestedScale )
    {
        double valueToReturn = -0.0;
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
            std::cerr<<"Error when getting time value in CurrentTimes, found time scale "<<requestedScale<<std::endl;
        }
        return valueToReturn;
    }
};

class EarthSiderealTimeCalculator
{
public:
    //! Constructor of time scale conversion object.
    /*! Constructor of time scale conversion object. Input values are required for conversions from TT<->UT1.
     *  \param dailyUtcUt1CorrectionInterpolator. Object to interpolate daily corrections between UT1 and UTC,
     *  typically from IERS products.
     *  \param shortPeriodUt1CorrectionCalculator. Object to calculate (sub-)diurnal variations in UT1 not captured by
     *  (typically daily) values given by interpolator.
     *  \param groundStationCartesianPositionFunction. Position in geocentric system where conversion is to be performed.
     *  The TDB<-> TT conversion is sensitive to this at the microSecond level.
     */
    EarthSiderealTimeCalculator(
            const boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > > dailyUtcUt1CorrectionInterpolator,
            const boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator ):
        dailyUtcUt1CorrectionInterpolator_( dailyUtcUt1CorrectionInterpolator ),
        shortPeriodUt1CorrectionCalculator_( shortPeriodUt1CorrectionCalculator ),
        previousEarthFixedPosition_( Eigen::Vector3d::Zero( ) ), updateUt1_( 1 )
    { }

    EarthSiderealTimeCalculator( ):
        dailyUtcUt1CorrectionInterpolator_( boost::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > >( ) ),
        shortPeriodUt1CorrectionCalculator_( getDefaultUT1CorrectionCalculator( ) ),
        previousEarthFixedPosition_( Eigen::Vector3d::Zero( ) ), updateUt1_( 0 )
    { }


    //! Function to convert a time value from the input to the output scale.
    /*!
     *  This function converts a time value from the input to the output scale.
     *  The availabel time scales are defined in the TimeScales enum.
     *  \param inputScale Time scale of inputTimeValue.
     *  \param outputScale Desired time scale for output value.
     *  \param inputTimeValue Time value that is to be converted.
     *  \return Converted time value.
     */
    template< typename TimeType >
    TimeType getCurrentTime( const basic_astrodynamics::TimeScales inputScale, const basic_astrodynamics::TimeScales outputScale, const TimeType& inputTimeValue,
                             const Eigen::Vector3d& earthFixedPosition )
    {
        TimeType convertedTime;
        if( inputScale == outputScale )
        {
            convertedTime = inputTimeValue;
        }
        else
        {
            if( ( static_cast< TimeType >( currentTimes_.getTimeValue( inputScale ) ) !=
                  static_cast< TimeType >( inputTimeValue ) ) ||
                    ( previousEarthFixedPosition_ != earthFixedPosition ) )
            {
                updateTimes( inputScale, inputTimeValue, earthFixedPosition );
            }
            convertedTime = currentTimes_.getTimeValue( outputScale );
        }
        return convertedTime;
    }

    //! Function to recalculate time-values at all time scales from given unput values.
    /*! Function to recalculate time-values at all time scales from given unput values.
    *  \param inputScale Time scale of inputTimeValue.
    *  \param inputTimeValue Time value from which there is to be converted.
    */
    template< typename TimeType >
    void updateTimes( const basic_astrodynamics::TimeScales inputScale, const TimeType& inputTimeValue,
                      const Eigen::Vector3d& earthFixedPosition )
    {
        CurrentTimes< TimeType >& timesToUpdate = getCurrentTimeList< TimeType >( );

        previousEarthFixedPosition_ = earthFixedPosition;
        double siteLongitude = std::atan2( earthFixedPosition.y( ), earthFixedPosition.x( ) );
        double distanceFromSpinAxis = std::sqrt( earthFixedPosition.x( ) * earthFixedPosition.x( ) +
                                                 earthFixedPosition.y( ) * earthFixedPosition.y( ) );
        double distanceFromEquatorialPlane = earthFixedPosition.z( );

        TimeType tdbMinusTt;

        switch( inputScale )
        {
        case basic_astrodynamics::tdb_scale:
            timesToUpdate.tdb = inputTimeValue;
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        inputTimeValue, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tt = timesToUpdate.tdb - tdbMinusTt;
            timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );

            calculateUniversalTimes< TimeType >( );
            break;

        case basic_astrodynamics::tt_scale:
            timesToUpdate.tt = inputTimeValue;
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        inputTimeValue, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;
            timesToUpdate.tai = basic_astrodynamics::convertTTtoTAI< TimeType >( timesToUpdate.tt );

            calculateUniversalTimes< TimeType >( );
            break;

        case basic_astrodynamics::tai_scale:
            timesToUpdate.tai = inputTimeValue;
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;
            calculateUniversalTimes< TimeType >( );

            break;
        case basic_astrodynamics::utc_scale:
            timesToUpdate.utc = inputTimeValue;
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            if( updateUt1_ )
            {
                timesToUpdate.ut1 = dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.utc ) + timesToUpdate.utc;
                timesToUpdate.ut1 += shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tdb );
            }

            break;
        case basic_astrodynamics::ut1_scale:
            if( updateUt1_ == 0 )
            {
                std::cerr<<"Error, requested time conversion from UT1, but cannot update UT1 values"<<std::endl;
            }
            timesToUpdate.ut1 = inputTimeValue;
            timesToUpdate.utc = timesToUpdate.ut1 - dailyUtcUt1CorrectionInterpolator_->interpolate( timesToUpdate.ut1 );
            timesToUpdate.utc -= shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.utc );
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            timesToUpdate.utc = currentTimes_.ut1 - dailyUtcUt1CorrectionInterpolator_->interpolate( currentTimes_.ut1 );
            timesToUpdate.utc -= shortPeriodUt1CorrectionCalculator_->getCorrections( timesToUpdate.tdb );
            timesToUpdate.tai = sofa_interface::convertUTCtoTAI< TimeType >( timesToUpdate.utc );
            timesToUpdate.tt = basic_astrodynamics::convertTAItoTT< TimeType >( timesToUpdate.tai );
            tdbMinusTt = sofa_interface::getTDBminusTT(
                        timesToUpdate.tt, siteLongitude, distanceFromSpinAxis, distanceFromEquatorialPlane );
            timesToUpdate.tdb = timesToUpdate.tt + tdbMinusTt;

            break;

        default:
            std::cerr<<"Currently can only start time converstion at TT or TDB"<<std::endl;
            break;
        }
    }

    void setUpdateUt1( const bool updateUt1 )
    {
        if( updateUt1_ == 1 && ( dailyUtcUt1CorrectionInterpolator_ == NULL || shortPeriodUt1CorrectionCalculator_ == NULL ) )
        {
            std::cerr<<"Error, resetting update ut 1 variable in time converter, but cannot update UT1"<<std::endl;
        }
        updateUt1_ = updateUt1;
    }

private:

    template< typename TimeType >
    CurrentTimes< TimeType >& getCurrentTimeList( );

    template< typename TimeType >
    void calculateUniversalTimes( )
    {
        currentTimes_.utc = sofa_interface::convertTAItoUTC< TimeType >( currentTimes_.tai );

        if( updateUt1_ )
        {
            currentTimes_.ut1 = dailyUtcUt1CorrectionInterpolator_->interpolate( currentTimes_.utc ) + currentTimes_.utc;
            // NOTE: THIS SHOULD ITERATE TO CONVERGENCE!
            currentTimes_.ut1 += shortPeriodUt1CorrectionCalculator_->getCorrections( currentTimes_.tt );
        }
    }

    boost::shared_ptr< interpolators::OneDimensionalInterpolator
    < double, double > > dailyUtcUt1CorrectionInterpolator_;

    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator_;

    CurrentTimes< double > currentTimes_;
    CurrentTimes< Time > currentTimesSplit_;


    Eigen::Vector3d previousEarthFixedPosition_;

    bool updateUt1_;



};

boost::shared_ptr< EarthSiderealTimeCalculator > createDefaultTimeConverter( boost::shared_ptr< EOPReader > eopReader =
        boost::make_shared< EOPReader >( ) );

}

}

#endif // EARTHSIDEREALTIMECALCULATOR_H
