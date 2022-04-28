/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Data file:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *                        http://celestrak.com/SpaceData/sw20110101.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#ifndef TUDAT_SOLAR_ACTIVITY_DATA_H
#define TUDAT_SOLAR_ACTIVITY_DATA_H

#include <string>
#include <map>

#include <Eigen/Core>

#include <memory>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{
namespace input_output
{
namespace solar_activity
{

//! Solar activity data class.
/*!
 * Class containing all variables provided in the http://celestrak.com/SpaceData/sw19571001.txt
 * file, describing "space weather" on each day since October 1957. See reference for details on
 * variables.
 */
struct SolarActivityData
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SolarActivityData( );

    //! Year.
    unsigned int year;

    //! Month.
    unsigned int month;

    //! Day.
    unsigned int day;

    //! Bartels Solar Rotation Number.
    unsigned int bartelsSolarRotationNumber;

    //! Day within the Bartels 27-day cycle (01-27).
    unsigned int dayOfBartelsCycle;

    //! Sum of the 8 planetary range indices.
    /*! Sum of the 8 planetary range indices (Kp) for the day expressed to the nearest third of a
    *   unit.
    */
    unsigned int planetaryRangeIndexSum;

    //! Arithmetic average of the 8 planetary equivalent amplitudes (Ap) indices for the day.
    unsigned int planetaryEquivalentAmplitudeAverage;

    //! Cp or Planetary Daily Character Figure.
    /*! A qualitative estimate of overall level of magnetic activity for the day determined from
     *  the sum of the 8 Ap indices. planetaryDailyCharacterFigure ranges, in steps of one-tenth,
     *  from 0 (quiet) to 2.5 (highly disturbed).
     */
    double planetaryDailyCharacterFigure;

    //! planetaryDailyCharacterFigureConverted.
    /*! A conversion of the 0-to-2.5 range of the planetaryDailyCharacterFigure index to one digit
     *  between 0 and 9.
     */
    unsigned int planetaryDailyCharacterFigureConverted;

    //! International sunspot number.
    /*! International sunspot number. Records contain the Zurich number through 1980 Dec 31 and
     * the International Brussels
     * number thereafter.
     */
    unsigned int internationalSunspotNumber;

    //! 10.7-cm Solar Radio Flux (F10.7) Adjusted to 1 AU.
    /*! 10.7-cm Solar Radio Flux (F10.7) measured at Ottawa at 1700 UT daily from 1947 Feb 14
     * until 1991 May 31 and measured at Penticton at 2000 UT from 1991 Jun 01 on.
     * Expressed in units of 10-22 W/m2/Hz.
     */
    double solarRadioFlux107Adjusted;

    //! Flux Qualifier.
    /*! Flux Qualifier with one of following values:
     * 0 indicates flux required no adjustment;
     * 1 indicates flux required adjustment for burst in progress at time of measurement;
     * 2 indicates a flux approximated by either interpolation or extrapolation;
     * 3 indicates no observation; and
     * 4 indicates CSSI interpolation of missing data.
     */
    unsigned int fluxQualifier;

    //! Centered 81-day arithmetic average of F10.7 (adjusted).
    /*! Centered 81-day arithmetic average of F10.7, adjusted to 1 AU.
     */
    double centered81DaySolarRadioFlux107Adjusted;

    //! Last 81-day arithmetic average of F10.7 (adjusted).
    /*! Last 81-day arithmetic average of F10.7, adjusted to 1 AU.
     */
    double last81DaySolarRadioFlux107Adjusted;

    //! Observed (unadjusted) value of F10.7.
    double solarRadioFlux107Observed;

    //! Centered 81-day arithmetic average of F10.7 (observed).
    double  centered81DaySolarRadioFlux107Observed;

    //! Last 81-day arithmetic average of F10.7 (observed).
    double  last81DaySolarRadioFlux107Observed;

    //! Vector containing 3-hourly planetary range indices (Kp).
    /*! Vector of eight 3-hourly planetary range indices.
     */
    Eigen::VectorXd planetaryRangeIndexVector;

    //! Vector containing 3-hourly planetary equivalent amplitude indices (Ap).
    /*! Vector of eight 3-hourly planetary equivalent amplitude indices.
     */
    Eigen::VectorXd planetaryEquivalentAmplitudeVector;

    //! String defining line's data-type
    /*! 1: "Observed"
     *  2: "Daily Predicted"
     *  3: "Monthly Predicted"
     *  4: "Monthly Fit"
     */
    unsigned int dataType;

    //! Overloaded ostream to print class information.
    /*! Overloaded ostream to print class information; prints all converted Solar Activity
     * variables listed in the http://celestrak.com/SpaceData/sw19571001.txt file.
     */
    friend std::ostream& operator << ( std::ostream& stream,
                                     SolarActivityData& solarActivityData );

protected:

private:
};

//! Pointer to a SolarActivityData structure
typedef std::shared_ptr< SolarActivityData > SolarActivityDataPtr;

//! Data map of SolarActivityData structure Pointers
typedef std::map< double , SolarActivityDataPtr >  SolarActivityDataMap ;

struct SolarActivityContainer
{
    SolarActivityContainer(
            const std::map< double, SolarActivityDataPtr >& solarActivityDataMap ):
        solarActivityDataMap_( solarActivityDataMap )
    {

        lookUpScheme_ = std::make_shared< interpolators::BinarySearchLookupScheme< double > >(
                    utilities::createVectorFromMapKeys( solarActivityDataMap ) );
    }

    std::shared_ptr< SolarActivityData > getSolarActivityData( const double time )
    {
        double julianDay = basic_astrodynamics::convertSecondsSinceEpochToJulianDay( time );
        return getSolarActivityDataAtJulianDay( julianDay );
    }

    std::shared_ptr< SolarActivityData > getSolarActivityDataAtJulianDay( const double julianDay )
    {
        int nearestJulianDay = lookUpScheme_->getIndependentVariableValue(
                    lookUpScheme_->findNearestLowerNeighbour( julianDay ) );
        return solarActivityDataMap_.at( nearestJulianDay );
    }


private:

     std::map< double, SolarActivityDataPtr > solarActivityDataMap_;

     std::shared_ptr< interpolators::LookUpScheme< double > > lookUpScheme_;

};

//! Function that reads a SpaceWeather data file
/*!
 * This function reads a SpaceWeather data file and returns a map with SolarActivityData
 *
 * \param filePath std::string
 * \return solarActivityDataMap std::map< double , SolarActivityDataPtr >
 */
SolarActivityDataMap readSolarActivityData( std::string filePath ) ;

} // namespace solar_activity
} // namespace input_output
} // namespace tudat

#endif // TUDAT_SOLAR_ACTIVITY_DATA_H
