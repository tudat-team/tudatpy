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

#ifndef TUDAT_MARS_DTM_ATMOSPHERE_MODEL_H
#define TUDAT_MARS_DTM_ATMOSPHERE_MODEL_H

#include <iostream>
#include <vector>
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/dateTime.h"

namespace tudat
{

namespace aerodynamics
{
    class marsDate {
    public:
        int year;
        int month;
        int day;
        double hours;
        double minutes;
        double seconds;

        //constructor
        //structure of the marsDate to be defined later
        marsDate(int y, int m, int d, double hh, double mm, double ss) : year(y), month(m), day(d), hours(hh), minutes(mm), seconds(ss) {}

        //distructor
        virtual ~marsDate( ) { }

        // takes the absolute difference between the two dates after converting them to Julian date.
        double dateDifference(const marsDate& date1, const marsDate& date2) {
            return std::abs( basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    date1.year, date1.month, date1.day, date1.hours, date1.minutes, date1.seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000) -
                             basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                            date2.year, date2.month, date2.day, date2.hours, date2.minutes, date2.seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000));
        }

        // function to find the nearest date in the list of dates provided in the next variable
        marsDate findNearestDate(const marsDate &targetDate);

        double marsDayofYear(const marsDate &targetDate);
    };

    //from: https://www.planetary.org/articles/mars-calendar
    static std::vector<marsDate> marsYears = {
            marsDate(1955, 4, 11, 0, 0, 0),
            marsDate(1957, 2, 26, 0, 0, 0),
            marsDate(1959, 1, 14, 0, 0, 0),
            marsDate(1960, 12, 1, 0, 0, 0),
            marsDate(1962, 10, 19, 0, 0, 0),
            marsDate(1964, 9, 5, 0, 0, 0),
            marsDate(1966, 7, 24, 0, 0, 0),
            marsDate(1968, 6, 10, 0, 0, 0),
            marsDate(1970, 4, 28, 0, 0, 0),
            marsDate(1972, 3, 15, 0, 0, 0),
            marsDate(1974, 1, 31, 0, 0, 0),
            marsDate(1975, 12, 19, 0, 0, 0),
            marsDate(1977, 11, 5 ,0, 0, 0),
            marsDate(1979, 9, 23, 0, 0, 0),
            marsDate(1981, 8, 10, 0, 0, 0),
            marsDate(1983, 6, 28, 0, 0, 0),
            marsDate(1985, 5, 15, 0, 0, 0),
            marsDate(1987, 4, 1, 0, 0, 0),
            marsDate(1989, 2, 16, 0, 0, 0),
            marsDate(1991, 1, 4, 0, 0, 0),
            marsDate(1992, 11, 21 , 0, 0, 0),
            marsDate(1994, 10, 9 , 0, 0, 0),
            marsDate(1996, 8, 26 , 0, 0, 0),
            marsDate(1998, 7, 14 , 0, 0, 0),
            marsDate(2000, 5, 31 , 0, 0, 0),
            marsDate(2002, 4, 18 , 0, 0, 0),
            marsDate(2004, 3, 5 , 0, 0, 0),
            marsDate(2006, 1, 21 , 0, 0, 0),
            marsDate(2007, 12, 9 , 0, 0, 0),
            marsDate(2009, 10, 26 , 0, 0, 0),
            marsDate(2011, 9, 13 , 0, 0, 0),
            marsDate(2013, 7, 31 , 0, 0, 0),
            marsDate(2015, 6, 18 , 0, 0, 0),
            marsDate(2017, 5, 5 , 0, 0, 0),
            marsDate(2019, 3, 23 , 0, 0, 0),
            marsDate(2021, 2, 7 , 0, 0, 0),
            marsDate(2022, 12, 26 , 0, 0, 0),
            marsDate(2024, 11, 12 , 0, 0, 0),
            marsDate(2026, 9, 30 , 0, 0, 0),
            marsDate(2028, 8, 17 , 0, 0, 0)
    };

// Function to compute the solar longitude of Mars
std::tuple<double, double> computeSolarLongitude(const double longitude, const int day ,const int month,const int year, const int hours, const int minutes);

class MarsDtmAtmosphereModel: public AtmosphereModel
{
public:


    MarsDtmAtmosphereModel(
        const std::function< double( const double ) > f107Function = [](const double){return 65.0;} );

    //! Default destructor.
    /*!
    * Default destructor.
    */
    virtual ~MarsDtmAtmosphereModel( ) { }


    //! Get local density.
    /*!
    * Returns the local density parameter of the atmosphere in kg per meter^3.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric density.
    */
    virtual double getDensity( const double altitude, const double longitude,
                               const double latitude, const double time )
    {
        computeProperties( altitude, longitude, latitude, time );
        return currentDensity_;
    }

    //! Get local pressure.
    /*!
    * Returns the local pressure of the atmosphere parameter in Newton per meter^2.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric pressure.
    */
    virtual double getPressure( const double altitude, const double longitude,
                                const double latitude, const double time )
    {
        throw std::runtime_error( "Error wen getting Mars DTM pressure; model not yet implemented." );
        return TUDAT_NAN;
    }

    //! Get local temperature.
    /*!
    * Returns the local temperature of the atmosphere parameter in Kelvin.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric temperature.
    */
    virtual double getTemperature( const double altitude, const double longitude,
                                   const double latitude, const double time )
    {
        computeProperties( altitude, longitude, latitude, time );
        return currentTemperature_z;
    }

    //! Get local speed of sound.
    /*!
    * Returns the local speed of sound of the atmosphere in m/s.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Atmospheric speed of sound.
    */
    virtual double getSpeedOfSound( const double altitude, const double longitude,
                                    const double latitude, const double time )
    {
        throw std::runtime_error( "Error wen getting Mars DTM speed of sounds; model not yet implemented." );
        return TUDAT_NAN;
    }



    const double Omega = 2.0*mathematical_constants::PI/686.98;//(686.98*24.63);//rad/h//2*M_PI/(659355072); //rad/s
    const double omega = 2.0*mathematical_constants::PI/88668;//24.63;//2*M_PI/(88668); // rad/s

    double computeLocalSolarTime(const double longitude, const int day, const int month, const int year, const int hours,
                                 const int minutes);
    void updateLegendrePolynomials( const double latitude );
    double computeGl(const double latitude, const double longitude, const int minutes, const int hours_, const int day_,
                     const int month_, const int year_, const int indexg);
    double computeGeopotentialAltitude( const double altitude );

    double computeCurrentTemperature( const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexg);

    double computeGamma(const double latitude, const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_, const int indexm);

    double heightDistributionFunction(const double altitude, const double latitude,
                                      const double longitude, const int minutes_, const int hours_, const int day_ ,const int month_,const int year_,const int indexm);

    double getTotalDensity( const double altitude, const double latitude,
                            const double longitude, const int minutes_, const int hours_, const int day_ , const int month_, const int year_);

    double getCurrentLegendrePolynomial( const int degree )
    {
       return currentLegendrePolynomials_.at( degree );
    }

    double computeGl_Subr(const double latitude, const double longitude, const int minutes_, const int hours_,
                          const int day_, const int month_, const int year_, const int indexg);

    void defineDustPars(const int rows);

    double computeDustStorm(const double Ls);

    double getSolarFluxIndex() const{
        return currentF107_;
    }
    double getSolarLongitude() const{
        std::cout << "Solar Longitude: " << Ls_ << std::endl;
        return Ls_;
    }
protected:

private:

    void computeProperties(  const double altitude, const double longitude,
                             const double latitude, const double time );

    double polarRadius_;

    std::basic_string<char, std::char_traits<char>, std::allocator<char>> filename_;

    std::function< double( const double ) > f107Function_;

    std::vector< double > alpha_;
    std::vector< double > mmass_;

    std::vector<std::vector<double>> coefficients_;
    std::vector<double> taus;

    double Ls_;
    double currentF107_;
    std::vector< double > currentLegendrePolynomials_;
    double currentGeopotentialAltitude_;
    double currentTemperature_138;
    double currentTemperature_inf;
    double currentdTemperature_;
    double currentDensity_;
    double sigma;
    double currentTemperature_z;

};


} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_MARS_DTM_ATMOSPHERE_MODEL_H
