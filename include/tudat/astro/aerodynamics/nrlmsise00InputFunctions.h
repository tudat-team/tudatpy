/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H
#define TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H

#include <vector>
#include <cmath>

#include "tudat/io/solarActivityData.h"


namespace tudat
{

namespace aerodynamics
{


//! Struct for Solar Activity data.
/*!
 *
 */
struct NRLMSISE00Input
{
    //! Input for computation of NRLMSISE00 atmospheric conditions at current time and position.
    /*!
     * Input for computation of NRLMSISE00 atmospheric conditions at current time and position. The computation of
     * this struct may be reperformed every time step, to reflect the changes in atmospheric condition.
     * \param year Current year
     * \param dayOfTheYear Day in the current year
     * \param secondOfTheDay Number of seconds into the current day.
     * \param localSolarTime Local solar time at the computation position
     * \param f107 Current daily F10.7 flux for previous day
     * \param f107a 81 day average of F10.7 flux (centered on current dayOfTheYear).
     * \param apDaily Current daily magnetic index
     * \param apVector Current magnetic index data vector: \sa ap_array
     * \param switches List of NRLMSISE-specific flags: \sa nrlmsise_flags
     */
    NRLMSISE00Input(
            int year = 0, int dayOfTheYear = 0, double secondOfTheDay = 0.0,
            double localSolarTime = 0.0, double f107 = 0.0, double f107a = 0.0,
            double apDaily = 0.0,
            std::vector< double > apVector = std::vector< double>( 7, 0.0 ),
            std::vector< int > switches = std::vector< int >( ) )
        : year( year ), dayOfTheYear( dayOfTheYear ), secondOfTheDay( secondOfTheDay) ,
          localSolarTime( localSolarTime ), f107( f107 ), f107a( f107a ),
          apDaily( apDaily ), apVector( apVector ), switches( switches )
    {
        if( switches.empty( ) )
        {
            this->switches = std::vector< int >( 24, 1 );
            this->switches[ 0 ] = 0;
        }
    }

    //! Current year
    int year;

    //! Day in the current year
    int dayOfTheYear;

    //! Number of seconds into the current day.
    double secondOfTheDay;\

    //! Local solar time at the computation position
    double localSolarTime;

    //! Current daily F10.7 flux for previous day
    double f107;

    //! 81 day average of F10.7 flux (centered on current dayOfTheYear).
    double f107a;

    //! Current daily magnetic index
    double apDaily;

    //! Current magnetic index data vector: \sa ap_array
    std::vector< double > apVector;

    //! List of NRLMSISE-specific flags: \sa nrlmsise_flag
    std::vector< int > switches;
};


//! NRLMSISE00 Input function
/*!
 * This function is used to define the input for the NRLMSISE model.
 * This function reads solar activity data and defines the input using this data.
 * \param altitude Altitude at which output is to be computed [m].
 * \param longitude Longitude at which output is to be computed [rad].
 * \param latitude Latitude at which output is to be computed [rad].
 * \param time Time at which output is to be computed (seconds since J2000).
 * \param solarActivityMap SolarActivityData structure
 * \param adjustSolarTime Boolean denoting whether the computed local solar time should be overidden with localSolarTime
 * input.
 * \param localSolarTime Local solar time that is used when adjustSolarTime is set to true.
 * \return NRLMSISE00Input nrlmsiseInputFunction
 */
NRLMSISE00Input nrlmsiseInputFunction( const double altitude, const double longitude,
                                       const double latitude, const double time,
                                       const tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap,
                                       const bool adjustSolarTime = false, const double localSolarTime = 0.0 );

}  // namespace aerodynamics
}  // namespace tudat

#endif // TUDAT_NRLMSISE00_ATMOSPHERE_H_
