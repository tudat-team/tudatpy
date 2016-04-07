/*!   Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined function.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      151104    J. Geul           Complete rewrite
 *
 *    References
 *
 */

#ifndef TUDAT_NRLMSISE00_ATMOSPHERE_H
#define TUDAT_NRLMSISE00_ATMOSPHERE_H

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

extern "C" {
    #include "nrlmsise-00.h"
}


namespace tudat {
namespace aerodynamics {

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

//! NRLMSISE-00 atmosphere model class.
/*!
 *  NRLMSISE-00 atmosphere model class. This class uses the NRLMSISE00 atmosphere model to calculate atmospheric
 *  density and temperature. The GTD7 function is used: Neutral Atmosphere Empircial Model from the surface to lower
 *  exosphere.
 *  Compositional data is currently not used (computation of pressure; speed of sound, etc.).
 *
 */
class NRLMSISE00Atmosphere : public AtmosphereModel {
 public:

    //! Solar activity function
    /*!
     * Boost function that accepts (altitude, longitude, latitude, time) and returns solar activity data.
     */
    typedef boost::function< NRLMSISE00Input( double, double, double, double ) >
        NRLMSISE00InputFunction;

    //! Default constructor.
    /*!
     * Default constructor.
     * \param nrlmsise00InputFunction Function which provides the NRLMSISE00 model input as a function of
     * (altitude, longitude, latitude, time).
     * \param useIdealGasLaw Variable denoting whether to use the ideal gas law for computation of pressure.
     */
    NRLMSISE00Atmosphere( const NRLMSISE00InputFunction nrlmsise00InputFunction,
                          const bool useIdealGasLaw = false )
        :nrlmsise00InputFunction_( nrlmsise00InputFunction ), useIdealGasLaw_( useIdealGasLaw )
    {
        if( useIdealGasLaw )
        {
            std::cerr<<"Warning, using ideal gas law with R at standard atmospheric conditions in pressure computation of NRLMSISE00InputFunction."<<std::endl;
        }
        resetHashKey( );
    }

    //! Get local density.
    /*!
     * Returns the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude [km].
     * \param longitude Longitude [deg].
     * \param latitude Latitude [deg].
     * \param time time since simulation start epoch.
     * \return Atmospheric density [kg/m^3].
     */
    double getDensity( const double altitude, const double longitude,
                       const double latitude, const double time)
    {
        computeProperties( altitude, longitude, latitude, time );
        return density_;
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * \param altitude Altitude.
     * \param longitude Longitude. (optional)
     * \param latitude Latitude. (optional)
     * \param time Time. (optional)
     * \return Atmospheric pressure.
     */
    double getPressure( const double altitude, const double longitude,
                        const double latitude, const double time )
    {
        if( useIdealGasLaw_ )
        {
            computeProperties( altitude, longitude, latitude, time );
        }
        else
        {
            throw std::runtime_error( "Error, non-ideal gas-law pressure-computation not yet implemented in NRLMSISE00Atmosphere." );
        }
        return pressure_;
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
    double getTemperature( const double altitude, const double longitude,
                           const double latitude, const double time )
    {
        computeProperties( altitude, longitude, latitude, time );
        return temperature_;
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
    double getSpeedOfSound( const double altitude, const double longitude,
                            const double latitude, const double time )
    {
        if( useIdealGasLaw_ )
        {
            computeProperties( altitude, longitude, latitude, time );
        }
        else
        {
            throw std::runtime_error( "Error, non-ideal gas-law pressure-computation not yet implemented in NRLMSISE00Atmosphere." );
        }

        return aerodynamics::computeSpeedOfSound( temperature_, 1.4, physical_constants::SPECIFIC_GAS_CONSTANT_AIR );
    }

    //! Get the full model output
    /*!
     * Gets the output directly from the model. This will return a
     * pair of double vectors containing density and temperature
     * values.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return Full density and temperature values
     */
    std::pair< std::vector< double >, std::vector< double > > getFullOutput(
        const double altitude, const double longitude,
        const double latitude, const double time);

    //! Reset the hash key
    /*!
     * Resets the hash key, this allows re-computation even if the
     * independent parameters haven't changed. Such as in the case of
     * changes to the model.
     */
    void resetHashKey( )
    {
        hashKey_ = 0;
    }

 private:
    //! Shared pointer to solar activity function
    NRLMSISE00InputFunction nrlmsise00InputFunction_;

    //! Current key hash
    size_t hashKey_;

    //! Current local density
    double density_;

    //! Current local temperature
    double temperature_;

    //! Current local pressure
    double pressure_;

    //! Variable denoting whether to use the ideal gas law for computation of pressure.
    bool useIdealGasLaw_;

    //! The file name of the solar activity data.
    struct nrlmsise_flags flags_;

    //!  Magnetic index structure of 7d array.
    /*!
     *   Magnetic index structure of 7d array:
     *   0 : daily AP
     *   1 : 3 hr AP index for current time
     *   2 : 3 hr AP index for 3 hrs before current time
     *   3 : 3 hr AP index for 6 hrs before current time
     *   4 : 3 hr AP index for 9 hrs before current time
     *   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
     *           prior to current time
     *   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
     *           prior to current time 
     */
    struct ap_array aph_;

    //! Input structure with 10 scalar settings and magnetic values struct.
    /*!
     *  Input structure with 10 scalar settings and magnetic values struct.
     *      year   - year, currently ignored           (int)
     *      doy    - day of the year                   (int)
     *      sec    - seconds in the day (UT)           (double)
     *      alt    - altitude in kilometers            (double)
     *      g_lat  - geodetic latitude                 (double)
     *      g_long - geodetic longitude                (double)
     *      lst    - local apparent solar time (hours) (double)
     *      f107A  - 81 day average of F10.7 flux (centered on doy) (double)
     *      f107   - daily F10.7 flux for previous day (double)
     *      ap     - magnetic index (daily)            (double)
     *      ap_a   - magnetic index struct (see above) (ap_array)       
     */
    struct nrlmsise_input input_;

    //! Ouput structure of densities (9d) and temperature (2d) arrays.
    /*!
     *  Ouput structure of densities (9d) and temperature (2d) arrays:
     *      d[0] - HE NUMBER DENSITY(CM-3)
     *      d[1] - O NUMBER DENSITY(CM-3)
     *      d[2] - N2 NUMBER DENSITY(CM-3)
     *      d[3] - O2 NUMBER DENSITY(CM-3)
     *      d[4] - AR NUMBER DENSITY(CM-3)                       
     *      d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
     *      d[6] - H NUMBER DENSITY(CM-3)
     *      d[7] - N NUMBER DENSITY(CM-3)
     *      d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
     *      t[0] - EXOSPHERIC TEMPERATURE
     *      t[1] - TEMPERATURE AT ALT
     * 
     *
     *      O, H, and N are set to zero below 72.5 km
     *
     *      t[0], Exospheric temperature, is set to global average for
     *      altitudes below 120 km. The 120 km gradient is left at global
     *      average value for altitudes below 72 km.
     *
     *      d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
     *      and GTD7D
     *
     *        SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
     *        species labeled by indices 0-4 and 6-7 in output variable d.
     *        This includes He, O, N2, O2, Ar, H, and N but does NOT include
     *        anomalous oxygen (species index 8).
     *
     *        SUBROUTINE GTD7D -- d[5] is the "effective total mass density
     *        for drag" and is the sum of the mass densities of all species
     *        in this model, INCLUDING anomalous oxygen.
     */
    struct nrlmsise_output output_;

    //! Get Hash Key.
    /*!
     * Returns hash key value based on a vector of keys
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     * \return hash key value.
     */
    size_t hashFunc( const double altitude, const double longitude,
                     const double latitude, const double time )
    {
        size_t seed = 0;
        boost::hash_combine( seed, boost::hash< double >( )( altitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( longitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( latitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( time ) );
        return seed;
    }

    //! Compute the local atmospheric properties.
    /*!
     * Computes the local atmospheric density, pressure and temperature.
     * \param altitude Altitude.
     * \param longitude Longitude.
     * \param latitude Latitude.
     * \param time Time.
     */
    void computeProperties( const double altitude, const double longitude,
                            const double latitude, const double time );
};

}  // namespace aerodynamics
}  // namespace tudat

#endif // TUDAT_NRLMSISE00_ATMOSPHERE_H_
