/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NRLMSISE00_ATMOSPHERE_H
#define TUDAT_NRLMSISE00_ATMOSPHERE_H

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

#include <functional>
#include <boost/functional/hash.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/io/solarActivityData.h"

extern "C" {
#include <nrlmsise00/nrlmsise-00.h>
}

namespace tudat
{

namespace aerodynamics
{

//! NRLMSISE-00 atmosphere model class.
/*!
 *  NRLMSISE-00 atmosphere model class. This class uses the NRLMSISE00 atmosphere model to calculate atmospheric
 *  density and temperature. The GTD7 function is used: Neutral Atmosphere Empircial Model from the surface to lower
 *  exosphere.
 *  Currently the ideal gas law is used to compute the speed of sound.
 *  The specific heat ratio is assumed to be constant and equal to 1.4.
 */
class NRLMSISE00Atmosphere : public AtmosphereModel
{
public:
    //! NRLMSISEInput function
    /*!
     * Boost function that accepts (altitude, longitude, latitude, time ) and returns NRLMSISEInput data.
     */
    typedef std::function< NRLMSISE00Input( double, double, double, double ) > NRLMSISE00InputFunction;

    //! Default constructor.
    /*!
     * Default constructor.
     * \param nrlmsise00InputFunction Function which provides the NRLMSISE00 model input as a function of
     * (altitude, longitude, latitude, time ).
     * \param useIdealGasLaw Variable denoting whether to use the ideal gas law for computation of pressure.
     */
    NRLMSISE00Atmosphere( const NRLMSISE00InputFunction nrlmsise00InputFunction,
                          const bool useIdealGasLaw = true,
                          const bool useAnomalousOxygen = true ):
        AtmosphereModel( true, true ), nrlmsise00InputFunction_( nrlmsise00InputFunction ), useIdealGasLaw_( useIdealGasLaw ),
        useAnomalousOxygen_( useAnomalousOxygen )
    {
        resetHashKey( );
        specificHeatRatio_ = 1.4;
    }

    NRLMSISE00Atmosphere( const tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData,
                          const bool useIdealGasLaw = true,
                          const bool useStormConditions = false,
                          const bool useAnomalousOxygen = true ):
        AtmosphereModel( true, true ), solarActivityContainer_( solarActivityData ), useIdealGasLaw_( useIdealGasLaw ),
        useAnomalousOxygen_( useAnomalousOxygen )
    {
        nrlmsise00InputFunction_ = std::bind( &tudat::aerodynamics::nrlmsiseInputFunction,
                                              std::placeholders::_1,
                                              std::placeholders::_2,
                                              std::placeholders::_3,
                                              std::placeholders::_4,
                                              solarActivityContainer_,
                                              false,
                                              TUDAT_NAN,
                                              useStormConditions ? -1 : 1 );

        resetHashKey( );
        specificHeatRatio_ = 1.4;
    }

    //! Get local density.
    /*!
     * Returns the local density of the atmosphere in kg per meter^3.
     * \param altitude Altitude at which density is to be computed [m].
     * \param longitude Longitude at which density is to be computed [rad].
     * \param latitude Latitude at which density is to be computed [rad].
     * \param time Time at which density is to be computed (seconds since J2000).
     * \return Atmospheric density [kg/m^3].
     */
    double getDensity( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return density_;
    }

    //! Get local pressure.
    /*!
     * Returns the local pressure of the atmosphere in Newton per meter^2.
     * The pressure is not implemented in the current version (returns NaN).
     * \param altitude Altitude at which pressure is to be computed [m].
     * \param longitude Longitude at which pressure is to be computed [rad].
     * \param geodeticLatitude Latitude at which pressure is to be computed [rad].
     * \param time Time at which pressure is to be computed (seconds since J2000).
     * \return Atmospheric pressure.
     */
    double getPressure( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        if( useIdealGasLaw_ )
        {
            computeProperties( altitude, longitude, geodeticLatitude, time );
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
     * \param altitude Altitude at which temperature is to be computed [m].
     * \param longitude Longitude at which temperature is to be computed [rad].
     * \param geodeticLatitude Latitude at which temperature is to be computed [rad].
     * \param time Time at which temperature is to be computed (seconds since J2000).
     * \return Atmospheric temperature.
     */
    double getTemperature( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return temperature_;
    }

    //! Get local speed of sound.
    /*!
     * Returns the local speed of sound in m/s.
     * \param altitude Altitude at which speed of sound is to be computed [m].
     * \param longitude Longitude at which speed of sound is to be computed [rad].
     * \param geodeticLatitude Latitude at whichspeed of sound is to be computed [rad].
     * \param time Time at which speed of sound is to be computed (seconds since J2000).
     * \return Speed of sound.
     */
    double getSpeedOfSound( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return speedOfSound_;
    }

    //! Get local mean free path.
    /*!
     * Returns the local mean free path in m.
     * \param altitude Altitude at which mean free path is to be computed [m].
     * \param longitude Longitude at which mean free path is to be computed [rad].
     * \param geodeticLatitude Latitude at which mean free path is to be computed [rad].
     * \param time Time at which mean free path is to be computed (seconds since J2000).
     * \return Mean free path.
     */
    double getMeanFreePath( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return meanFreePath_;
    }

    //! Get local mean molar mass.
    /*!
     * Returns the local mean molar mass in kg/mol.
     * \param altitude Altitude at which mean molar mass is to be computed [m].
     * \param longitude Longitude at which mean molar mass  is to be computed [rad].
     * \param geodeticLatitude Latitude at which mean molar mass  is to be computed [rad].
     * \param time Time at which mean molar mass  is to be computed (seconds since J2000).
     * \return mean molar mass.
     */
    double getMeanMolarMass( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return meanMolarMass_;
    }

    //! get local number density of the gas components.
    /*!
     * Returns the number density of each gas component in a vector.
     * \param altitude Altitude at which number density is to be computed [m].
     * \param longitude Longitude at which number density is to be computed [rad].
     * \param geodeticLatitude Latitude at which number density is to be computed [rad].
     * \param time Time at which number density is to be computed (seconds since J2000).
     * \return Number densities of gas components
     */
    std::vector< double > getNumberDensities( const double altitude,
                                              const double longitude,
                                              const double geodeticLatitude,
                                              const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return numberDensities_;
    }

    virtual double getNumberDensity( const AtmosphericCompositionSpecies species,
                                     const double altitude,
                                     const double longitude,
                                     const double geodeticLatitude,
                                     const double time )
    {
        if( !( speciesIndices.count( species ) > 0 ) )
        {
            throw std::runtime_error( "Error, NRLMSISE00 has no dependency on species " + std::to_string( species ) + "." );
        }
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return numberDensities_.at( speciesIndices.at( species ) );
    }

    bool doesModelDefineSpeciesNumberDensity( const AtmosphericCompositionSpecies species )
    {
        if( speciesIndices.count( species ) > 0 )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Get local average number density.
    /*!
     * Returns the local average number density in m^-3
     * \param altitude Altitude at which density is to be computed [m].
     * \param longitude Longitude at which density is to be computed [rad].
     * \param geodeticLatitude Latitude at which density is to be computed [rad].
     * \param time Time at which density is to be computed (seconds since J2000).
     * \return average number density.
     */
    double getAverageNumberDensity( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return averageNumberDensity_;
    }

    //! Get local weighted average collision diameter.
    /*!
     * Returns the local weighted average collision diameter using the number densities as weights.
     * \param altitude Altitude at which average collision diameter is to be computed [m].
     * \param longitude Longitude at which average collision diameter is to be computed [rad].
     * \param geodeticLatitude Latitude at which average collision diameter is to be computed [rad].
     * \param time Time at which average collision diameter is to be computed (seconds since J2000).
     * \return weighted average collision diameter.
     */
    double getWeightedAverageCollisionDiameter( const double altitude,
                                                const double longitude,
                                                const double geodeticLatitude,
                                                const double time )
    {
        computeProperties( altitude, longitude, geodeticLatitude, time );
        return weightedAverageCollisionDiameter_;
    }

    //! Get the full model output
    /*!
     * Gets the output directly from the model. This will return a
     * pair of double vectors containing density and temperature
     * values.
     * \param altitude Altitude at which output is to be computed [m].
     * \param longitude Longitude at which output is to be computed [rad].
     * \param geodeticLatitude Latitude at which output is to be computed [rad].
     * \param time Time at which output is to be computed (seconds since J2000).
     * \return Full density and temperature values
     */
    std::pair< std::vector< double >, std::vector< double > > getFullOutput( const double altitude,
                                                                             const double longitude,
                                                                             const double geodeticLatitude,
                                                                             const double time );

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

    input_output::solar_activity::SolarActivityContainer& getSolarActivityContainer( )
    {
        return solarActivityContainer_;
    }

    //! Function to get  Input data to NRLMSISE00 atmosphere model
    /*!
     *  Function to get input data to NRLMSISE00 atmosphere model
     *  \return Input data to NRLMSISE00 atmosphere model
     */
    NRLMSISE00Input getNRLMSISE00Input( )
    {
        return inputData_;
    }

    NRLMSISE00InputFunction getNrlmsise00InputFunction( )
    {
        return nrlmsise00InputFunction_;
    }

    nrlmsise_input getNRLMSISE00InputStruct( )
    {
        return input_;
    }

    void setInputStruct( const double altitude, const double longitude, const double geodeticLatitude, const double time );

private:
    //! Shared pointer to solar activity function
    NRLMSISE00InputFunction nrlmsise00InputFunction_;

    //! Current key hash
    size_t hashKey_;

    //!  Current local density (kg/m3)
    double density_;

    //! Current local temperature (K)
    double temperature_;

    //! Current local pressure (Implemented with ideal gass law only!)
    double pressure_;

    //! Current speed of sound (m/s)
    double speedOfSound_;

    //! Current mean free path (m)
    double meanFreePath_;

    /*!
     *  Current number densities of gas components
     *      numberDensities_[0] - HE NUMBER DENSITY     (M-3)
     *      numberDensities_[1] - O NUMBER DENSITY      (M-3)
     *      numberDensities_[2] - N2 NUMBER DENSITY     (M-3)
     *      numberDensities_[3] - O2 NUMBER DENSITY     (M-3)
     *      numberDensities_[4] - AR NUMBER DENSITY     (M-3)
     *      numberDensities_[5] - H NUMBER DENSITY      (M-3)
     *      numberDensities_[6] - N NUMBER DENSITY      (M-3)
     *      numberDensities_[7] - Anomalous oxygen NUMBER DENSITY   (M-3)
     */
    std::vector< double > numberDensities_;

    //! Current average number density (M-3)
    double averageNumberDensity_;

    //! Current weighted average of the collision diameter using the number density as weights in (M)
    double weightedAverageCollisionDiameter_;

    //! mean molar mass (kg/mole)
    double meanMolarMass_;

    //! Specific heat ratio
    double specificHeatRatio_;

    //! Flags set for NRLMSISE computations
    nrlmsise_flags flags_;

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
    ap_array aph_;

    //! Input structure with 10 scalar settings and magnetic values struct.
    /*!
     *  Input structure with 10 scalar settings and magnetic values struct.
     *      year   - year, currently ignored           (int)
     *      doy    - day of the year                   (int)
     *      sec    - seconds in the day (UT)           (double)
     *      alt    - altitude in kilometers            (double)
     *      g_lat  - geodetic latitude in deg          (double)
     *      g_long - geodetic longitude in deg         (double)
     *      lst    - local apparent solar time (hours) (double)
     *      f107A  - 81 day average of F10.7 flux (centered on doy) (double)
     *      f107   - daily F10.7 flux for previous day (double)
     *      ap     - magnetic index (daily)            (double)
     *      ap_a   - magnetic index struct (see above) (ap_array)
     */
    nrlmsise_input input_;

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
    nrlmsise_output output_;

    //! Get Hash Key.
    /*!
     * Returns hash key value based on a vector of keys
     * \param altitude Altitude [m].
     * \param longitude Longitude [rad].
     * \param geodeticLatitude Latitude [rad].
     * \param time time.
     * \return hash key value.
     */
    size_t hashFunc( const double altitude, const double longitude, const double geodeticLatitude, const double time )
    {
        size_t seed = 0;
        boost::hash_combine( seed, boost::hash< double >( )( altitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( longitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( geodeticLatitude ) );
        boost::hash_combine( seed, boost::hash< double >( )( time ) );
        return seed;
    }

    //! Compute the local atmospheric properties.
    /*!
     * Computes the local atmospheric density, pressure and temperature.
     * \param altitude Altitude at which output is to be computed [m].
     * \param longitude Longitude at which output is to be computed [rad].
     * \param geodeticLatitude Latitude at which output is to be computed [rad].
     * \param time Time at which output is to be computed (seconds since J2000).
     */
    void computeProperties( const double altitude, const double longitude, const double geodeticLatitude, const double time );

    //! Input data to NRLMSISE00 atmosphere model
    NRLMSISE00Input inputData_;

    input_output::solar_activity::SolarActivityContainer solarActivityContainer_;

    //! Use the ideal gas law for the computation of the pressure.
    bool useIdealGasLaw_;

    bool useAnomalousOxygen_;

    std::map< AtmosphericCompositionSpecies, int > speciesIndices = { { he_species, 0 }, { o_species, 1 },          { n2_species, 2 },
                                                                      { o2_species, 3 }, { ar_species, 4 },         { h_species, 5 },
                                                                      { n_species, 6 },  { anomalous_o_species, 7 } };

};

Eigen::VectorXd getNrlmsiseInputAsVector( const nrlmsise_input& input );

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_NRLMSISE00_ATMOSPHERE_H_
