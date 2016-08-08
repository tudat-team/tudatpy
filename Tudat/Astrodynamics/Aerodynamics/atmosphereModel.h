/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Set time to zero by default.
 *      110412    F.M. Engelen      Included flag boolean isAltitudeModelDefault_.
 *      110427    F.M. Engelen      Changed input arguments to altitude, longitude, latitude.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      130120    K. Kumar          Made function calls const-correct; added shared-pointer
 *                                  typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_ATMOSPHERE_MODEL_H
#define TUDAT_ATMOSPHERE_MODEL_H

#include <boost/shared_ptr.hpp>
#include "tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace aerodynamics
{

//! Gas component properties data structure
/*!
 * This data structure contains the molar mass and collision diameter of
 * the most frequently observed gasses in the atmosphere.
 */
struct GasComponentProperties
{
    //! default constructor
    /*!
     * Constructs a GasComponentProperties data structure with values obtained from the following sources:
     * Collision diameters from: Tables of Physical & Chemical Constants Kaye & Laby Online, Kaye & Laby Online, 2016
     * (http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_4.html)
     * Molar mass from: NIST,2016 (http://www.nist.gov/pml/data/images/illo_for_2014_PT_1.PNG)
     *
     * Data to be verified
     *
     */
    GasComponentProperties():
    diameterArgon(340E-12), diameterAtomicHydrogen(260E-12), diameterHelium(256E-12),
      diameterNitrogen(370E-12), diameterOxygen(358E-12), diameterAtomicNitrogen(290E-12),
      diameterAtomicOxygen(280E-12), molarMassArgon(39.948E-3), molarMassAtomicHydrogen(1.008E-3),
      molarMassHelium(4.002602E-3), molarMassNitrogen(2.0*14.007E-3), molarMassOxygen(2.0*15.999E-3),
      molarMassAtomicNitrogen(14.007E-3), molarMassAtomicOxygen(15.999E-3)
    { }

    //! Molecular colision diameter of Argon in m
    double diameterArgon;

    //! Molecular colision diameter of Atomic Hydrogen in m
    double diameterAtomicHydrogen;

    //! Molecular colision diameter of Helium in m
    double diameterHelium;

    //! Molecular colision diameter of Nitrogen in m
    double diameterNitrogen;

    //! Molecular colision diameter of Oxygen in m
    double diameterOxygen;

    //! Molecular colision diameter of Atomic Nitrogen in m
    double diameterAtomicNitrogen;

    //! Molecular colision diameter of Atomic Oxygen in m
    double diameterAtomicOxygen;


    //! molar mass of Argon in kg/mole
    double molarMassArgon;

    //! Molar mass of Atomic Hydrogen in kg/mole
    double molarMassAtomicHydrogen;

    //! Molar mass of Helium in kg/mole
    double molarMassHelium;

    //! Molar mass of Nitrogen in kg/mole
    double molarMassNitrogen;

    //! Molar mass of Oxygen in kg/mole
    double molarMassOxygen;

    //! Molar mass of Atomic Nitrogen in kg/mole
    double molarMassAtomicNitrogen;

    //! Molar mass of Atomic Oxygen in kg/mole
    double molarMassAtomicOxygen;
};

//! Atmosphere model class.
/*!
 * Base class for all atmosphere models.
 * To keep the function generic for both reference atmospheres and standard
 * atmospheres, altitude, longitude, latitude, and time are inputs for all functions.
 */
class AtmosphereModel
{
public:

    //! Default destructor.
    /*!
    * Default destructor.
    */
    virtual ~AtmosphereModel( ) { }

    //! Set gas component properties.
    /*!
    * Sets the gas component properties, which contains the molar mass and collision diameter of the molecules.
    * \param GasComponentProperties .
    * \return void.
    */
    virtual void setGasComponentProperties( GasComponentProperties GasComponentProperties ){}

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
                               const double latitude, const double time ) = 0;

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
                                const double latitude, const double time ) = 0;

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
                                   const double latitude, const double time ) = 0;

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
                                    const double latitude, const double time ) = 0;

    //! Get mean free path.
    /*!
    * This function returns the mean free path, which can be calculated using the number density and the collision diameter.
    * The function is not implemented by default and returns NaN if the atmosphere model doesn't implement this function.
    * \param altitude Altitude.
    * \param longitude Longitude.
    * \param latitude Latitude.
    * \param time Time.
    * \return double mean free path.
    */
    virtual double getMeanFreePath(const double altitude, const double longitude,
                                   const double latitude, const double time)
    {
        return TUDAT_NAN;
    }

    //! Reset Hash key
    /*!
    * This function resets the hashKey, which is used to check if recalculation of the atmospheric properties is required or not.
    * \return void
    */
    virtual void resetHashKey(){}

protected:

private:
};

//! Typedef for shared-pointer to AtmosphereModel object.
typedef boost::shared_ptr< AtmosphereModel > AtmosphereModelPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_ATMOSPHERE_MODEL_H
