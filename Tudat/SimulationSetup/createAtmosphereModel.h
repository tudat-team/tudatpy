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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_CREATEATMOSPHEREMODEL_H
#define TUDAT_CREATEATMOSPHEREMODEL_H

#include <string>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"

namespace tudat
{

namespace simulation_setup
{

//! List of atmosphere models available in simulations
/*!
 *  List of atmosphere models available in simulations. Atmosphere models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum AtmosphereTypes
{
    exponential_atmosphere,
    tabulated_atmosphere
};

//! Class for providing settings for atmosphere model.
/*!
 *  Class for providing settings for automatic atmosphere model creation. This class is a
 *  functional (base) class for settings of atmosphere models that require no information in
 *  addition to their type. Atmosphere model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class AtmosphereSettings
{
public:

    //! Constructor, sets type of atmosphere model.
    /*!
     *  Constructor, sets type of atmosphere model. Settings for atmosphere models requiring
     *  additional information should be defined in a derived class.
     *  \param atmosphereType Type of atmosphere model that is to be created.
     */
    AtmosphereSettings( const AtmosphereTypes atmosphereType ):
        atmosphereType_( atmosphereType ){ }

    //! Destructor
    virtual ~AtmosphereSettings( ){ }

    //! Function to return type of atmosphere model that is to be created.
    /*!
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( ){ return atmosphereType_; }

private:

    //!  Type of atmosphere model that is to be created.
    AtmosphereTypes atmosphereType_;
};

//! AtmosphereSettings for defining an exponential atmosphere.
class ExponentialAtmosphereSettings: public AtmosphereSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param densityScaleHeight Scale heigh for density profile of atmosphere.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at ground level.
     *  \param specificGasConstant Specific gas constant for (constant) atmospheric chemical
     *  composition.
     */

    ExponentialAtmosphereSettings(
            const double densityScaleHeight, const double constantTemperature,
            const double densityAtZeroAltitude, const double specificGasConstant ):
        AtmosphereSettings( exponential_atmosphere ),
        densityScaleHeight_( densityScaleHeight ), constantTemperature_( constantTemperature ),
        densityAtZeroAltitude_( densityAtZeroAltitude ), specificGasConstant_( specificGasConstant )
    { }

    //! Function to return scale heigh for density profile of atmosphere.
    /*!
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( ){ return densityScaleHeight_; }

    //! Function to return constant atmospheric temperature.
    /*!
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( ){ return constantTemperature_; }

    //! Function to return atmospheric density at ground level.
    /*!
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( ){ return densityAtZeroAltitude_; }

    //! Function to return specific gas constant for (constant) atmospheric chemical
    /*!
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( ){ return specificGasConstant_; }
private:

    //! Scale heigh for density profile of atmosphere.
    double densityScaleHeight_;

    //! Constant atmospheric temperature.
    double constantTemperature_;

    //! Atmospheric density at ground level.
    double densityAtZeroAltitude_;

    //! Specific gas constant for (constant) atmospheric chemical
    double specificGasConstant_;
};

//! AtmosphereSettings for defining an atmosphere with tabulated data from file.
class TabulatedAtmosphereSettings: public AtmosphereSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param atmosphereFile File containing atmospheric properties, file should contain
     *  four columns of atmospheric data with altitude, density, pressure and temperature,
     *  respecrively.
     */
    TabulatedAtmosphereSettings( const std::string& atmosphereFile ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereFile ){ }

    //! Function to return file containing atmospheric properties.
    /*!
     *  Function to return file containing atmospheric properties.
     *  \return File containing atmospheric properties.
     */
    std::string getAtmosphereFile( ){ return atmosphereFile_; }

private:

    //! File containing atmospheric properties.
    /*!
     *  File containing atmospheric properties, file should contain
     *  four columns of atmospheric data with altitude, density, pressure and temperature,
     *  respecrively.
     */
    std::string atmosphereFile_;
};

//! Function to create an atmosphere model.
/*!
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
boost::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel(
        const boost::shared_ptr< AtmosphereSettings > atmosphereSettings,
        const std::string& body );
}

}


#endif // TUDAT_CREATEATMOSPHEREMODEL_H
