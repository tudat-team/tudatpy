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

#ifndef TUDAT_CREATEGRAVITYFIELD_H
#define TUDAT_CREATEGRAVITYFIELD_H

#include <iostream>

#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! List of gravity field models available in simulations
/*!
 *  List of gravity field models available in simulations. Gravity field models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum GravityFieldType
{
    central,
    central_spice,
    spherical_harmonic
};

//! Class for providing settings for gravity field model.
/*!
 *  Class for providing settings for automatic gravity field model creation. This class is a
 *  functional (base) class for settings of gravity field models that require no information in
 *  addition to their type. Gravity field model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class GravityFieldSettings
{
public:

    //! Constructor, sets type of gravity field model.
    /*!
     *  Constructor, sets type of gravity field model. Settings for gravity field models requiring
     *  additional information should be defined in a derived class.
     *  \param gravityFieldType Type of gravity field model that is to be created.
     */
    GravityFieldSettings( const GravityFieldType gravityFieldType ):
        gravityFieldType_( gravityFieldType ){ }

    //! Destructor
    virtual ~GravityFieldSettings( ){ }

    //! Function to return type of gravity field model that is to be created.
    /*!
     *  Function to return type of gravity field model that is to be created.
     *  \return Type of gravity field model that is to be created.
     */
    GravityFieldType getGravityFieldType( ){ return gravityFieldType_; }

protected:

    //! Type of gravity field model that is to be created.
    GravityFieldType gravityFieldType_;
};

//! Derived class of GravityFieldSettings defining settings of point mass gravity field.
class CentralGravityFieldSettings: public GravityFieldSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param gravitationalParameter Gravitational parameter of gravity field.
     */
    CentralGravityFieldSettings( double gravitationalParameter ):GravityFieldSettings( central ),
        gravitationalParameter_( gravitationalParameter ){ }

    //! Function to return gravitational parameter for gravity field.
    /*!
     *  Function to return gravitational parameter for gravity field.
     *  \return Gravitational parameter for gravity field.
     */
    double getGravitationalParameter( ){ return gravitationalParameter_; }
private:

    //! Gravitational parameter for gravity field.
    double gravitationalParameter_;
};


//! Derived class of GravityFieldSettings defining settings of spherical harmonic gravity
//! field representation.
class SphericalHarmonicsGravityFieldSettings: public GravityFieldSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param gravitationalParameter Gravitational parameter of gravity field.
     *  \param referenceRadius Reference radius of spherical harmonic field expansion.
     *  \param cosineCoefficients Cosine spherical harmonic coefficients (geodesy normalized).
     *  \param sineCoefficients Sine spherical harmonic coefficients (geodesy normalized).
     *  \param associatedReferenceFrame Identifier for body-fixed reference frame to which
     *  the coefficients are referred.
     */
    SphericalHarmonicsGravityFieldSettings( const double gravitationalParameter,
                                            const double referenceRadius,
                                            const Eigen::MatrixXd cosineCoefficients,
                                            const Eigen::MatrixXd sineCoefficients,
                                            const std::string associatedReferenceFrame ):
        GravityFieldSettings( spherical_harmonic ),
        gravitationalParameter_( gravitationalParameter ),
        referenceRadius_( referenceRadius ),
        cosineCoefficients_( cosineCoefficients ),
        sineCoefficients_( sineCoefficients ),
        associatedReferenceFrame_( associatedReferenceFrame )
    {  }

    //! Function to return gravitational parameter for gravity field.
    /*!
     *  Function to return gravitational parameter for gravity field.
     *  \return Gravitational parameter for gravity field.
     */
    double getGravitationalParameter( ){ return gravitationalParameter_; }

    //! Function to return reference radius of spherical harmonic field expansion
    /*!
     *  Function to return reference radius of spherical harmonic field expansion
     *  \return Reference radius of spherical harmonic field expansion
     */
    double getReferenceRadius( ){ return referenceRadius_; }

    //! Function to return cosine spherical harmonic coefficients (geodesy normalized).
    /*!
     *  Function to return cosine spherical harmonic coefficients (geodesy normalized).
     *  \return Cosine spherical harmonic coefficients (geodesy normalized).
     */
    Eigen::MatrixXd getCosineCoefficients( ){ return cosineCoefficients_; }

    //! Function to return sine spherical harmonic coefficients (geodesy normalized).
    /*!
     *  Function to return sine spherical harmonic coefficients (geodesy normalized).
     *  \return Sine spherical harmonic coefficients (geodesy normalized).
     */
    Eigen::MatrixXd getSineCoefficients( ){ return sineCoefficients_; }

    //! Function to return identifier for body-fixed reference frame.
    /*!
     *  Function to return identifier for body-fixed reference frame to which the coefficients
     *  are referred.
     *  \return Identifier for body-fixed reference frame to which the coefficients are referred.
     */
    std::string getAssociatedReferenceFrame( ){ return associatedReferenceFrame_; }

private:


    //! Gravitational parameter for gravity field that is to be created.
    double gravitationalParameter_;

    //! Reference radius of spherical harmonic field expansion
    double referenceRadius_;

    //! Cosine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd cosineCoefficients_;

    //! Sine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd sineCoefficients_;

    //! Identifier for body-fixed reference frame to which the coefficients are referred.
    std::string associatedReferenceFrame_;
};

//! Function to create a gravity field model.
/*!
 *  Function to create a gravity field model based on model-specific settings for the gravity field.
 *  \param gravityFieldSettings Settings for the gravity field model that is to be created, defined
 *  a pointer to an object of class (derived from) GravityFieldSettings.
 *  \param body Name of the body for which the gravity field model is to be created.
 *  \return Gravity field model created according to settings in gravityFieldSettings.
 */
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body);
}

}
#endif // TUDAT_CREATEGRAVITYFIELD_H
