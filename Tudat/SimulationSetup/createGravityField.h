/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGRAVITYFIELD_H
#define TUDAT_CREATEGRAVITYFIELD_H

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
#include "Tudat/SimulationSetup/createGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"

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
                                            const std::string& associatedReferenceFrame ):
        GravityFieldSettings( spherical_harmonic ),
        gravitationalParameter_( gravitationalParameter ),
        referenceRadius_( referenceRadius ),
        cosineCoefficients_( cosineCoefficients ),
        sineCoefficients_( sineCoefficients ),
        associatedReferenceFrame_( associatedReferenceFrame ),
        createTimeDependentField_( 0 )
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

    bool getCreateTimeDependentField( ){ return createTimeDependentField_; }

    void setCreateTimeDependentField( const bool createTimeDependentField )
    {
        createTimeDependentField_ = createTimeDependentField;
    }

private:


    //! Gravitational parameter for gravity field that is to be created.
    double gravitationalParameter_;

    //! Reference radius of spherical harmonic field expansion.
    double referenceRadius_;

    //! Cosine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd cosineCoefficients_;

    //! Sine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd sineCoefficients_;

    //! Identifier for body-fixed reference frame to which the coefficients are referred.
    std::string associatedReferenceFrame_;

    bool createTimeDependentField_;

};

//! Function to read a spherical harmonic gravity field file
/*!
 *  Function to read a spherical harmonic gravity field file, returns (by reference) cosine and sine
 *  spherical harmomic coefficients.
 *  The file structure should be as follows: The first line may be a file header with metadata. If this is the case,
 *  the gravitationalParameterIndex and referenceRadiusIndex should be used as input, to communicate which entries in
 *  the list of metadata represents these quantities. If both these variables are NaN, the file is assumed to have no
 *  header. The following lines of the file must have the following structure:
 *  Degree, Order, Cosine Coefficient, Sine Coefficients
 *  Subsequent columns may be present in the file, but are ignored when parsing.
 *  All coefficients not defined in the file are set to zero (except C(0,0) which is always 1.0)
 *  \param fileName Name of PDS gravity field file to be loaded.
 *  \param maximumDegree Maximum degree of gravity field to be loaded.
 *  \param maximumOrder Maximum order of gravity field to be loaded.
 *  \param gravitationalParameterIndex
 *  \param referenceRadiusIndex
 *  \param coefficients Spherical harmonics coefficients (first is cosine, second is sine).
 *  \return Pair of gravitational parameter and reference radius, values are non-NaN if
 *  gravitationalParameterIndex and referenceRadiusIndex are >=0.
 */
std::pair< double, double > readGravityFieldFile(
        const std::string& fileName, const int maximumDegree, const int maximumOrder,
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd >& coefficients,
        const int gravitationalParameterIndex = -1, const int referenceRadiusIndex = -1 );

//! Function to create a gravity field model.
/*!
 *  Function to create a gravity field model based on model-specific settings for the gravity field.
 *  \param gravityFieldSettings Settings for the gravity field model that is to be created, defined
 *  a pointer to an object of class (derived from) GravityFieldSettings.
 *  \param body Name of the body for which the gravity field model is to be created.
 *  \param bodyMap List of body objects, as currently created (used when setting
 *  gravityFieldVariationSettings)
 *  \param gravityFieldVariationSettings List of settings for the variations of the gravity field
 *  that are to be used (but not immediately set!) by current body under consideration.
 *  \return Gravity field model created according to settings in gravityFieldSettings.
 */
boost::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const boost::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const NamedBodyMap& bodyMap = NamedBodyMap( ),
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings =
        std::vector< boost::shared_ptr< GravityFieldVariationSettings > >( ) );

} // namespace simulation_setup

} // namespace tudat
#endif // TUDAT_CREATEGRAVITYFIELD_H
