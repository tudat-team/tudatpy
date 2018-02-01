/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h"
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

    //! Function to reset identifier for body-fixed reference frame to which the coefficients are referred.
    /*!
     *  Function to reset identifier for body-fixed reference frame to which the coefficients are referred.
     *  \param associatedReferenceFrame Identifier for body-fixed reference frame to which the coefficients are referred.
     */
    void resetAssociatedReferenceFrame( const std::string& associatedReferenceFrame )
    {
        associatedReferenceFrame_ = associatedReferenceFrame;
    }

    //! Function to retrieve boolean that denotes whether the field should be created as time-dependent
    /*!
     *  Function to retrieve boolean that denotes whether the field should be created as time-dependent
     *  \return Boolean that denotes whether the field should be created as time-dependent
     */
    bool getCreateTimeDependentField( )
    {
        return createTimeDependentField_;
    }

    //! Function to reset boolean that denotes whether the field should be created as time-dependent
    /*!
     *  Function to reset boolean that denotes whether the field should be created as time-dependent
     *  \param createTimeDependentField Boolean that denotes whether the field should be created as time-dependent
     */
    void setCreateTimeDependentField( const bool createTimeDependentField )
    {
        createTimeDependentField_ = createTimeDependentField;
    }

protected:


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

    //! Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed intially)
    bool createTimeDependentField_;

};


//! Spherical harmonics models supported by Tudat.
enum SphericalHarmonicsModel
{
    customModel,
    egm96,
    ggm02c,
    ggm02s,
    glgm3150,
    lpe200,
    jgmro120d
};

//! Get the path of the SH file for a SH model.
/*!
 * @copybrief getPathForSphericalHarmonicsModel
 * \param sphericalHarmonicsModel The spherical harmonics model.
 * \return The path of the SH file for a SH model.
 */
std::string getPathForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel );

//! Get the associated reference frame for a SH model.
/*!
 * @copybrief getReferenceFrameForSphericalHarmonicsModel
 * \param sphericalHarmonicsModel The spherical harmonics model.
 * \return The associated reference frame for a SH model.
 */
std::string getReferenceFrameForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel );

//! Derived class of SphericalHarmonicsGravityFieldSettings defining settings of spherical harmonic gravity
//! field representation to be loaded from a spherical harmonics model file.
class FromFileSphericalHarmonicsGravityFieldSettings: public SphericalHarmonicsGravityFieldSettings
{
public:
    //! Constructor with custom model.
    /*!
     * Constructor with custom model.
     * \param filePath Path of PDS gravity field file to be loaded.
     * \param associatedReferenceFrame Identifier for body-fixed reference frame to which the coefficients are referred.
     * \param maximumDegree Maximum degree of gravity field to be loaded.
     * \param maximumOrder Maximum order of gravity field to be loaded.
     * \param gravitationalParameterIndex Index at which the gravitational parameter can be found in the header
     * (first line of the file). Set to -1 if the file has no header.
     * \param referenceRadiusIndex Index at which the reference radius can be found in the header
     * (first line of the file). Set to -1 if the file has no header.
     * \param gravitationalParameter Gravitational parameter of gravity field to be used if file has no header.
     * \param referenceRadius Reference radius of gravity field to be used if file has no header.
     */
    FromFileSphericalHarmonicsGravityFieldSettings( const std::string& filePath,
                                                 const std::string& associatedReferenceFrame,
                                                 const int maximumDegree,
                                                 const int maximumOrder,
                                                 const int gravitationalParameterIndex,
                                                 const int referenceRadiusIndex,
                                                 const double gravitationalParameter = TUDAT_NAN,
                                                 const double referenceRadius = TUDAT_NAN );

    //! Constructor with model included in Tudat.
    /*!
     * Constructor with model included in Tudat.
     * \param sphericalHarmonicsModel Spherical harmonics model to be used.
     */
    FromFileSphericalHarmonicsGravityFieldSettings( const SphericalHarmonicsModel sphericalHarmonicsModel );

    //! Get the sphericals harmonics model.
    /*!
     * @copybrief getSphericalHarmonicsModel
     * \return The sphericals harmonics model.
     */
    SphericalHarmonicsModel getSphericalHarmonicsModel( )
    {
        return sphericalHarmonicsModel_;
    }

    //! Get the sphericals harmonics model.
    /*!
     * @copybrief getSphericalHarmonicsModel
     * \return The sphericals harmonics model.
     */
    std::string getFilePath( )
    {
        return filePath_;
    }

    //! Get the maximum degree.
    /*!
     * @copybrief getMaximumDegree
     * \return The maximum degree.
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Get the maximum order.
    /*!
     * @copybrief getMaximumOrder
     * \return The maximum order.
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

    //! Get the gravitational parameter index.
    /*!
     * @copybrief getGravitationalParameterIndex
     * \return The gravitational parameter index.
     */
    int getGravitationalParameterIndex( )
    {
        return gravitationalParameterIndex_;
    }

    //! Get the reference radius index.
    /*!
     * @copybrief getReferenceRadiusIndex
     * \return The reference radius index.
     */
    int getReferenceRadiusIndex( )
    {
        return referenceRadiusIndex_;
    }

protected:
    //! Spherical harmonics model.
    SphericalHarmonicsModel sphericalHarmonicsModel_ = customModel;

    //! Path of loaded PDS gravity field file.
    std::string filePath_;

    //! Maximum loaded degree from file.
    int maximumDegree_;

    //! Maximum loaded order from file.
    int maximumOrder_;

    //! Index at which the gravitational parameter can be found in the first line of the file.
    //! -1 if this information is not available in the file.
    int gravitationalParameterIndex_;

    //! Index at which the reference radius can be found in the first line of the file.
    //! -1 if this information is not available in the file.
    int referenceRadiusIndex_;

};

//! Function to create gravity field settings for a homogeneous triaxial ellipsoid
/*!
 * Function to create gravity field settings for a homogeneous triaxial ellipsoid. The gravity field is expressed in
 * normalized spherical harmonic coefficients.  X-axis is alligned
 * with largest axis, y-axis with middle axis and z-axis with smallest axis
 * \param axisA Largest axis of triaxial ellipsoid
 * \param axisB Middle axis of triaxial ellipsoid
 * \param axisC Smallest axis of triaxial ellipsoid
 * \param ellipsoidDensity Density of ellipsoid.
 * \param maximumDegree Maximum degree of expansion
 * \param maximumOrder Maximum oredr of expansion
 * \param associatedReferenceFrame Identifier for body-fixed reference frame to which
 * the coefficients are referred.
 * \return Gravity field settings for a homogeneous triaxial ellipsoid of given properties.
 */
boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA, const double axisB, const double axisC, const double ellipsoidDensity,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedReferenceFrame  );

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
