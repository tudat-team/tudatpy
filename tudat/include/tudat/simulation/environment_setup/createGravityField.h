/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <memory>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/astro/gravitation/ringGravityField.h"

namespace tudat
{

namespace simulation_setup
{

// List of gravity field models available in simulations
/*
 *  List of gravity field models available in simulations. Gravity field models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
//! @get_docstring(GravityFieldType.__docstring__)
enum GravityFieldType { central, central_spice, spherical_harmonic, polyhedron, one_dimensional_ring };

// Class for providing settings for gravity field model.
/*
 *  Class for providing settings for automatic gravity field model creation. This class is a
 *  functional (base) class for settings of gravity field models that require no information in
 *  addition to their type. Gravity field model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */

//! @get_docstring(GravityFieldSettings.__docstring__)
class GravityFieldSettings
{
public:
    // Constructor, sets type of gravity field model.
    /*
     *  Constructor, sets type of gravity field model. Settings for gravity field models requiring
     *  additional information should be defined in a derived class.
     *  \param gravityFieldType Type of gravity field model that is to be created.
     */
    GravityFieldSettings( const GravityFieldType gravityFieldType ): gravityFieldType_( gravityFieldType ) { }

    // Destructor
    virtual ~GravityFieldSettings( ) { }

    // Function to return type of gravity field model that is to be created.
    /*
     *  Function to return type of gravity field model that is to be created.
     *  \return Type of gravity field model that is to be created.
     */
    GravityFieldType getGravityFieldType( )
    {
        return gravityFieldType_;
    }

protected:
    // Type of gravity field model that is to be created.
    GravityFieldType gravityFieldType_;
};

// Derived class of GravityFieldSettings defining settings of point mass gravity field.
//! @get_docstring(CentralGravityFieldSettings.__docstring__)
class CentralGravityFieldSettings : public GravityFieldSettings
{
public:
    // Constructor.
    /*
     *  Constructor.
     *  \param gravitationalParameter Gravitational parameter of gravity field.
     */
    CentralGravityFieldSettings( double gravitationalParameter ):
        GravityFieldSettings( central ), gravitationalParameter_( gravitationalParameter )
    { }

    virtual ~CentralGravityFieldSettings( ) { }

    // Function to return gravitational parameter for gravity field.
    /*
     *  Function to return gravitational parameter for gravity field.
     *  \return Gravitational parameter for gravity field.
     */
    double getGravitationalParameter( )
    {
        return gravitationalParameter_;
    }

    void resetGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
    }

private:
    // Gravitational parameter for gravity field.
    double gravitationalParameter_;
};

// Derived class of GravityFieldSettings defining settings of spherical harmonic gravity
// field representation.
//! @get_docstring(SphericalHarmonicsGravityFieldSettings.__docstring__)
class SphericalHarmonicsGravityFieldSettings : public GravityFieldSettings
{
public:
    // Constructor.
    /*
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
                                            const Eigen::MatrixXd& cosineCoefficients,
                                            const Eigen::MatrixXd& sineCoefficients,
                                            const std::string& associatedReferenceFrame,
                                            const double scaledMeanMomentOfInertia = TUDAT_NAN ):
        GravityFieldSettings( spherical_harmonic ), gravitationalParameter_( gravitationalParameter ), referenceRadius_( referenceRadius ),
        inertiaTensor_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ), cosineCoefficients_( cosineCoefficients ),
        sineCoefficients_( sineCoefficients ), associatedReferenceFrame_( associatedReferenceFrame ), createTimeDependentField_( 0 ),
        scaledMeanMomentOfInertia_( scaledMeanMomentOfInertia )
    { }

    SphericalHarmonicsGravityFieldSettings( const double gravitationalParameter,
                                            const double referenceRadius,
                                            const Eigen::Matrix3d& inertiaTensor,
                                            const Eigen::MatrixXd& cosineCoefficients,
                                            const Eigen::MatrixXd& sineCoefficients,
                                            const std::string& associatedReferenceFrame ):
        GravityFieldSettings( spherical_harmonic ), gravitationalParameter_( gravitationalParameter ), referenceRadius_( referenceRadius ),
        inertiaTensor_( inertiaTensor ), cosineCoefficients_( cosineCoefficients ), sineCoefficients_( sineCoefficients ),
        associatedReferenceFrame_( associatedReferenceFrame ), createTimeDependentField_( 0 )
    {
        std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > degreeTwoField =
                gravitation::getDegreeTwoSphericalHarmonicCoefficients( inertiaTensor_, gravitationalParameter_, referenceRadius_ );

        cosineCoefficients_.block( 2, 0, 1, 3 ) = std::get< 0 >( degreeTwoField ).block( 2, 0, 1, 3 );
        sineCoefficients_.block( 2, 1, 1, 2 ) = std::get< 1 >( degreeTwoField ).block( 2, 1, 1, 2 );
        scaledMeanMomentOfInertia_ = std::get< 2 >( degreeTwoField );
    }

    virtual ~SphericalHarmonicsGravityFieldSettings( ) { }

    // Function to return gravitational parameter for gravity field.
    /*
     *  Function to return gravitational parameter for gravity field.
     *  \return Gravitational parameter for gravity field.
     */
    double getGravitationalParameter( )
    {
        return gravitationalParameter_;
    }

    // Function to reset gravitational parameter for gravity field.
    /*
     *  Function to reset gravitational parameter for gravity field.
     *  \param gravitationalParameter New gravitational parameter for gravity field.
     */
    void resetGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
    }

    // Function to return reference radius of spherical harmonic field expansion
    /*
     *  Function to return reference radius of spherical harmonic field expansion
     *  \return Reference radius of spherical harmonic field expansion
     */
    double getReferenceRadius( )
    {
        return referenceRadius_;
    }

    // Function to return cosine spherical harmonic coefficients (geodesy normalized).
    /*
     *  Function to return cosine spherical harmonic coefficients (geodesy normalized).
     *  \return Cosine spherical harmonic coefficients (geodesy normalized).
     */
    Eigen::MatrixXd getCosineCoefficients( )
    {
        return cosineCoefficients_;
    }

    Eigen::Matrix3d getInertiaTensor( )
    {
        return inertiaTensor_;
    }

    double getScaledMeanMomentOfInertia( )
    {
        return scaledMeanMomentOfInertia_;
    }

    void setScaledMeanMomentOfInertia( const double scaledMeanMomentOfInertia )
    {
        scaledMeanMomentOfInertia_ = scaledMeanMomentOfInertia;
    }

    void resetCosineCoefficients( const Eigen::MatrixXd cosineCoefficients )
    {
        cosineCoefficients_ = cosineCoefficients;
    }

    // Function to return sine spherical harmonic coefficients (geodesy normalized).
    /*
     *  Function to return sine spherical harmonic coefficients (geodesy normalized).
     *  \return Sine spherical harmonic coefficients (geodesy normalized).
     */
    Eigen::MatrixXd getSineCoefficients( )
    {
        return sineCoefficients_;
    }

    void resetSineCoefficients( const Eigen::MatrixXd sineCoefficients )
    {
        sineCoefficients_ = sineCoefficients;
    }

    // Function to return identifier for body-fixed reference frame.
    /*
     *  Function to return identifier for body-fixed reference frame to which the coefficients
     *  are referred.
     *  \return Identifier for body-fixed reference frame to which the coefficients are referred.
     */
    std::string getAssociatedReferenceFrame( )
    {
        return associatedReferenceFrame_;
    }

    // Function to reset identifier for body-fixed reference frame to which the coefficients are referred.
    /*
     *  Function to reset identifier for body-fixed reference frame to which the coefficients are referred.
     *  \param associatedReferenceFrame Identifier for body-fixed reference frame to which the coefficients are referred.
     */
    void resetAssociatedReferenceFrame( const std::string& associatedReferenceFrame )
    {
        associatedReferenceFrame_ = associatedReferenceFrame;
    }

    // Function to retrieve boolean that denotes whether the field should be created as time-dependent
    /*
     *  Function to retrieve boolean that denotes whether the field should be created as time-dependent
     *  \return Boolean that denotes whether the field should be created as time-dependent
     */
    bool getCreateTimeDependentField( )
    {
        return createTimeDependentField_;
    }

    // Function to reset boolean that denotes whether the field should be created as time-dependent
    /*
     *  Function to reset boolean that denotes whether the field should be created as time-dependent
     *  \param createTimeDependentField Boolean that denotes whether the field should be created as time-dependent
     */
    void setCreateTimeDependentField( const bool createTimeDependentField )
    {
        createTimeDependentField_ = createTimeDependentField;
    }

protected:
    // Gravitational parameter for gravity field that is to be created.
    double gravitationalParameter_;

    // Reference radius of spherical harmonic field expansion.
    double referenceRadius_;

    Eigen::Matrix3d inertiaTensor_;

    // Cosine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd cosineCoefficients_;

    // Sine spherical harmonic coefficients (geodesy normalized).
    Eigen::MatrixXd sineCoefficients_;

    // Identifier for body-fixed reference frame to which the coefficients are referred.
    std::string associatedReferenceFrame_;

    // Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed intially)
    bool createTimeDependentField_;

    double scaledMeanMomentOfInertia_;
};

// Derived class of GravityFieldSettings defining settings of polyhedron gravity
// field representation.
// References: "EXTERIOR GRAVITATION OF A POLYHEDRON DERIVED AND COMPARED WITH HARMONIC AND MASCON GRAVITATION REPRESENTATIONS
//              OF ASTEROID 4769 CASTALIA", Werner and Scheeres (1997), Celestial Mechanics and Dynamical Astronomy
class PolyhedronGravityFieldSettings : public GravityFieldSettings
{
public:
    // Constructor.
    /*
     *  Constructor.
     *  \param gravitationalConstant Gravitational constant of the gravity field.
     *  \param density Density of polyhedron.
     *  \param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
     *  \param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
     *  \param associatedReferenceFrame Identifier for body-fixed reference frame to which the polyhedron is referred.
     */
    PolyhedronGravityFieldSettings( const double gravitationalConstant,
                                    const double density,
                                    const Eigen::MatrixXd& verticesCoordinates,
                                    const Eigen::MatrixXi& verticesDefiningEachFacet,
                                    const std::string& associatedReferenceFrame ):
        GravityFieldSettings( polyhedron ), density_( density ), verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ), associatedReferenceFrame_( associatedReferenceFrame ),
        gravitationalConstant_( gravitationalConstant )
    {
        volume_ = basic_astrodynamics::computePolyhedronVolume( verticesCoordinates, verticesDefiningEachFacet );
        gravitationalParameter_ = gravitationalConstant_ * density_ * volume_;
    }

    /*! Constructor.
     *
     * Constructor.
     * @param gravitationalParameter Gravitational parameter of the polyhedron.
     * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
     * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
     * @param associatedReferenceFrame Identifier for body-fixed reference frame to which the polyhedron is referred.
     * @param density Density of the polyhedron. It is later used to set the function to update the inertia tensor.
     */
    PolyhedronGravityFieldSettings( const double gravitationalParameter,
                                    const Eigen::MatrixXd& verticesCoordinates,
                                    const Eigen::MatrixXi& verticesDefiningEachFacet,
                                    const std::string& associatedReferenceFrame,
                                    const double gravitationalConstant = physical_constants::GRAVITATIONAL_CONSTANT ):
        GravityFieldSettings( polyhedron ), gravitationalParameter_( gravitationalParameter ), verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ), associatedReferenceFrame_( associatedReferenceFrame ),
        gravitationalConstant_( gravitationalConstant )
    {
        volume_ = basic_astrodynamics::computePolyhedronVolume( verticesCoordinates, verticesDefiningEachFacet );
        density_ = gravitationalParameter_ / ( gravitationalConstant_ * volume_ );
    }

    //! Destructor
    virtual ~PolyhedronGravityFieldSettings( ) { }

    // Function to return the gravitational parameter.
    double getGravitationalParameter( )
    {
        return gravitationalParameter_;
    }

    // Function to reset the gravitational parameter.
    void resetGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
        density_ = gravitationalParameter_ / ( gravitationalConstant_ * volume_ );
    }

    // Function to return the density.
    double getDensity( )
    {
        return density_;
    }

    // Function to reset the density.
    void resetDensity( double density )
    {
        density_ = density;
        gravitationalParameter_ = gravitationalConstant_ * density_ * volume_;
    }

    // Function to return identifier for body-fixed reference frame.
    std::string getAssociatedReferenceFrame( )
    {
        return associatedReferenceFrame_;
    }

    // Function to reset identifier for body-fixed reference frame to which the polyhedron is referred.
    void resetAssociatedReferenceFrame( const std::string& associatedReferenceFrame )
    {
        associatedReferenceFrame_ = associatedReferenceFrame;
    }

    //! Function to return the vertices coordinates.
    const Eigen::MatrixXd& getVerticesCoordinates( )
    {
        return verticesCoordinates_;
    }

    //! Function to reset the vertices coordinates.
    void resetVerticesCoordinates( const Eigen::MatrixXd& verticesCoordinates )
    {
        verticesCoordinates_ = verticesCoordinates;
    }

    //! Function to return the vertices defining each facet.
    const Eigen::MatrixXi& getVerticesDefiningEachFacet( )
    {
        return verticesDefiningEachFacet_;
    }

    //! Function to reset the vertices defining each facet.
    void resetVerticesDefiningEachFacet( const Eigen::MatrixXi& verticesDefiningEachFacet )
    {
        verticesDefiningEachFacet_ = verticesDefiningEachFacet;
    }

protected:
    // Gravitational parameter
    double gravitationalParameter_;

    // Density of the polyhedron
    double density_;

    // Cartesian coordinates of each vertex.
    Eigen::MatrixXd verticesCoordinates_;

    // Indices of the 3 vertices describing each facet.
    Eigen::MatrixXi verticesDefiningEachFacet_;

    // Identifier for body-fixed reference frame to which the polyhedron is referred.
    std::string associatedReferenceFrame_;

    double gravitationalConstant_;

    double volume_;
};

// Derived class of GravityFieldSettings defining settings of polyhedron gravity
// field representation.
// References: Precise computation of acceleration due to uniform ring or disk, Toshio Fukushima (2010), Celestial Mechanics
//             and Dynamical Astronomy, 108:339–356.
class RingGravityFieldSettings : public GravityFieldSettings
{
public:
    /*! Constructor.
     *
     * Constructor.
     * @param gravitationalParameter Gravitational parameter of the ring.
     * @param ringRadius Radius of the ring.
     * @param associatedReferenceFrame Identifier for body-fixed reference frame to which the ring is referred.
     * @param ellipticIntegralSFromDAndB Flag indicating whether to compute S(m) from D(m) and B(m) (if true),
     *      or from K(m) and E(m) (if false). The former has a lower loss of accuracy due to numerical cancellation.
     */
    RingGravityFieldSettings( const double gravitationalParameter,
                              const double ringRadius,
                              const std::string& associatedReferenceFrame,
                              const bool ellipticIntegralSFromDAndB = true ):
        GravityFieldSettings( one_dimensional_ring ), gravitationalParameter_( gravitationalParameter ), ringRadius_( ringRadius ),
        associatedReferenceFrame_( associatedReferenceFrame ), ellipticIntegralSFromDAndB_( ellipticIntegralSFromDAndB )
    { }

    //! Destructor
    virtual ~RingGravityFieldSettings( ) { }

    // Function to return the gravitational parameter.
    double getGravitationalParameter( )
    {
        return gravitationalParameter_;
    }

    // Function to reset the gravitational parameter.
    void resetGravitationalParameter( const double gravitationalParameter )
    {
        gravitationalParameter_ = gravitationalParameter;
    }

    // Function to return the ring radius
    double getRingRadius( )
    {
        return ringRadius_;
    }

    // Function to reset the density.
    void resetRingRadius( double ringRadius )
    {
        ringRadius_ = ringRadius;
    }

    // Function to return identifier for body-fixed reference frame.
    std::string getAssociatedReferenceFrame( )
    {
        return associatedReferenceFrame_;
    }

    // Function to reset identifier for body-fixed reference frame to which the ring is referred.
    void resetAssociatedReferenceFrame( const std::string& associatedReferenceFrame )
    {
        associatedReferenceFrame_ = associatedReferenceFrame;
    }

    // Function to get the flag indicating whether to compute S(m) from D(m) and B(m)
    bool getEllipticIntegralSFromDAndB( )
    {
        return ellipticIntegralSFromDAndB_;
    }

    // Function to reset the flag indicating whether to compute S(m) from D(m) and B(m)
    void resetEllipticIntegralSFromDAndB( bool ellipticIntegralSFromDAndB )
    {
        ellipticIntegralSFromDAndB_ = ellipticIntegralSFromDAndB;
    }

protected:
    // Gravitational parameter
    double gravitationalParameter_;

    // Radius of the ring
    double ringRadius_;

    // Identifier for body-fixed reference frame to which the ring is referred
    std::string associatedReferenceFrame_;

    // Flag indicating whether to compute S(m) from D(m) and B(m) (if true), or from K(m) and E(m) (if false)
    bool ellipticIntegralSFromDAndB_;
};

// Spherical harmonics models supported by Tudat.
//! @get_docstring(SphericalHarmonicsModel.__docstring__)
enum SphericalHarmonicsModel { customModel, egm96, ggm02c, ggm02s, goco05c, glgm3150, lpe200, gggrx1200, jgmro120d, jgmess160a, shgj180u };

// Get the path of the SH file for a SH model.
/*
 * @copybrief getPathForSphericalHarmonicsModel
 * \param sphericalHarmonicsModel The spherical harmonics model.
 * \return The path of the SH file for a SH model.
 */
std::string getPathForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel );

int getMaximumGravityFieldDegreeOrder( const SphericalHarmonicsModel sphericalHarmonicsModel );

// Get the associated reference frame for a SH model.
/*
 * @copybrief getReferenceFrameForSphericalHarmonicsModel
 * \param sphericalHarmonicsModel The spherical harmonics model.
 * \return The associated reference frame for a SH model.
 */
std::string getReferenceFrameForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel );

// Derived class of SphericalHarmonicsGravityFieldSettings defining settings of spherical harmonic gravity
// field representation to be loaded from a spherical harmonics model file.
class FromFileSphericalHarmonicsGravityFieldSettings : public SphericalHarmonicsGravityFieldSettings
{
public:
    // Constructor with custom model.
    /*
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
    // Constructor with model included in Tudat.
    /*
     * Constructor with model included in Tudat.
     * \param sphericalHarmonicsModel Spherical harmonics model to be used.
     */
    FromFileSphericalHarmonicsGravityFieldSettings( const SphericalHarmonicsModel sphericalHarmonicsModel, const int maximumDegree = -1 );

    virtual ~FromFileSphericalHarmonicsGravityFieldSettings( ) { }
    // Get the sphericals harmonics model.
    /*
     * @copybrief getSphericalHarmonicsModel
     * \return The sphericals harmonics model.
     */
    SphericalHarmonicsModel getSphericalHarmonicsModel( )
    {
        return sphericalHarmonicsModel_;
    }

    // Get the sphericals harmonics model.
    /*
     * @copybrief getSphericalHarmonicsModel
     * \return The sphericals harmonics model.
     */
    std::string getFilePath( )
    {
        return filePath_;
    }

    // Get the maximum degree.
    /*
     * @copybrief getMaximumDegree
     * \return The maximum degree.
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    // Get the maximum order.
    /*
     * @copybrief getMaximumOrder
     * \return The maximum order.
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

    // Get the gravitational parameter index.
    /*
     * @copybrief getGravitationalParameterIndex
     * \return The gravitational parameter index.
     */
    int getGravitationalParameterIndex( )
    {
        return gravitationalParameterIndex_;
    }

    // Get the reference radius index.
    /*
     * @copybrief getReferenceRadiusIndex
     * \return The reference radius index.
     */
    int getReferenceRadiusIndex( )
    {
        return referenceRadiusIndex_;
    }

protected:
    // Spherical harmonics model.
    SphericalHarmonicsModel sphericalHarmonicsModel_ = customModel;

    // Path of loaded PDS gravity field file.
    std::string filePath_;

    // Maximum loaded degree from file.
    int maximumDegree_;

    // Maximum loaded order from file.
    int maximumOrder_;

    // Index at which the gravitational parameter can be found in the first line of the file.
    // -1 if this information is not available in the file.
    int gravitationalParameterIndex_;

    // Index at which the reference radius can be found in the first line of the file.
    // -1 if this information is not available in the file.
    int referenceRadiusIndex_;
};

// Function to create gravity field settings for a homogeneous triaxial ellipsoid with density as argument.
/*
 * Function to create gravity field settings for a homogeneous triaxial ellipsoid with density as argument. The gravity
 * field is expressed in normalized spherical harmonic coefficients.  X-axis is alligned
 * with largest axis, y-axis with middle axis and z-axis with smallest axis
 * \param axisA Largest axis of triaxial ellipsoid
 * \param axisB Middle axis of triaxial ellipsoid
 * \param axisC Smallest axis of triaxial ellipsoid
 * \param ellipsoidDensity Density of ellipsoid.
 * \param maximumDegree Maximum degree of expansion
 * \param maximumOrder Maximum oredr of expansion
 * \param associatedReferenceFrame Identifier for body-fixed reference frame to which
 * the coefficients are referred.
 * \param gravitationalConstant Gravitational constant.
 * \return Gravity field settings for a homogeneous triaxial ellipsoid of given properties.
 */
std::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA,
        const double axisB,
        const double axisC,
        const double ellipsoidDensity,
        const int maximumDegree,
        const int maximumOrder,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant = physical_constants::GRAVITATIONAL_CONSTANT );

// Function to create gravity field settings for a homogeneous triaxial ellipsoid with the gravitational parameter as argument.
/*
 * Function to create gravity field settings for a homogeneous triaxial ellipsoid with the gravitational parameters as
 * argument. The gravity field is expressed in normalized spherical harmonic coefficients.  X-axis is alligned
 * with largest axis, y-axis with middle axis and z-axis with smallest axis
 * \param axisA Largest axis of triaxial ellipsoid
 * \param axisB Middle axis of triaxial ellipsoid
 * \param axisC Smallest axis of triaxial ellipsoid
 * \param ellipsoidGravitationalParameter Gravitational parameter of the gravity field.
 * \param maximumDegree Maximum degree of expansion
 * \param maximumOrder Maximum oredr of expansion
 * \param associatedReferenceFrame Identifier for body-fixed reference frame to which
 * the coefficients are referred.
 * \return Gravity field settings for a homogeneous triaxial ellipsoid of given properties.
 */
std::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA,
        const double axisB,
        const double axisC,
        const int maximumDegree,
        const int maximumOrder,
        const std::string& associatedReferenceFrame,
        const double ellipsoidGravitationalParameter );

// Function to read a spherical harmonic gravity field file
/*
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
std::pair< double, double > readGravityFieldFile( const std::string& fileName,
                                                  const int maximumDegree,
                                                  const int maximumOrder,
                                                  std::pair< Eigen::MatrixXd, Eigen::MatrixXd >& coefficients,
                                                  const int gravitationalParameterIndex = -1,
                                                  const int referenceRadiusIndex = -1 );

// Function to create a gravity field model.
/*
 *  Function to create a gravity field model based on model-specific settings for the gravity field.
 *  \param gravityFieldSettings Settings for the gravity field model that is to be created, defined
 *  a pointer to an object of class (derived from) GravityFieldSettings.
 *  \param body Name of the body for which the gravity field model is to be created.
 *  \param bodies List of body objects, as currently created (used when setting
 *  gravityFieldVariationSettings)
 *  \param gravityFieldVariationSettings List of settings for the variations of the gravity field
 *  that are to be used (but not immediately set!) by current body under consideration.
 *  \return Gravity field model created according to settings in gravityFieldSettings.
 */
std::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const std::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const SystemOfBodies& bodies = SystemOfBodies( ),
        const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings =
                std::vector< std::shared_ptr< GravityFieldVariationSettings > >( ) );

//! @get_docstring(centralGravitySettings)
inline std::shared_ptr< GravityFieldSettings > centralGravitySettings( const double gravitationalParameter )
{
    return std::make_shared< CentralGravityFieldSettings >( gravitationalParameter );
}

//! @get_docstring(centralGravityFromSpiceSettings)
inline std::shared_ptr< GravityFieldSettings > centralGravityFromSpiceSettings( )
{
    return std::make_shared< GravityFieldSettings >( central_spice );
}

//! @get_docstring(sphericalHarmonicsGravitySettings)
inline std::shared_ptr< GravityFieldSettings > sphericalHarmonicsGravitySettings(
        const double gravitationalParameter,
        const double referenceRadius,
        const Eigen::MatrixXd normalizedCosineCoefficients,  // NOTE: entry (i,j) denotes coefficient at degree i and order j
        const Eigen::MatrixXd normalizedSineCoefficients,    // NOTE: entry (i,j) denotes coefficient at degree i and order j
        const std::string& associatedReferenceFrame )
{
    return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
            gravitationalParameter, referenceRadius, normalizedCosineCoefficients, normalizedSineCoefficients, associatedReferenceFrame );
}

inline std::shared_ptr< GravityFieldSettings > sphericalHarmonicsGravitySettings(
        const double gravitationalParameter,
        const double referenceRadius,
        const Eigen::Matrix3d& inertiaTensor,
        const Eigen::MatrixXd normalizedCosineCoefficients,  // NOTE: entry (i,j) denotes coefficient at degree i and order j
        const Eigen::MatrixXd normalizedSineCoefficients,    // NOTE: entry (i,j) denotes coefficient at degree i and order j
        const std::string& associatedReferenceFrame )
{
    return std::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter,
                                                                       referenceRadius,
                                                                       inertiaTensor,
                                                                       normalizedCosineCoefficients,
                                                                       normalizedSineCoefficients,
                                                                       associatedReferenceFrame );
}

inline std::shared_ptr< GravityFieldSettings > polyhedronGravitySettings(
        const double density,
        const Eigen::MatrixXd verticesCoordinates,
        const Eigen::MatrixXi verticesDefiningEachFacet,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant = physical_constants::GRAVITATIONAL_CONSTANT )
{
    return std::make_shared< PolyhedronGravityFieldSettings >(
            gravitationalConstant, density, verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame );
}

inline std::shared_ptr< GravityFieldSettings > polyhedronGravitySettingsFromMu(
        const double gravitationalParameter,
        const Eigen::MatrixXd verticesCoordinates,
        const Eigen::MatrixXi verticesDefiningEachFacet,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant = physical_constants::GRAVITATIONAL_CONSTANT )
{
    return std::make_shared< PolyhedronGravityFieldSettings >(
            gravitationalParameter, verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame, gravitationalConstant );
}

inline std::shared_ptr< GravityFieldSettings > ringGravitySettings( const double gravitationalParameter,
                                                                    const double ringRadius,
                                                                    const std::string& associatedReferenceFrame,
                                                                    const bool ellipticIntegralSFromDAndB = true )
{
    return std::make_shared< RingGravityFieldSettings >(
            gravitationalParameter, ringRadius, associatedReferenceFrame, ellipticIntegralSFromDAndB );
}

enum RigidBodyPropertiesType {
    constant_rigid_body_properties,
    from_function_rigid_body_properties,
    from_gravity_field_rigid_body_properties,
    mass_dependent_rigid_body_properties
};

class RigidBodyPropertiesSettings
{
public:
    RigidBodyPropertiesSettings( const RigidBodyPropertiesType rigidBodyPropertiesType ):
        rigidBodyPropertiesType_( rigidBodyPropertiesType )
    { }

    virtual ~RigidBodyPropertiesSettings( ) { }

    RigidBodyPropertiesType getRigidBodyPropertiesType( )
    {
        return rigidBodyPropertiesType_;
    }

protected:
    RigidBodyPropertiesType rigidBodyPropertiesType_;
};

class ConstantRigidBodyPropertiesSettings : public RigidBodyPropertiesSettings
{
public:
    ConstantRigidBodyPropertiesSettings( const double mass,
                                         const Eigen::Vector3d& centerOfMass = Eigen::Vector3d::Constant( TUDAT_NAN ),
                                         const Eigen::Matrix3d& inertiaTensor = Eigen::Matrix3d::Constant( TUDAT_NAN ) ):
        RigidBodyPropertiesSettings( constant_rigid_body_properties ), mass_( mass ), centerOfMass_( centerOfMass ),
        inertiaTensor_( inertiaTensor )
    { }

    virtual ~ConstantRigidBodyPropertiesSettings( ) { }

    double getMass( )
    {
        return mass_;
    }

    Eigen::Vector3d getCenterOfMass( )
    {
        return centerOfMass_;
    }

    Eigen::Matrix3d getInertiaTensor( )
    {
        return inertiaTensor_;
    }

protected:
    double mass_;

    Eigen::Vector3d centerOfMass_;

    Eigen::Matrix3d inertiaTensor_;
};

class FromFunctionRigidBodyPropertiesSettings : public RigidBodyPropertiesSettings
{
public:
    FromFunctionRigidBodyPropertiesSettings( const std::function< double( const double ) > massFunction,
                                             const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction = nullptr,
                                             const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction = nullptr ):
        RigidBodyPropertiesSettings( from_function_rigid_body_properties ), massFunction_( massFunction ),
        centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction )
    { }

    virtual ~FromFunctionRigidBodyPropertiesSettings( ) { }

    std::function< double( const double ) > getMassFunction( )
    {
        return massFunction_;
    }

    std::function< Eigen::Vector3d( const double ) > getCenterOfMassFunction( )
    {
        return centerOfMassFunction_;
    }

    std::function< Eigen::Matrix3d( const double ) > getInertiaTensorFunction( )
    {
        return inertiaTensorFunction_;
    }

protected:
    std::function< double( const double ) > massFunction_;

    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;
};

class MassDependentMassDistributionSettings : public RigidBodyPropertiesSettings
{
public:
    MassDependentMassDistributionSettings( const double currentMass,
                                           const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
                                           const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction ):
        RigidBodyPropertiesSettings( mass_dependent_rigid_body_properties ), currentMass_( currentMass ),
        centerOfMassFunction_( centerOfMassFunction ), inertiaTensorFunction_( inertiaTensorFunction )
    { }

    virtual ~MassDependentMassDistributionSettings( ) { }

    double getCurrentMass( )
    {
        return currentMass_;
    }

    std::function< Eigen::Vector3d( const double ) > getCenterOfMassFunction( )
    {
        return centerOfMassFunction_;
    }

    std::function< Eigen::Matrix3d( const double ) > getInertiaTensorFunction( )
    {
        return inertiaTensorFunction_;
    }

protected:
    double currentMass_;

    std::function< Eigen::Vector3d( const double ) > centerOfMassFunction_;

    std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction_;
};

inline std::shared_ptr< RigidBodyPropertiesSettings > constantRigidBodyPropertiesSettings(
        const double mass,
        const Eigen::Vector3d& centerOfMass = Eigen::Vector3d::Constant( TUDAT_NAN ),
        const Eigen::Matrix3d& inertiaTensor = Eigen::Matrix3d::Constant( TUDAT_NAN ) )
{
    return std::make_shared< ConstantRigidBodyPropertiesSettings >( mass, centerOfMass, inertiaTensor );
}

inline std::shared_ptr< RigidBodyPropertiesSettings > fromFunctionRigidBodyPropertiesSettings(
        const std::function< double( const double ) > massFunction,
        const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction = nullptr,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction = nullptr )
{
    return std::make_shared< FromFunctionRigidBodyPropertiesSettings >( massFunction, centerOfMassFunction, inertiaTensorFunction );
}

inline std::shared_ptr< RigidBodyPropertiesSettings > massDependentMassDistributionSettings(
        const double currentMass,
        const std::function< Eigen::Vector3d( const double ) > centerOfMassFunction,
        const std::function< Eigen::Matrix3d( const double ) > inertiaTensorFunction )
{
    return std::make_shared< MassDependentMassDistributionSettings >( currentMass, centerOfMassFunction, inertiaTensorFunction );
}

std::shared_ptr< RigidBodyProperties > createRigidBodyProperties(
        const std::shared_ptr< RigidBodyPropertiesSettings > rigidBodyPropertiesSettings,
        const std::string& body,
        const SystemOfBodies& bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEGRAVITYFIELD_H
