/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODYSHAPEMODEL_H
#define TUDAT_CREATEBODYSHAPEMODEL_H

#include <memory>

#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

//  Types of body shape models that can be used.
enum BodyShapeTypes { spherical, spherical_spice, oblate_spheroid, oblate_spice, polyhedron_shape, hybrid_shape };

//  Class for providing settings for body shape model.
/*
 *  Class for providing settings for automatic body shape model creation. This class is a functional
 *  (base) class for settings of body shapels models that require no information in addition to
 *  their type. Types requiring additional information must be created using an object derived from
 *  this class.
 */

//! @get_docstring(BodyShapeSettings.__docstring__)
class BodyShapeSettings
{
public:
    //  Constructor
    /*
     *Constructor
     * \param bodyShapeType Type of body shape model that is to be created.
     */
    BodyShapeSettings( BodyShapeTypes bodyShapeType ): bodyShapeType_( bodyShapeType ) { }

    //  Virtual destructor
    virtual ~BodyShapeSettings( ) { }

    //  Function to return the type of body shape model that is to be created.
    /*
     *  Function to return the type of body shape model that is to be created.
     *  \return Type of body shape model that is to be created.
     */
    BodyShapeTypes getBodyShapeType( )
    {
        return bodyShapeType_;
    }

protected:
    //  Type of body shape model that is to be created.
    BodyShapeTypes bodyShapeType_;
};

//  BodyShapeSettings derived class for defining settings of a spherical shape model
//! @get_docstring(SphericalBodyShapeSettings.__docstring__)
class SphericalBodyShapeSettings : public BodyShapeSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param radius Radius of spherical shape model.
     */
    SphericalBodyShapeSettings( const double radius ): BodyShapeSettings( spherical ), radius_( radius ) { }

    //  Function to return the radius of spherical shape model.
    /*
     *  Function to return the radius of spherical shape model.
     *  \return Radius of spherical shape model.
     */
    double getRadius( )
    {
        return radius_;
    }

    void resetRadius( const double radius )
    {
        radius_ = radius;
    }

private:
    //  Radius of spherical shape model.
    double radius_;
};

//  BodyShapeSettings derived class for defining settings of an oblate spheroid (flattened sphere)
//  shape model
class OblateSphericalBodyShapeSettings : public BodyShapeSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param equatorialRadius Equatorial radius of spheroid shape model.
     * \param flattening Flattening of spheroid shape model.
     */
    OblateSphericalBodyShapeSettings( const double equatorialRadius, const double flattening ):
        BodyShapeSettings( oblate_spheroid ), equatorialRadius_( equatorialRadius ), flattening_( flattening )
    { }

    //  Function to return the equatorial radius of spheroid shape model.
    /*
     *  Function to return the equatorial radius of spheroid shape model.
     *  \return Flattening of spheroid shape model.
     */
    double getEquatorialRadius( )
    {
        return equatorialRadius_;
    }

    void resetEquatorialRadius( const double equatorialRadius )
    {
        equatorialRadius_ = equatorialRadius;
    }

    //  Function to return the flattening of spheroid shape model.
    /*
     *  Function to return the flattening of spheroid shape model.
     *  \return Flattening of spheroid shape model.
     */
    double getFlattening( )
    {
        return flattening_;
    }

    void resetFlattening( const double flattening )
    {
        flattening_ = flattening;
    }

private:
    //  Equatorial radius of spheroid shape model.
    double equatorialRadius_;

    //  Flattening of spheroid shape model.
    double flattening_;
};

//  BodyShapeSettings derived class for defining settings of a polyhedron shape model
class PolyhedronBodyShapeSettings : public BodyShapeSettings
{
public:
    /*! Constructor.
     *
     * Constructor. Defaults define a shape model for which the distance is computed wrt to all polyhedron
     * features (vertices, facets and edges) and for which the altitude is signed.
     *
     * The computeAltitudeWithSign indicates whether the altitude should be computed with sign (i.e. >0 if above surface,
     * <0 otherwise) or having always a positive value. If the the sign of the altitude is not relevant, then setting
     * it to "false" is recommended, as that reduces the CPU time for the computation of the altitude.
     *
     * Setting the justComputeDistanceToVertices flag to true defines a polyhedron for which the distance is computed
     * solely wrt the vertices. Useful for medium to high altitudes, as it allows reducing the CPU time with only small
     * errors in the computed altitude ("small" depends on the resolution of the used polyhedron and altitude).
     * To further reduce the CPU time, this simplified computation of the altitude might be applied to a low-resolution
     * polyhedron. For further discussion see Avillez (2022), MSc thesis (TU Delft).
     *
     * @param verticesCoordinates Matrix with coordinates of the polyhedron vertices. Each row represents the (x,y,z)
     * coordinates of one vertex.
     * @param verticesDefiningEachFacet Matrix with the indices (0 indexed) of the vertices defining each facet. Each
     * row contains 3 indices, which must be provided in counterclockwise order when seen from outside the polyhedron.
     * @param computeAltitudeWithSign Flag indicating whether the altitude should be computed with sign (i.e. >0 if
     * above surface, <0 otherwise) or having always a positive value.
     * @param justComputeDistanceToVertices Flag indicating whether the distance should be computed wrt to all the
     * polyhedron features or wrt to just the vertices.
     */
    PolyhedronBodyShapeSettings( const Eigen::MatrixXd& verticesCoordinates,
                                 const Eigen::MatrixXi& verticesDefiningEachFacet,
                                 const bool computeAltitudeWithSign = true,
                                 const bool justComputeDistanceToVertices = false ):
        BodyShapeSettings( polyhedron_shape ), verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ), computeAltitudeWithSign_( computeAltitudeWithSign ),
        justComputeDistanceToVertices_( justComputeDistanceToVertices )
    { }

    //! Destructor
    ~PolyhedronBodyShapeSettings( ) { }

    // Function to return the vertices coordinates.
    const Eigen::MatrixXd& getVerticesCoordinates( )
    {
        return verticesCoordinates_;
    }

    // Function to reset the vertices coordinates.
    void resetVerticesCoordinates( const Eigen::MatrixXd& verticesCoordinates )
    {
        verticesCoordinates_ = verticesCoordinates;
    }

    // Function to return the vertices defining each facet.
    const Eigen::MatrixXi& getVerticesDefiningEachFacet( )
    {
        return verticesDefiningEachFacet_;
    }

    // Function to reset the vertices defining each facet.
    void resetVerticesDefiningEachFacet( const Eigen::MatrixXi& verticesDefiningEachFacet )
    {
        verticesDefiningEachFacet_ = verticesDefiningEachFacet;
    }

    // Function to return the computeAltitudeWithSign flag.
    bool getComputeAltitudeWithSign( )
    {
        return computeAltitudeWithSign_;
    }

    // Function to reset the computeAltitudeWithSign flag.
    void resetComputeAltitudeWithSign( bool computeAltitudeWithSign )
    {
        computeAltitudeWithSign_ = computeAltitudeWithSign;
    }

    // Function to return the justComputeDistanceToVertices flag.
    bool getJustComputeDistanceToVertices( )
    {
        return justComputeDistanceToVertices_;
    }

    // Function to reset the justComputeDistanceToVertices flag.
    void resetJustComputeDistanceToVertices( bool justComputeDistanceToVertices )
    {
        justComputeDistanceToVertices_ = justComputeDistanceToVertices;
    }

private:
    // Matrix with coordinates of the polyhedron vertices.
    Eigen::MatrixXd verticesCoordinates_;

    // Matrix with the indices (0 indexed) of the vertices defining each facet.
    Eigen::MatrixXi verticesDefiningEachFacet_;

    // Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
    // having always a positive value
    bool computeAltitudeWithSign_;

    // Flag indicating whether the distance should be computed wrt to all the polyhedron features or wrt to just the
    // vertices.
    bool justComputeDistanceToVertices_;
};

//  BodyShapeSettings derived class for defining settings of a hybrid shape model.
class HybridBodyShapeSettings : public BodyShapeSettings
{
public:
    /*! Constructor.
     *
     * Constructor. Settings to define a hybrid model, consisting of a low-resolution model which is used at high altitudes
     * (above the switchover altitude) and a high-resolution model used at low altitudes (below the switchover altitude).
     * Useful when the evaluation of the high-resolution model is computationally expensive (e.g. polyhedron model).
     * @param lowResolutionBodyShapeSettings Body shape settings of the low-resolution model.
     * @param highResolutionBodyShapeSettings Body shape settings of the high-resolution model.
     * @param switchoverAltitude Altitude at which the model used to compute the altitude is changed.
     */
    HybridBodyShapeSettings( const std::shared_ptr< BodyShapeSettings > lowResolutionBodyShapeSettings,
                             const std::shared_ptr< BodyShapeSettings > highResolutionBodyShapeSettings,
                             double switchoverAltitude ):
        BodyShapeSettings( hybrid_shape ), lowResolutionBodyShapeSettings_( lowResolutionBodyShapeSettings ),
        highResolutionBodyShapeSettings_( highResolutionBodyShapeSettings ), switchoverAltitude_( switchoverAltitude )
    { }

    //! Destructor
    ~HybridBodyShapeSettings( ) { }

    // Function to return the body shape settings of the low-resolution model.
    std::shared_ptr< BodyShapeSettings > getLowResolutionBodyShapeSettings( )
    {
        return lowResolutionBodyShapeSettings_;
    }

    // Function to reset the body shape settings of the low-resolution model.
    void resetLowResolutionBodyShapeSettings( std::shared_ptr< BodyShapeSettings > lowResolutionBodyShapeSettings )
    {
        lowResolutionBodyShapeSettings_ = lowResolutionBodyShapeSettings;
    }

    // Function to return the body shape settings of the high-resolution model.
    std::shared_ptr< BodyShapeSettings > getHighResolutionBodyShapeSettings( )
    {
        return highResolutionBodyShapeSettings_;
    }

    // Function to reset the body shape settings of the high-resolution model.
    void resetHighResolutionBodyShapeSettings( std::shared_ptr< BodyShapeSettings > highResolutionBodyShapeSettings )
    {
        highResolutionBodyShapeSettings_ = highResolutionBodyShapeSettings;
    }

    // Function to return the switchover altitude
    double getSwitchoverAltitude( )
    {
        return switchoverAltitude_;
    }

    // Function to reset the switchover altitude
    void resetSwitchoverAltitude( double switchoverAltitude )
    {
        switchoverAltitude_ = switchoverAltitude;
    }

private:
    // Settings of the model used to compute the altitude for altitudes higher than switchoverAltitude.
    std::shared_ptr< BodyShapeSettings > lowResolutionBodyShapeSettings_;

    // Settings of the model used to compute the altitude for altitudes lower than switchoverAltitude.
    std::shared_ptr< BodyShapeSettings > highResolutionBodyShapeSettings_;

    // Altitude at which the model used to compute the altitude is changed.
    double switchoverAltitude_;
};

//  Function to create a body shape model.
/*
 *  Function to create a body shape model based on model-specific settings for the shape.
 *  \param shapeSettings Settings for the shape model that is to be created, defined
 *  a pointer to an object of class (derived from) BodyShapeSettings.
 *  \param body Name of the body for which the shape model is to be created.
 *  \return Shape model created according to settings in shapeSettings.
 */
std::shared_ptr< basic_astrodynamics::BodyShapeModel > createBodyShapeModel( const std::shared_ptr< BodyShapeSettings > shapeSettings,
                                                                             const std::string& body );

//! @get_docstring(sphericalBodyShapeSettings)
inline std::shared_ptr< BodyShapeSettings > sphericalBodyShapeSettings( const double radius )
{
    return std::make_shared< SphericalBodyShapeSettings >( radius );
}

//! @get_docstring(fromSpiceSphericalBodyShapeSettings)
inline std::shared_ptr< BodyShapeSettings > fromSpiceSphericalBodyShapeSettings( )
{
    return std::make_shared< BodyShapeSettings >( BodyShapeTypes::spherical_spice );
}

//! @get_docstring(oblateSphericalBodyShapeSettings)
inline std::shared_ptr< BodyShapeSettings > oblateSphericalBodyShapeSettings( const double equatorialRadius, const double flattening )
{
    return std::make_shared< OblateSphericalBodyShapeSettings >( equatorialRadius, flattening );
}

inline std::shared_ptr< BodyShapeSettings > fromSpiceOblateSphericalBodyShapeSettings( )
{
    return std::make_shared< BodyShapeSettings >( oblate_spice );
}

inline std::shared_ptr< BodyShapeSettings > polyhedronBodyShapeSettings( const Eigen::MatrixXd& verticesCoordinates,
                                                                         const Eigen::MatrixXi& verticesDefiningEachFacet,
                                                                         const bool computeAltitudeWithSign = true,
                                                                         const bool justComputeDistanceToVertices = false )
{
    return std::make_shared< PolyhedronBodyShapeSettings >(
            verticesCoordinates, verticesDefiningEachFacet, computeAltitudeWithSign, justComputeDistanceToVertices );
}

inline std::shared_ptr< BodyShapeSettings > hybridBodyShapeSettings( std::shared_ptr< BodyShapeSettings > lowResolutionBodyShapeSettings,
                                                                     std::shared_ptr< BodyShapeSettings > highResolutionBodyShapeSettings,
                                                                     double switchoverAltitude )
{
    return std::make_shared< HybridBodyShapeSettings >(
            lowResolutionBodyShapeSettings, highResolutionBodyShapeSettings, switchoverAltitude );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEBODYSHAPEMODEL_H
