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
enum BodyShapeTypes
{
    spherical,
    spherical_spice,
    oblate_spheroid,
    polyhedron
};

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
    BodyShapeSettings( BodyShapeTypes bodyShapeType ):bodyShapeType_( bodyShapeType ){ }

    //  Virtual destructor
    virtual ~BodyShapeSettings( ){ }

    //  Function to return the type of body shape model that is to be created.
    /* 
     *  Function to return the type of body shape model that is to be created.
     *  \return Type of body shape model that is to be created.
     */
    BodyShapeTypes getBodyShapeType( ){ return bodyShapeType_; }

protected:

    //  Type of body shape model that is to be created.
    BodyShapeTypes bodyShapeType_;
};

//  BodyShapeSettings derived class for defining settings of a spherical shape model
//! @get_docstring(SphericalBodyShapeSettings.__docstring__)
class SphericalBodyShapeSettings: public BodyShapeSettings
{
public:

    //  Constructor
    /* 
     * Constructor
     * \param radius Radius of spherical shape model.
     */
    SphericalBodyShapeSettings( const double radius ) :
         BodyShapeSettings( spherical ), radius_( radius ){ }

    //  Function to return the radius of spherical shape model.
    /* 
     *  Function to return the radius of spherical shape model.
     *  \return Radius of spherical shape model.
     */
    double getRadius( ){ return radius_; }

    void resetRadius( const double radius ){ radius_ = radius; }

private:

    //  Radius of spherical shape model.
    double radius_;
};

//  BodyShapeSettings derived class for defining settings of an oblate spheroid (flattened sphere)
//  shape model
class OblateSphericalBodyShapeSettings: public BodyShapeSettings
{
public:

    //  Constructor
    /* 
     * Constructor
     * \param equatorialRadius Equatorial radius of spheroid shape model.
     * \param flattening Flattening of spheroid shape model.
     */
    OblateSphericalBodyShapeSettings( const double equatorialRadius,
                                      const double flattening ):
        BodyShapeSettings( oblate_spheroid ), equatorialRadius_( equatorialRadius ),
        flattening_( flattening ){ }


    //  Function to return the equatorial radius of spheroid shape model.
    /* 
     *  Function to return the equatorial radius of spheroid shape model.
     *  \return Flattening of spheroid shape model.
     */
    double getEquatorialRadius( ){ return equatorialRadius_; }

    void resetEquatorialRadius( const double equatorialRadius ){ equatorialRadius_ = equatorialRadius; }

    //  Function to return the flattening of spheroid shape model.
    /* 
     *  Function to return the flattening of spheroid shape model.
     *  \return Flattening of spheroid shape model.
     */
    double getFlattening( ){ return flattening_; }

    void resetFlattening( const double flattening ){ flattening_ = flattening; }

private:

    //  Equatorial radius of spheroid shape model.
    double equatorialRadius_;

    //  Flattening of spheroid shape model.
    double flattening_;
};

class PolyhedronBodyShapeSettings: public BodyShapeSettings
{
public:

    /*! Constructor.
     *
     * Constructor.
     * @param verticesCoordinates Matrix with coordinates of the polyhedron vertices. Each row represents the (x,y,z)
     * coordinates of one vertex.
     * @param verticesDefiningEachFacet Matrix with the indices (0 indexed) of the vertices defining each facet. Each
     * row contains 3 indices, which must be provided in counterclockwise order when seen from outside the polyhedron.
     * @param computeAltitudeWithSign Flag indicating whether the altitude should be computed with sign (i.e. >0 if
     * above surface, <0 otherwise) or having always a positive value.
     * @param justComputeDistanceToVertices Flag indicating whether the distance should be computed wrt to all the
     * polyhedron features or wrt to just the vertices.
     */
    PolyhedronBodyShapeSettings(
            const Eigen::MatrixXd& verticesCoordinates,
            const Eigen::MatrixXi& verticesDefiningEachFacet,
            const bool computeAltitudeWithSign = true,
            const bool justComputeDistanceToVertices = false):
        BodyShapeSettings( polyhedron ),
        verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        computeAltitudeWithSign_( computeAltitudeWithSign ),
        justComputeDistanceToVertices_( justComputeDistanceToVertices )
    { }

    /*! Constructor.
     *
     * Constructor.
     * @param verticesCoordinates Matrix with coordinates of the polyhedron vertices. Each row represents the (x,y,z)
     * coordinates of one vertex.
     */
    PolyhedronBodyShapeSettings(
            const Eigen::MatrixXd& verticesCoordinates):
        BodyShapeSettings( polyhedron ),
        verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( Eigen::Matrix3d::Constant( TUDAT_NAN ) ),
        computeAltitudeWithSign_( false ),
        justComputeDistanceToVertices_( true )
    { }

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

//  Function to create a body shape model.
/* 
 *  Function to create a body shape model based on model-specific settings for the shape.
 *  \param shapeSettings Settings for the shape model that is to be created, defined
 *  a pointer to an object of class (derived from) BodyShapeSettings.
 *  \param body Name of the body for which the shape model is to be created.
 *  \return Shape model created according to settings in shapeSettings.
 */
std::shared_ptr< basic_astrodynamics::BodyShapeModel > createBodyShapeModel(
        const std::shared_ptr< BodyShapeSettings > shapeSettings,
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
inline std::shared_ptr< BodyShapeSettings > oblateSphericalBodyShapeSettings( const double equatorialRadius,
																			  const double flattening )
{
	return std::make_shared< OblateSphericalBodyShapeSettings >( equatorialRadius, flattening );
}

inline std::shared_ptr< BodyShapeSettings > fullPolyhedronBodyShapeSettings(
        const Eigen::MatrixXd& verticesCoordinates,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const bool computeAltitudeWithSign = true,
        const bool justComputeDistanceToVertices = false )
{
    return std::make_shared< PolyhedronBodyShapeSettings >( verticesCoordinates, verticesDefiningEachFacet,
                                                            computeAltitudeWithSign, justComputeDistanceToVertices);
}

inline std::shared_ptr< BodyShapeSettings > simplifiedPolyhedronBodyShapeSettings(
        const Eigen::MatrixXd& verticesCoordinates )
{
    return std::make_shared< PolyhedronBodyShapeSettings >( verticesCoordinates );
}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODYSHAPEMODEL_H
