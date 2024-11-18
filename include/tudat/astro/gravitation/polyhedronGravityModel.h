/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_POLYHEDRONGRAVITYMODEL_H
#define TUDAT_POLYHEDRONGRAVITYMODEL_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include "tudat/math/basic/polyhedron.h"

namespace tudat
{
namespace gravitation
{

class PolyhedronGravitationalAccelerationModel: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{

protected:
    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector3d& ) > StateFunction;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of polyhedron paramers.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and volume of the
     * body exerting the acceleration, polyhedron parameters, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter and volume. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param aVolume A (constant) volume [m^3].
     * \param aVerticesCoordinatesMatrix A (constant) vertices coordinates matrix [m].
     * \param aVerticesDefiningEachFacetMatrix A (constant) vertices defining each facet matrix.
     * \param aVerticesDefiningEachEdgeMatrix A (constant) vertices defining each edge matrix.
     * \param aFacetDyadsVector A (constant) facet dyads vector.
     * \param aEdgeDyadsVector A (constant) edge dyads vector.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param polyhedronCache Cache object for computing/retrieving repeated terms in polyhedron potential
     *          gradient calculation.
     * \param updatePotential Flag indicating whether to update the gravitational potential when calling
     * the updateMembers function.
     * \param updateLaplacianOfPotential Flag indicating whether to update the laplacian of the
     * gravitational potential when calling the updateMembers function.
     */
    PolyhedronGravitationalAccelerationModel (
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double aVolume,
            const Eigen::MatrixXd& aVerticesCoordinatesMatrix,
            const Eigen::MatrixXi& aVerticesDefiningEachFacetMatrix,
            const Eigen::MatrixXi& aVerticesDefiningEachEdgeMatrix,
            const std::vector< Eigen::MatrixXd >& aFacetDyadsVector,
            const std::vector< Eigen::MatrixXd >& aEdgeDyadsVector,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                    [ ]( Eigen::Vector3d& input) { input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond ( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    [ ] ( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            const bool updateGravitationalPotential = false,
            const bool updateLaplacianOfGravitationalPotential = false)
        : subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction_( [ = ]( ){ return aGravitationalParameter; } ),
          volumeFunction_( [ = ]( ){ return aVolume; } ),
          getVerticesCoordinates_( [ = ]( ){ return aVerticesCoordinatesMatrix; } ),
          getVerticesDefiningEachFacet_( [ = ]( ){ return aVerticesDefiningEachFacetMatrix; } ),
          getVerticesDefiningEachEdge_( [ = ]( ){ return aVerticesDefiningEachEdgeMatrix; } ),
          getFacetDyads_( [ = ]( ){ return aFacetDyadsVector; } ),
          getEdgeDyads_( [ = ]( ){ return aEdgeDyadsVector; } ),
          sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          polyhedronCache_( std::make_shared< PolyhedronGravityCache >(
                 aVerticesCoordinatesMatrix, aVerticesDefiningEachFacetMatrix, aVerticesDefiningEachEdgeMatrix) ),
          currentPotential_( TUDAT_NAN ),
          currentLaplacianOfPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential ),
          updateLaplacianOfPotential_( updateLaplacianOfGravitationalPotential )
    { }

    //! Constructor taking functions for position of bodies, and parameters of polyhedron.
    /*!
     * Constructor taking pointer to functions returning the position of the body subject to
     * gravitational acceleration, the gravitational parameter and volume of the body exerting the
     * acceleration (central body), polyhedron parameters, and the position of the central body.
     * The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     *
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param polyhedronCache Cache object for computing/retrieving repeated terms in polyhedron potential
     *          gradient calculation.
     * \param updatePotential Flag indicating whether to update the gravitational potential when calling
     * the updateMembers function.
     * \param updateLaplacianOfPotential Flag indicating whether to update the laplacian of the
     * gravitational potential when calling the updateMembers function.
     */
    PolyhedronGravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const std::function< double() > gravitationalParameterFunction,
            const std::function< double() > volumeFunction,
            const std::function< Eigen::MatrixXd() > verticesCoordinatesFunction,
            const std::function< Eigen::MatrixXi() > verticesDefiningEachFacetFunction,
            const std::function< Eigen::MatrixXi() > verticesDefiningEachEdgeFunction,
            const std::function< std::vector< Eigen::MatrixXd >() > facetDyadsFunction,
            const std::function< std::vector< Eigen::MatrixXd >() > edgeDyadsFunction,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                [ ]( Eigen::Vector3d& input ){ input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            const bool updateGravitationalPotential = false,
            const bool updateLaplacianOfGravitationalPotential = false)
        : subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction_( gravitationalParameterFunction ),
          volumeFunction_( volumeFunction ),
          getVerticesCoordinates_( verticesCoordinatesFunction ),
          getVerticesDefiningEachFacet_( verticesDefiningEachFacetFunction ),
          getVerticesDefiningEachEdge_( verticesDefiningEachEdgeFunction ),
          getFacetDyads_( facetDyadsFunction ),
          getEdgeDyads_( edgeDyadsFunction ),
          sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
          rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          polyhedronCache_( std::make_shared< PolyhedronGravityCache >(
                 verticesCoordinatesFunction(), verticesDefiningEachFacetFunction(),
                 verticesDefiningEachEdgeFunction() ) ),
          currentPotential_( TUDAT_NAN ),
          currentLaplacianOfPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential ),
          updateLaplacianOfPotential_( updateLaplacianOfGravitationalPotential )
    { }

    ~PolyhedronGravitationalAccelerationModel( ){ }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class members of this class.
     * The potential and laplacian of potential are only updated if the associated flags indicate so.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {

            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );

            subjectPositionFunction_( positionOfBodySubjectToAcceleration_ );
            sourcePositionFunction_( positionOfBodyExertingAcceleration_ );
            currentInertialRelativePosition_ =
                positionOfBodySubjectToAcceleration_ - positionOfBodyExertingAcceleration_ ;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;

            polyhedronCache_->update( currentRelativePosition_ );

            // Compute the current acceleration
            currentAccelerationInBodyFixedFrame_ = basic_mathematics::calculatePolyhedronGradientOfGravitationalPotential(
                gravitationalParameterFunction_() / volumeFunction_(),
                polyhedronCache_->getVerticesCoordinatesRelativeToFieldPoint(),
                getVerticesDefiningEachFacet_(),
                getVerticesDefiningEachEdge_(),
                getFacetDyads_(),
                getEdgeDyads_(),
                polyhedronCache_->getPerFacetFactor(),
                polyhedronCache_->getPerEdgeFactor() );

            currentAcceleration_ = rotationToIntegrationFrame_ * currentAccelerationInBodyFixedFrame_;

            // Compute the current gravitational potential
            if ( updatePotential_ )
            {
                currentPotential_ = basic_mathematics::calculatePolyhedronGravitationalPotential(
                    gravitationalParameterFunction_( ) / volumeFunction_( ),
                    polyhedronCache_->getVerticesCoordinatesRelativeToFieldPoint( ),
                    getVerticesDefiningEachFacet_( ),
                    getVerticesDefiningEachEdge_( ),
                    getFacetDyads_( ),
                    getEdgeDyads_( ),
                    polyhedronCache_->getPerFacetFactor( ),
                    polyhedronCache_->getPerEdgeFactor( ) );
            }

            // Compute the current laplacian
            if ( updateLaplacianOfPotential_ )
            {
                currentLaplacianOfPotential_ = basic_mathematics::calculatePolyhedronLaplacianOfGravitationalPotential(
                    gravitationalParameterFunction_( ) / volumeFunction_( ),
                    polyhedronCache_->getPerFacetFactor( ) );
            }
        }
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
    //! fixed to body undergoing acceleration
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
    //! frame
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

    //! Function to return the function returning the relevant gravitational parameter.
    std::function< double( ) > getGravitationalParameterFunction( )
    {
        return gravitationalParameterFunction_;
    }

    //! Function to return the function returning the relevant volume.
    std::function< double( ) > getVolumeFunction( )
    {
        return volumeFunction_;
    }

    //! Function to return the function returning the facet dyads.
    std::function< std::vector< Eigen::MatrixXd >( ) > getFacetDyadsFunction( )
    {
        return getFacetDyads_;
    }

    //! Function to return the function returning the facet dyads.
    std::function< std::vector< Eigen::MatrixXd >( ) > getEdgeDyadsFunction( )
    {
        return getEdgeDyads_;
    }

    //! Function to return current position vector of body exerting gravitational acceleration in inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodySubjectToAcceleration( )
    {
        return positionOfBodySubjectToAcceleration_;
    }

    //! Function to return current position vector of body undergoing gravitational acceleration in inertial frame.
    Eigen::Vector3d getCurrentPositionOfBodyExertingAcceleration( )
    {
        return positionOfBodyExertingAcceleration_;
    }

    //! Function to return the function returning position of body exerting acceleration.
    /*!
     * Function to return the function returning position of body exerting acceleration.
     * \return Function returning position of body exerting acceleration.
     */
    StateFunction getStateFunctionOfBodyExertingAcceleration( )
    { return sourcePositionFunction_; }

    //! Function to return the function returning position of body subject to acceleration.
    /*!
     * Function to return the function returning position of body subject to acceleration.
     * \return Function returning position of body subject to acceleration.
     */
    StateFunction getStateFunctionOfBodyUndergoingAcceleration( )
    { return subjectPositionFunction_; }

    //! Function to retrieve the spherical harmonics cache for this acceleration.
    std::shared_ptr< PolyhedronGravityCache > getPolyhedronCache( )
    {
        return polyhedronCache_;
    }

    //! Function to return the value of the current gravitational potential.
    double getCurrentPotential ( )
    { return currentPotential_; }

    //! Function to return the value of the current laplacian of the gravitational potential.
    double getCurrentLaplacianOfPotential ( )
    { return currentLaplacianOfPotential_; }

    //! Function to return the update potential flag.
    bool getUpdatePotential ( )
    { return updatePotential_; }

    //! Function to reset the update potential flag.
    void resetUpdatePotential ( bool updatePotential )
    { updatePotential_ = updatePotential; }

    //! Function to return the update laplacian of potential flag.
    bool getUpdateLaplacianOfPotential ( )
    { return updateLaplacianOfPotential_; }

    //! Function to reset the update laplacian of potential flag.
    void resetUpdateLaplacianOfPotential ( bool updateLaplacianOfPotential )
    { updateLaplacianOfPotential_ = updateLaplacianOfPotential; }

private:

    //! Pointer to function returning position of body subject to acceleration.
    const StateFunction subjectPositionFunction_;

    //! Function returning a gravitational parameter [m^3 s^-2].
    const std::function< double( ) > gravitationalParameterFunction_;

    //! Function returning a volume [m^3].
    const std::function< double( ) > volumeFunction_;

    //! Pointer to function returning the vertices coordinates
    const std::function< Eigen::MatrixXd( ) > getVerticesCoordinates_;

    //! Pointer to function returning the vertices defining each facet
    const std::function< Eigen::MatrixXi( ) > getVerticesDefiningEachFacet_;

    //! Pointer to function returning the vertices defining each edge
    const std::function< Eigen::MatrixXi( ) > getVerticesDefiningEachEdge_;

    //! Pointer to function returning the facet dyads
    const std::function< std::vector< Eigen::MatrixXd >( ) > getFacetDyads_;

    //! Pointer to function returning the edge dyads
    const std::function< std::vector< Eigen::MatrixXd >( ) > getEdgeDyads_;

    //! Pointer to function returning position of body exerting acceleration.
    const StateFunction sourcePositionFunction_;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Variable denoting whether mutual acceleration between bodies is included.
    bool isMutualAttractionUsed_;

    //!  Polyhedron cache for this acceleration
    std::shared_ptr< PolyhedronGravityCache > polyhedronCache_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    //! Current acceleration in frame fixed to body undergoing acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    //! Position of body subject to acceleration.
    Eigen::Vector3d positionOfBodySubjectToAcceleration_;

    //! Position of body exerting acceleration.
    Eigen::Vector3d positionOfBodyExertingAcceleration_;

    //! Current gravitational potential acting on the body undergoing acceleration, as computed by last call to
    //! updateMembers function
    double currentPotential_;

    //! Current laplacian of the gravitational potential acting on the body undergoing acceleration, as computed by last
    //! call to updateMembers function
    double currentLaplacianOfPotential_;

    //! Flag indicating whether to update the gravitational potential when calling the updateMembers function.
    bool updatePotential_;

    //! Flag indicating whether to update the laplacian of the gravitational potential when calling the updateMembers function.
    bool updateLaplacianOfPotential_;
};


} // namespace gravitation

} // namespace tudat

#endif //TUDAT_POLYHEDRONGRAVITYMODEL_H
