/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POLYHEDRONACCELERATIONPARTIAL_H
#define TUDAT_POLYHEDRONACCELERATIONPARTIAL_H

#include <tudat/astro/gravitation/polyhedronGravityField.h>
#include "tudat/astro/gravitation/polyhedronGravityModel.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"

#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{
namespace acceleration_partials
{

class PolyhedronGravityPartial: public AccelerationPartial
{
public:

    //! Contructor.
    /*!
     *  Constructor, requires input on the acceleration model as of which partials are to be computed.
     *  If any partials of parameters of the rotation model of the body exerting acceleration are to be calculated,
     *  RotationMatrixPartial objects must be pre-constructed and passed here as a map, with one object for each parameter
     *  wrt which a partial is to be taken.
     *  \param acceleratedBody Name of body undergoing acceleration.
     *  \param acceleratingBody Name of body exerting acceleration.
     *  \param accelerationModel Spherical harmonic gravity acceleration model from which acceleration is calculated wrt
     *  which the object being constructed is to calculate partials.
     *  \param rotationMatrixPartials Map of RotationMatrixPartial, one for each paramater representing a property of the
     *  rotation of the body exerting the acceleration wrt which an acceleration partial will be calculated.
     */
    PolyhedronGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
            observation_partials::RotationMatrixPartialNamedList( ) );

    //! Destructor
    ~PolyhedronGravityPartial( ){ }

    //! Function for updating the partial object to current state and time.
    /*!
     *  Function for updating the partial object to current state and time. Calculates the variables that are
     *  used for the calculation of multple partials, to prevent multiple calculations of same function.
     *  \param currentTime Time to which object is to be updated (note that most update functions are time-independent,
     *  since the 'current' state of the bodies is typically updated globally by the NBodyStateDerivative class).
     */
    void update( const double currentTime = TUDAT_NAN );

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
    }

    //! Function to retrieve partial of acceleration wrt the position of body undergoing acceleration, in inertial coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the position of the body undergoing the acceleration,
     * in inertial coordinates
     * \return Current partial of the acceleration wrt the position of the body undergoing the acceleration, in inertial coordinates.
     */
    Eigen::Matrix3d getCurrentPartialWrtPosition( )
    {
        return currentPartialWrtPosition_;
    }

    //! Function to retrieve partial of acceleration wrt the position of body undergoing acceleration, in body-fixed coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the position of the body undergoing the acceleration,
     * in body-fixed coordinates
     * \return Current partial of the acceleration wrt the position of the body undergoing the acceleration, in body-fixed coordinates.
     */
    Eigen::Matrix3d getCurrentBodyFixedPartialWrtPosition( )
    {
        return currentBodyFixedPartialWrtPosition_;
    }

    //! Function to retrieve partial of acceleration wrt the velocity of body undergoing acceleration, in inertial coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the velocity of the body undergoing the acceleration,
     * in inertial coordinates
     * \return Current partial of the acceleration wrt the velocity of the body undergoing the acceleration, in inertial coordinates.
     */
    Eigen::Matrix3d getCurrentPartialWrtVelocity( )
    {
        return currentPartialWrtVelocity_;
    }

private:

    //! Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration
    /*!
     *  Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration,
     *  set by update( time ) function.
     */
    Eigen::Vector3d bodyFixedPosition_;

    //! Current spherical coordinate of body undergoing acceleration
    /*!
     *  Current spherical coordinate of body undergoing acceleration, in reference frame fixed to body exerting acceleration.
     *  Order of components is radial distance (from center of body), latitude, longitude. Note that the the second entry
     *  differs from the direct output of the cartesian -> spherical coordinates, which produces a colatitude.
     */
    Eigen::Vector3d bodyFixedSphericalPosition_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  The partial wrt the position of the body exerting the acceleration is minus this value.
     *  Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration,
    //! with both acceleration and position in body-fixed frame.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  with both acceleration and position in body-fixed frame. Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentBodyFixedPartialWrtPosition_;

    //! The current partial of the acceleration wrt the velocity of the body undergoing the acceleration.
    /*!
     *  The current partial of the acceleration wrt the velocity of the body undergoing the acceleration.
     *  The partial wrt the velocity of the body exerting the acceleration is minus this value.
     * Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentPartialWrtVelocity_;

    //! Function to return the gravitational parameter used for calculating the acceleration.
    std::function< double( ) > gravitationalParameterFunction_;

    //! Function to return the volume used for calculating the acceleration.
    std::function< double( ) > volumeFunction_;

    //!  Polyhedron cache for this acceleration
    std::shared_ptr< gravitation::PolyhedronGravityCache > polyhedronCache_;

    //! Polyhedron facets dyads.
    std::vector< Eigen::MatrixXd > facetDyads_;

    //! Polyhedron edge dyads.
    std::vector< Eigen::MatrixXd > edgeDyads_;

    //! Function returning position of body undergoing acceleration.
    std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratedBody_;

    //! Function returning position of body exerting acceleration.
    std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratingBody_;

    //! Function return current rotation from inertial frame to frame fixed to body exerting acceleration.
    std::function< Eigen::Matrix3d( ) > fromBodyFixedToIntegrationFrameRotation_;

    //! Function to retrieve the current spherical harmonic acceleration.
    std::function< Eigen::Matrix< double, 3, 1 >( ) > accelerationFunction_;

    //! Function to update the acceleration to the current state and time.
    /*!
     *  Function to update the acceleration to the current state and time.
     *  Called when updating an object of this class with the update( time ) function,
     *  in case the partial is called before the acceleration model in the current iteration of the numerical integration.
     */
    std::function< void( const double ) > updateFunction_;

    //! Map of RotationMatrixPartial, one for each relevant rotation parameter
    /*!
     *  Map of RotationMatrixPartial, one for each parameter representing a property of the rotation of the
     *  body exerting the acceleration wrt which an acceleration partial will be calculated.
     *  Map is pre-created and set through the constructor.
     */
    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartials_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif //TUDAT_POLYHEDRONACCELERATIONPARTIAL_H
