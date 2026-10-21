/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_PANELGEOMETRYUTILS_H
#define TUDAT_PANELGEOMETRYUTILS_H

#include <Eigen/Core>
#include <vector>

namespace tudat
{
namespace system_models
{
/*!
 * Class containing all the information regarding a two-dimendional triangle.
 *
 * Class members are: 2d points (vertices of the triangle).
 */
class Triangle2d
{
private:
    Eigen::Vector2d vertexA_, vertexB_, vertexC_;

public:
    Triangle2d( );

    Triangle2d( const Eigen::Vector2d& vertexA, const Eigen::Vector2d& vertexB, const Eigen::Vector2d& vertexC );

    Eigen::Vector2d getVertexA( ) const
    {
        return vertexA_;
    }
    Eigen::Vector2d getVertexB( ) const
    {
        return vertexB_;
    }
    Eigen::Vector2d getVertexC( ) const
    {
        return vertexC_;
    }
};
/*!
 * Class containing all the information regarding a three-dimendional triangle.
 *
 * Class members are: 3d points (vertices of the triangle) and two unit vectors (base of the plane on which the triangle lies).
 */
class Triangle3d
{
private:
    Eigen::Vector3d vertexA_, vertexB_, vertexC_;
    Eigen::Vector3d versorL_, versorM_;

public:
    Triangle3d( );

    Triangle3d( const Eigen::Vector3d& vertexA, const Eigen::Vector3d& vertexB, const Eigen::Vector3d& vertexC );

    Eigen::Vector3d getVertexA( ) const
    {
        return vertexA_;
    }
    Eigen::Vector3d getVertexB( ) const
    {
        return vertexB_;
    }
    Eigen::Vector3d getVertexC( ) const
    {
        return vertexC_;
    }
    Eigen::Vector3d getVersorL( ) const
    {
        return versorL_;
    }
    Eigen::Vector3d getVersorM( ) const
    {
        return versorM_;
    }
};
/*!
 * Class used to compute the parallel projection of a 3d triangle wrt a given direction onto a given plane, identified by another 3d
 * triangle.
 *
 * Class members are: values of the coordinates using the decomposition via the (l, m) base, the resulting 2d triangle, the parallel
 * distances lambdas of the vectices from the plane, logic vectors to check if the lambdas are actually positive or negative (hence, not
 * close to zero due to numerical rounding).
 */
class ParallelProjection
{
private:
    double minimumCoordinateL_, maximumCoordinateL_, minimumCoordinateM_, maximumCoordinateM_;
    Triangle2d triangle2d_;
    std::vector< double > lambdas_;
    std::vector< bool > areLambdasActuallyPositive_;
    std::vector< bool > areLambdasActuallyNegative_;

public:
    ParallelProjection( );

    ParallelProjection( const Triangle3d& triangle, const Triangle3d& triangleToBeProjected, const Eigen::Vector3d& projectionDirection );

    Triangle2d getTriangle2d( ) const
    {
        return triangle2d_;
    }
    std::vector< double > getLambdas( ) const
    {
        return lambdas_;
    }
    double getMinimumL( ) const
    {
        return minimumCoordinateL_;
    }
    double getMaximumL( ) const
    {
        return maximumCoordinateL_;
    }
    double getMinimumM( ) const
    {
        return minimumCoordinateM_;
    }
    double getMaximumM( ) const
    {
        return maximumCoordinateM_;
    }

    void setAreLambdasActuallyPositive( const std::vector< bool > areLambdasActuallyPositive )
    {
        areLambdasActuallyPositive_ = areLambdasActuallyPositive;
    }

    std::vector< bool > getAreLambdasActuallyPositive( ) const
    {
        return areLambdasActuallyPositive_;
    }

    void setAreLambdasActuallyNegative( const std::vector< bool > areLambdasActuallyNegative )
    {
        areLambdasActuallyNegative_ = areLambdasActuallyNegative;
    }

    std::vector< bool > getAreLambdasActuallyNegative( ) const
    {
        return areLambdasActuallyNegative_;
    }
};
/*!
 * Function to clamp a value given two the boundary values.
 */
template< typename T >
constexpr const T& clamp( const T& value, const T& low, const T& high )
{
    return ( value < low ) ? low : ( value > high ) ? high : value;
}

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_GEOMETRYUTILS_H
