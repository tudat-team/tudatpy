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

#include <map>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
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
   Triangle2d( ):
   vertexA_( Eigen::Vector2d::Zero() ),
   vertexB_( Eigen::Vector2d::Zero() ),
   vertexC_( Eigen::Vector2d::Zero() )
   { }

   Triangle2d( const Eigen::Vector2d& vertexA,
               const Eigen::Vector2d& vertexB,
               const Eigen::Vector2d& vertexC ): 
               vertexA_( vertexA ), 
               vertexB_( vertexB ), 
               vertexC_( vertexC ) 
   { }

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
   Triangle3d( ) :
   vertexA_( Eigen::Vector3d::Zero( ) ),
   vertexB_( Eigen::Vector3d::Zero( ) ),
   vertexC_( Eigen::Vector3d::Zero( ) ),
   versorL_( Eigen::Vector3d::Zero( ) ),
   versorM_( Eigen::Vector3d::Zero( ) )
   { }
   
   Triangle3d( const Eigen::Vector3d& vertexA,
               const Eigen::Vector3d& vertexB,
               const Eigen::Vector3d& vertexC ): 
               vertexA_( vertexA ),
               vertexB_( vertexB ),
               vertexC_( vertexC ) 
   {  
        // compute unit vectors l, m
        // for simplicity, l is chosen as edge AB and m as the scalar product of the normal and l
        versorL_ = ( vertexB_ - vertexA_ ).normalized( ) ;
        Eigen::Vector3d edgeAB = vertexB_ - vertexA_;
        Eigen::Vector3d edgeAC = vertexC_ - vertexA_;
        Eigen::Vector3d normal = edgeAB.cross(edgeAC).normalized( );
        versorM_ = normal.cross( versorL_ ).normalized( ); 
    }

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
 * Class used to compute the parallel projection of a 3d triangle wrt a given direction onto a given plane, identified by another 3d triangle.
 *
 * Class members are: values of the coordinates using the decomposition via the (l, m) base, the resulting 2d triangle, the parallel distances lambdas
 * of the vectices from the plane, logic vectors to check if the lambdas are actually positive or negative (hence, not close to zero due to numerical rounding).
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
    ParallelProjection( ):  minimumCoordinateL_( 0.0 ), maximumCoordinateL_( 0.0 ), minimumCoordinateM_( 0.0 ), maximumCoordinateM_( 0.0 ),
                            triangle2d_( ), lambdas_( )
    { }

    ParallelProjection( const Triangle3d& triangle, const Triangle3d& triangleToBeProjected, 
                        const Eigen::Vector3d& projectionDirection )
    {
        Eigen::Matrix3d A(3, 3);
        A(0, 0) =  triangle.getVersorL( )[0]; A(1, 0) =  triangle.getVersorL( )[1]; A(2, 0) =  triangle.getVersorL( )[2];
        A(0, 1) =  triangle.getVersorM( )[0]; A(1, 1) =  triangle.getVersorM( )[1]; A(2, 1) =  triangle.getVersorM( )[2];
        A(0, 2) =  -projectionDirection[0]; A(1, 2) =  -projectionDirection[1]; A(2, 2) =  -projectionDirection[2];

        std::vector< Eigen::Vector3d > vertices = { triangleToBeProjected.getVertexA( ), 
                                                    triangleToBeProjected.getVertexB( ),
                                                    triangleToBeProjected.getVertexC( ) };
        std::vector< Eigen::Vector2d > vertices2d;
        std::vector< double > lambdas( 3 );
        std::vector< double > coordinatesL( 3 );
        std::vector< double > coordinatesM( 3 );

        for (int i = 0; i<3; i++)
        {
            Eigen::Vector3d b =  vertices[i] - triangle.getVertexA( );
            Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
            Eigen::Vector2d vertix2d;
            vertix2d(0) = x( 0 );
            vertix2d(1) = x( 1 );
            lambdas[i] = x( 2 );
            coordinatesL[i] = x( 0 );
            coordinatesM[i] = x( 1 );
            vertices2d.push_back( vertix2d );
        }
        lambdas_ = lambdas;
        Triangle2d triangle2d( vertices2d[0], vertices2d[1], vertices2d[2] );
        triangle2d_ = triangle2d;
        minimumCoordinateL_ = *std::min_element( coordinatesL.begin( ), coordinatesL.end( ) );
        maximumCoordinateL_ = *std::max_element( coordinatesL.begin( ), coordinatesL.end( ) );
        minimumCoordinateM_ = *std::min_element( coordinatesM.begin( ), coordinatesM.end( ) );
        maximumCoordinateM_ = *std::max_element( coordinatesM.begin( ), coordinatesM.end( ) );

    }

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
template < typename T >
constexpr const T& clamp( const T& value, const T& low, const T& high )
{
    return (value < low ) ? low : ( value > high ) ? high : value;
}

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_GEOMETRYUTILS_H
