/*    Copyright (c) 2010-2019, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/
#ifndef TUDAT_GEOMETRYUTILS_H
#define TUDAT_GEOMETRYUTILS_H

#include <map>
#include <iostream>
#include <memory>

namespace tudat
{

namespace system_models
{
class Triangle2d
{
private:
    Eigen::Vector2d vertexA_, vertexB_, vertexC_;

public:
   // default constructor
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

class ParallelProjection
{
private:
    double minCoordinateL_, maxCoordinateL_, minCoordinateM_, maxCoordinateM_;
    Triangle2d triangle2d_;
    std::vector< double > lambdas_;
    
public:
    ParallelProjection( ):  minCoordinateL_( 0.0 ), maxCoordinateL_( 0.0 ), minCoordinateM_( 0.0 ), maxCoordinateM_( 0.0 ),
                            triangle2d_( ), lambdas_( )
    { }

    ParallelProjection( const Triangle3d& triangle, const Triangle3d& triangleToBeProjected, 
                        const Eigen::Vector3d& projectionDirection_= Eigen::Vector3d::Zero( ) )
    {
        Eigen::Vector3d dummyDirection = Eigen::Vector3d::Zero( );
        Eigen::Vector3d projectionDirection = projectionDirection_;
        if ( projectionDirection_.isApprox(dummyDirection))
        {   
            projectionDirection = ( triangle.getVersorL( ) ).cross( triangle.getVersorM( ) );
        }
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
        minCoordinateL_ = *std::min_element( coordinatesL.begin( ), coordinatesL.end( ) );
        maxCoordinateL_ = *std::max_element( coordinatesL.begin( ), coordinatesL.end( ) );
        minCoordinateM_ = *std::min_element( coordinatesM.begin( ), coordinatesM.end( ) );
        maxCoordinateM_ = *std::max_element( coordinatesM.begin( ), coordinatesM.end( ) );
    }

    Triangle2d getTriangle2d( ) const
    {
        return triangle2d_;
    }
    std::vector< double > getLambdas( ) const
    {
        return lambdas_;
    }
    double getMinL( ) const
    {
        return minCoordinateL_;
    }
    double getMaxL( ) const
    {
        return maxCoordinateL_;
    }
    double getMinM( ) const
    {
        return minCoordinateM_;
    }
    double getMaxM( ) const
    {
        return maxCoordinateM_;
    }

};

template < typename T >
constexpr const T& clamp( const T& value, const T& low, const T& high )
{
    return (value < low ) ? low : ( value > high ) ? high : value;
}

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_GEOMETRYUTILS_H
