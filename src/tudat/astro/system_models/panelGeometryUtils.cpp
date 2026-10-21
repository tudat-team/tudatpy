/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/system_models/panelGeometryUtils.h"

#include <algorithm>

#include <Eigen/Geometry>
#include <Eigen/QR>

namespace tudat
{
namespace system_models
{

Triangle2d::Triangle2d( ): vertexA_( Eigen::Vector2d::Zero( ) ), vertexB_( Eigen::Vector2d::Zero( ) ), vertexC_( Eigen::Vector2d::Zero( ) )
{}

Triangle2d::Triangle2d( const Eigen::Vector2d& vertexA, const Eigen::Vector2d& vertexB, const Eigen::Vector2d& vertexC ):
    vertexA_( vertexA ), vertexB_( vertexB ), vertexC_( vertexC )
{}

Triangle3d::Triangle3d( ):
    vertexA_( Eigen::Vector3d::Zero( ) ), vertexB_( Eigen::Vector3d::Zero( ) ), vertexC_( Eigen::Vector3d::Zero( ) ),
    versorL_( Eigen::Vector3d::Zero( ) ), versorM_( Eigen::Vector3d::Zero( ) )
{}

Triangle3d::Triangle3d( const Eigen::Vector3d& vertexA, const Eigen::Vector3d& vertexB, const Eigen::Vector3d& vertexC ):
    vertexA_( vertexA ), vertexB_( vertexB ), vertexC_( vertexC )
{
    // For simplicity, l is chosen as edge AB and m as n x l.
    versorL_ = ( vertexB_ - vertexA_ ).normalized( );
    const Eigen::Vector3d edgeAB = vertexB_ - vertexA_;
    const Eigen::Vector3d edgeAC = vertexC_ - vertexA_;
    const Eigen::Vector3d normal = edgeAB.cross( edgeAC ).normalized( );
    versorM_ = normal.cross( versorL_ ).normalized( );
}

ParallelProjection::ParallelProjection( ):
    minimumCoordinateL_( 0.0 ), maximumCoordinateL_( 0.0 ), minimumCoordinateM_( 0.0 ), maximumCoordinateM_( 0.0 ), triangle2d_( ),
    lambdas_( )
{}

ParallelProjection::ParallelProjection( const Triangle3d& triangle,
                                        const Triangle3d& triangleToBeProjected,
                                        const Eigen::Vector3d& projectionDirection )
{
    Eigen::Matrix3d A( 3, 3 );
    A( 0, 0 ) = triangle.getVersorL( )[ 0 ];
    A( 1, 0 ) = triangle.getVersorL( )[ 1 ];
    A( 2, 0 ) = triangle.getVersorL( )[ 2 ];
    A( 0, 1 ) = triangle.getVersorM( )[ 0 ];
    A( 1, 1 ) = triangle.getVersorM( )[ 1 ];
    A( 2, 1 ) = triangle.getVersorM( )[ 2 ];
    A( 0, 2 ) = -projectionDirection[ 0 ];
    A( 1, 2 ) = -projectionDirection[ 1 ];
    A( 2, 2 ) = -projectionDirection[ 2 ];

    const std::vector< Eigen::Vector3d > vertices = { triangleToBeProjected.getVertexA( ),
                                                      triangleToBeProjected.getVertexB( ),
                                                      triangleToBeProjected.getVertexC( ) };
    std::vector< Eigen::Vector2d > vertices2d;
    std::vector< double > lambdas( 3 );
    std::vector< double > coordinatesL( 3 );
    std::vector< double > coordinatesM( 3 );

    for( int i = 0; i < 3; i++ )
    {
        const Eigen::Vector3d b = vertices[ i ] - triangle.getVertexA( );
        const Eigen::Vector3d x = A.colPivHouseholderQr( ).solve( b );
        Eigen::Vector2d vertix2d;
        vertix2d( 0 ) = x( 0 );
        vertix2d( 1 ) = x( 1 );
        lambdas[ i ] = x( 2 );
        coordinatesL[ i ] = x( 0 );
        coordinatesM[ i ] = x( 1 );
        vertices2d.push_back( vertix2d );
    }
    lambdas_ = lambdas;
    triangle2d_ = Triangle2d( vertices2d[ 0 ], vertices2d[ 1 ], vertices2d[ 2 ] );
    minimumCoordinateL_ = *std::min_element( coordinatesL.begin( ), coordinatesL.end( ) );
    maximumCoordinateL_ = *std::max_element( coordinatesL.begin( ), coordinatesL.end( ) );
    minimumCoordinateM_ = *std::min_element( coordinatesM.begin( ), coordinatesM.end( ) );
    maximumCoordinateM_ = *std::max_element( coordinatesM.begin( ), coordinatesM.end( ) );
}

}  // namespace system_models

}  // namespace tudat
