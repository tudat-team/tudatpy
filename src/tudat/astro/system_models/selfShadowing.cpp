/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>
#include <iostream>
#include <algorithm>
#include <memory>

#include <thread>
#include "tudat/astro/system_models/panelGeometryUtils.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/selfShadowing.h"

namespace tudat
{

namespace system_models
{

bool firstDiscrimination( const std::shared_ptr< system_models::VehicleExteriorPanel > panel, const Eigen::Vector3d& v )
{
    if( panel->getBodyFixedSurfaceNormal( )( ).dot( v ) < 0 )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
// algorithm from "SELF-SHADOWING OF A SPACECRAFT IN THE COMPUTATION OF SURFACE FORCES. AN EXAMPLE IN PLANETARY GEODESY" (Balmino et al.,
// 2018) section 4.3 "Algorithms of shadowing of S_0 by S"
bool secondDiscrimination( ParallelProjection& projection )
{
    std::vector< double > lambdas = projection.getLambdas( );

    const double EPSILON = 1e-12;

    int nonZeroCount = 0;
    if( std::abs( lambdas.at( 0 ) ) > EPSILON ) nonZeroCount++;
    if( std::abs( lambdas.at( 1 ) ) > EPSILON ) nonZeroCount++;
    if( std::abs( lambdas.at( 2 ) ) > EPSILON ) nonZeroCount++;

    if( nonZeroCount == 0 )
    {
        std::vector< bool > dummyVectorPositive = { true, true, true };
        projection.setAreLambdasActuallyPositive( dummyVectorPositive );
        std::vector< bool > dummyVectorNegative = { false, false, false };
        projection.setAreLambdasActuallyNegative( dummyVectorNegative );
        return 1;
    }

    bool hasPositive = ( lambdas.at( 0 ) > EPSILON ) || ( lambdas.at( 1 ) > EPSILON ) || ( lambdas.at( 2 ) > EPSILON );

    std::vector< bool > areLambdasActuallyPositive = { ( lambdas.at( 0 ) > EPSILON ),
                                                       ( lambdas.at( 1 ) > EPSILON ),
                                                       ( lambdas.at( 2 ) > EPSILON ) };
    projection.setAreLambdasActuallyPositive( areLambdasActuallyPositive );

    std::vector< bool > areLambdasActuallyNegative = { ( lambdas.at( 0 ) < -EPSILON ),
                                                       ( lambdas.at( 1 ) < -EPSILON ),
                                                       ( lambdas.at( 2 ) < -EPSILON ) };
    projection.setAreLambdasActuallyNegative( areLambdasActuallyNegative );

    return hasPositive ? 1 : 0;
}
// algorithm cited by "SELF-SHADOWING OF A SPACECRAFT IN THE COMPUTATION OF SURFACE FORCES. AN EXAMPLE IN PLANETARY GEODESY" (Balmino et
// al., 2018) appendix B "Point-in-polygon (PIP) algorithms", here is a custom implementation of the ray-casting PIP version to check an
// entire triangle
std::vector< bool > isTriangleInTriangle( const ParallelProjection& projection1, const ParallelProjection& projection2 )
{
    std::vector< Eigen::Vector2d > testPoints = { projection2.getTriangle2d( ).getVertexA( ),
                                                  projection2.getTriangle2d( ).getVertexB( ),
                                                  projection2.getTriangle2d( ).getVertexC( ) };
    std::vector< bool > testPointsLambdaActuallyPositive = projection2.getAreLambdasActuallyPositive( );
    std::vector< bool > results( 3 );
    Eigen::Vector2d testPoint;
    Eigen::Vector2d pointA = projection1.getTriangle2d( ).getVertexA( );
    Eigen::Vector2d pointB = projection1.getTriangle2d( ).getVertexB( );
    Eigen::Vector2d pointC = projection1.getTriangle2d( ).getVertexC( );
    double ymin, ymax, xmax;
    double xIntersect;
    int count;
    bool resultToBeChecked;

    for( int i = 0; i < 3; i++ )
    {
        count = 0;
        testPoint = testPoints[ i ];

        // edge from pointA to pointB
        {
            ymin = ( pointA( 1 ) < pointB( 1 ) ) ? pointA( 1 ) : pointB( 1 );
            ymax = ( pointA( 1 ) < pointB( 1 ) ) ? pointB( 1 ) : pointA( 1 );
            xmax = ( pointA( 0 ) > pointB( 0 ) ) ? pointA( 0 ) : pointB( 0 );
            if( ( testPoint( 1 ) > ymin ) && ( testPoint( 1 ) <= ymax ) && ( testPoint( 0 ) <= xmax ) )
            {
                xIntersect = pointA( 0 ) + ( testPoint( 1 ) - pointA( 1 ) ) * ( pointB( 0 ) - pointA( 0 ) ) / ( pointB( 1 ) - pointA( 1 ) );
                if( pointA( 0 ) == pointB( 0 ) || testPoint( 0 ) <= xIntersect ) count++;
            }
        }
        // edge from pointB to pointC
        {
            ymin = ( pointB( 1 ) < pointC( 1 ) ) ? pointB( 1 ) : pointC( 1 );
            ymax = ( pointB( 1 ) < pointC( 1 ) ) ? pointC( 1 ) : pointB( 1 );
            xmax = ( pointB( 0 ) > pointC( 0 ) ) ? pointB( 0 ) : pointC( 0 );
            if( ymax == ymin && testPoint( 1 ) == ymax )
            {
                count++;
            }
            if( ( testPoint( 1 ) > ymin ) && ( testPoint( 1 ) <= ymax ) && ( testPoint( 0 ) <= xmax ) )
            {
                xIntersect = pointB( 0 ) + ( testPoint( 1 ) - pointB( 1 ) ) * ( pointC( 0 ) - pointB( 0 ) ) / ( pointC( 1 ) - pointB( 1 ) );
                if( pointB( 0 ) == pointC( 0 ) || testPoint( 0 ) <= xIntersect ) count++;
            }
        }
        // edge from pointC to pointA
        {
            ymin = ( pointC( 1 ) < pointA( 1 ) ) ? pointC( 1 ) : pointA( 1 );
            ymax = ( pointC( 1 ) < pointA( 1 ) ) ? pointA( 1 ) : pointC( 1 );
            xmax = ( pointC( 0 ) > pointA( 0 ) ) ? pointC( 0 ) : pointA( 0 );
            if( ( testPoint( 1 ) > ymin ) && ( testPoint( 1 ) <= ymax ) && ( testPoint( 0 ) <= xmax ) )
            {
                xIntersect = pointC( 0 ) + ( testPoint( 1 ) - pointC( 1 ) ) * ( pointA( 0 ) - pointC( 0 ) ) / ( pointA( 1 ) - pointC( 1 ) );
                if( pointC( 0 ) == pointA( 0 ) || testPoint( 0 ) <= xIntersect ) count++;
            }
        }
        resultToBeChecked = ( count & 1 ) != 0;
        results[ i ] = ( testPointsLambdaActuallyPositive[ i ] ) ? resultToBeChecked : false;
    }
    return results;
}
// same as above, implementation to check only one point at the time
bool isPointInTriangle( const ParallelProjection& projection_, const double coordL, const double coordM )
{
    Eigen::Vector2d pointA = projection_.getTriangle2d( ).getVertexA( );
    Eigen::Vector2d pointB = projection_.getTriangle2d( ).getVertexB( );
    Eigen::Vector2d pointC = projection_.getTriangle2d( ).getVertexC( );
    int count = 0;
    double ymin, ymax, xmax;
    double xIntersect;
    // edge from pointA to pointB
    {
        ymin = ( pointA( 1 ) < pointB( 1 ) ) ? pointA( 1 ) : pointB( 1 );
        ymax = ( pointA( 1 ) < pointB( 1 ) ) ? pointB( 1 ) : pointA( 1 );
        xmax = ( pointA( 0 ) > pointB( 0 ) ) ? pointA( 0 ) : pointB( 0 );
        if( ( coordM > ymin ) && ( coordM <= ymax ) && ( coordL <= xmax ) )
        {
            xIntersect = pointA( 0 ) + ( coordM - pointA( 1 ) ) * ( pointB( 0 ) - pointA( 0 ) ) / ( pointB( 1 ) - pointA( 1 ) );
            if( pointA( 0 ) == pointB( 0 ) || coordL <= xIntersect )
            {
                count++;
            }
        }
    }
    // edge from pointB to pointC
    {
        ymin = ( pointB( 1 ) < pointC( 1 ) ) ? pointB( 1 ) : pointC( 1 );
        ymax = ( pointB( 1 ) < pointC( 1 ) ) ? pointC( 1 ) : pointB( 1 );
        xmax = ( pointB( 0 ) > pointC( 0 ) ) ? pointB( 0 ) : pointC( 0 );
        if( ymax == ymin && coordM == ymax )
        {
            count++;
        }
        if( ( coordM > ymin ) && ( coordM <= ymax ) && ( coordL <= xmax ) )
        {
            xIntersect = pointB( 0 ) + ( coordM - pointB( 1 ) ) * ( pointC( 0 ) - pointB( 0 ) ) / ( pointC( 1 ) - pointB( 1 ) );
            if( pointB( 0 ) == pointC( 0 ) || coordL <= xIntersect )
            {
                count++;
            }
        }
    }
    // edge from pointC to pointA
    {
        ymin = ( pointC( 1 ) < pointA( 1 ) ) ? pointC( 1 ) : pointA( 1 );
        ymax = ( pointC( 1 ) < pointA( 1 ) ) ? pointA( 1 ) : pointC( 1 );
        xmax = ( pointC( 0 ) > pointA( 0 ) ) ? pointC( 0 ) : pointA( 0 );
        if( ( coordM > ymin ) && ( coordM <= ymax ) && ( coordL <= xmax ) )
        {
            xIntersect = pointC( 0 ) + ( coordM - pointC( 1 ) ) * ( pointA( 0 ) - pointC( 0 ) ) / ( pointA( 1 ) - pointC( 1 ) );
            if( pointC( 0 ) == pointA( 0 ) || coordL <= xIntersect )
            {
                count++;
            }
        }
    }
    return ( count & 1 ) != 0;
}
// this is a completely new algorithm to check if two edges intersec, used as further testing (not present in literature)
// when the second discrimination does not yield clear result, it is mandatory to check if the edges of the projected 2d triangle intersect
// with the 2d triangle being tested
bool doEdgesIntersect( const Eigen::Vector2d& edge1Start,
                       const Eigen::Vector2d& edge1End,
                       const Eigen::Vector2d& edge2Start,
                       const Eigen::Vector2d& edge2End )
{
    // calculate direction vectors
    Eigen::Vector2d dir1 = edge1End - edge1Start;
    Eigen::Vector2d dir2 = edge2End - edge2Start;

    // calculate the determinant
    double det = dir1.x( ) * dir2.y( ) - dir1.y( ) * dir2.x( );

    // if determinant is zero, lines are parallel or collinear
    if( std::abs( det ) < 1e-9 )
    {
        // check if they are collinear and overlapping
        Eigen::Vector2d vec = edge2Start - edge1Start;
        double crossProduct = vec.x( ) * dir1.y( ) - vec.y( ) * dir1.x( );

        if( std::abs( crossProduct ) < 1e-9 )
        {
            // collinear - check for overlap
            double t0 = vec.dot( dir1 ) / dir1.dot( dir1 );
            double t1 = t0 + dir2.dot( dir1 ) / dir1.dot( dir1 );

            return ( t0 >= 0 && t0 <= 1 ) || ( t1 >= 0 && t1 <= 1 ) || ( t0 <= 0 && t1 >= 1 ) || ( t0 >= 1 && t1 <= 0 );
        }

        // parallel but not collinear
        return false;
    }

    // calculate parameters for the intersection point
    Eigen::Vector2d vec = edge2Start - edge1Start;
    double t = ( vec.x( ) * dir2.y( ) - vec.y( ) * dir2.x( ) ) / det;
    double s = ( vec.x( ) * dir1.y( ) - vec.y( ) * dir1.x( ) ) / det;

    // check if intersection point is within both line segments
    return ( t >= 0 && t <= 1 && s >= 0 && s <= 1 );
}
// accelerated version of the ray-cast PIP algorithm, optimized for large number of points, extremely useful during pixelation
// same exact logic but checks all the points on a line parallel to the ray direction, checking the intersection with the edges, instead of
// all the points one by one
void arePointsInTriangle( const ParallelProjection& projection_,
                          const std::vector< double >& gridCoordinatesL_,
                          const std::vector< double >& gridCoordinatesM_,
                          std::vector< int >& pixelationMatrix_,
                          int indexMinL_,
                          int indexMaxL_,
                          int indexMinM_,
                          int indexMaxM_ )
{
    int len = gridCoordinatesL_.size( );
    int value = 1;
    if( indexMaxL_ == -1 )
    {
        indexMaxL_ = len;
        value = 0;
    }
    if( indexMaxM_ == -1 )
    {
        indexMaxM_ = gridCoordinatesM_.size( );
    }
    double coordM;
    const Eigen::Vector2d pointA = projection_.getTriangle2d( ).getVertexA( );
    const Eigen::Vector2d pointB = projection_.getTriangle2d( ).getVertexB( );
    const Eigen::Vector2d pointC = projection_.getTriangle2d( ).getVertexC( );

    double minLab = ( pointA( 1 ) < pointB( 1 ) ) ? pointA( 1 ) : pointB( 1 );
    double maxLab = ( pointA( 1 ) < pointB( 1 ) ) ? pointB( 1 ) : pointA( 1 );
    double minLbc = ( pointB( 1 ) < pointC( 1 ) ) ? pointB( 1 ) : pointC( 1 );
    double maxLbc = ( pointB( 1 ) < pointC( 1 ) ) ? pointC( 1 ) : pointB( 1 );
    double minLca = ( pointC( 1 ) < pointA( 1 ) ) ? pointC( 1 ) : pointA( 1 );
    double maxLca = ( pointC( 1 ) < pointA( 1 ) ) ? pointA( 1 ) : pointC( 1 );
    double xIntersect;
    double minIntersection, maxIntersection;
    std::vector< double > intersections( 2 );
    int numberOfIntersections;

    for( int i = indexMinM_; i < indexMaxM_; i++ )
    {
        // loop over all y-parallel rays
        coordM = gridCoordinatesM_[ i ];
        numberOfIntersections = 0;
        // edge 1
        if( ( coordM > minLab ) && ( coordM <= maxLab ) )
        {
            xIntersect = pointA( 0 ) + ( coordM - pointA( 1 ) ) * ( pointB( 0 ) - pointA( 0 ) ) / ( pointB( 1 ) - pointA( 1 ) );
            intersections[ numberOfIntersections ] = xIntersect;
            numberOfIntersections += 1;
        }
        // edge 2
        if( ( coordM > minLbc ) && ( coordM <= maxLbc ) )
        {
            xIntersect = pointB( 0 ) + ( coordM - pointB( 1 ) ) * ( pointC( 0 ) - pointB( 0 ) ) / ( pointC( 1 ) - pointB( 1 ) );
            intersections[ numberOfIntersections ] = xIntersect;
            numberOfIntersections += 1;
        }
        // edge 3
        if( ( coordM > minLca ) && ( coordM <= maxLca ) )
        {
            xIntersect = pointC( 0 ) + ( coordM - pointC( 1 ) ) * ( pointA( 0 ) - pointC( 0 ) ) / ( pointA( 1 ) - pointC( 1 ) );
            intersections[ numberOfIntersections ] = xIntersect;
            numberOfIntersections += 1;
        }
        if( numberOfIntersections != 2 )
        {
            continue;
        }
        minIntersection = std::min( intersections[ 0 ], intersections[ 1 ] );
        maxIntersection = std::max( intersections[ 0 ], intersections[ 1 ] );

        if( value == 0 )
        {
            for( int j = indexMinL_; j < indexMaxL_; j++ )
            {
                if( gridCoordinatesL_[ j ] >= minIntersection && gridCoordinatesL_[ j ] <= maxIntersection )
                {
                    pixelationMatrix_[ i * len + j ] = 0;
                }
            }
        }
        else
        {  // value == 1
            for( int j = indexMinL_; j < indexMaxL_; j++ )
            {
                if( gridCoordinatesL_[ j ] >= minIntersection && gridCoordinatesL_[ j ] <= maxIntersection &&
                    pixelationMatrix_[ i * len + j ] == 0 )
                {
                    pixelationMatrix_[ i * len + j ] = 1;
                }
            }
        }
    }
}
// algorithm from "SELF-SHADOWING OF A SPACECRAFT IN THE COMPUTATION OF SURFACE FORCES. AN EXAMPLE IN PLANETARY GEODESY" (Balmino et al.,
// 2018) section 4.3.3 "General case:pixellation"
std::vector< double > computeFractionWithPixelation( const std::vector< std::vector< int > >& sigmaMatrix_,
                                                     const std::vector< std::shared_ptr< VehicleExteriorPanel > >& allPanels_,
                                                     const int maximumNumberOfPixels_,
                                                     const std::vector< std::vector< ParallelProjection > >& projections_,
                                                     const std::vector< int >& toBePixelated_ )
{
    std::vector< double > fractions( toBePixelated_.size( ) );
    std::vector< double > gridCoordinatesL, gridCoordinatesM;
    int index;
    double minL, maxL, minM, maxM;
    double deltaL, deltaM;
    double numberM, numberL;
    double sizeL, sizeM;
    int indexMinL, indexMaxL, indexMinM, indexMaxM;
    double fractionOnes, fractionZeros;
    std::vector< int > pixelationMatrix;
    ParallelProjection selfProjection;
    for( int i = 0; i < static_cast< int >( toBePixelated_.size( ) ); i++ )
    {
        index = toBePixelated_[ i ];
        selfProjection = allPanels_[ index ]->getSelfProjection( );
        minL = selfProjection.getMinimumL( );
        maxL = selfProjection.getMaximumL( );
        minM = selfProjection.getMinimumM( );
        maxM = selfProjection.getMaximumM( );

        deltaL = maxL - minL;
        deltaM = maxM - minM;

        // choosing number of pixels
        if( deltaM >= deltaL )
        {
            numberM = static_cast< double >( maximumNumberOfPixels_ );
            numberL = static_cast< double >( std::max( 2, static_cast< int >( std::round( numberM * deltaL / deltaM ) ) ) );
        }
        else
        {
            numberL = static_cast< double >( maximumNumberOfPixels_ );
            numberM = static_cast< double >( std::max( 2, static_cast< int >( std::round( numberL * deltaM / deltaL ) ) ) );
        }
        sizeL = deltaL / numberL;
        sizeM = deltaM / numberM;
        // creating the grid
        gridCoordinatesL = linspace( minL + sizeL / 2, maxL - sizeL / 2, numberL, gridCoordinatesL );
        gridCoordinatesM = linspace( minM + sizeM / 2, maxM - sizeM / 2, numberM, gridCoordinatesM );
        // resize and reset pixelation matrix to dummy value -1
        if( pixelationMatrix.size( ) != gridCoordinatesM.size( ) * gridCoordinatesL.size( ) )
        {
            pixelationMatrix.resize( gridCoordinatesM.size( ) * gridCoordinatesL.size( ) );
        }
        std::fill( pixelationMatrix.begin( ), pixelationMatrix.end( ), -1 );

        // apply PIP to shadowed plate to find its discretized surface (assign value of 0)
        arePointsInTriangle( allPanels_.at( index )->getSelfProjection( ), gridCoordinatesL, gridCoordinatesM, pixelationMatrix );

        ParallelProjection projection;
        for( int j = 0; j < static_cast< int >( sigmaMatrix_.size( ) ); j++ )
        {
            if( sigmaMatrix_[ index ][ j ] != 0 )
            {
                continue;
            }
            // resize pixelation grid to accelerate the process
            projection = projections_[ index ][ j ];
            indexMinL = clamp( static_cast< int >( std::floor( ( projection.getMinimumL( ) - minL ) / sizeL ) ) - 1,
                               0,
                               static_cast< int >( numberL ) - 1 );
            indexMaxL = clamp( static_cast< int >( std::ceil( ( projection.getMaximumL( ) - minL ) / sizeL ) ) + 1,
                               indexMinL + 1,
                               static_cast< int >( numberL ) );
            indexMinM = clamp( static_cast< int >( std::floor( ( projection.getMinimumM( ) - minM ) / sizeM ) ) - 1,
                               0,
                               static_cast< int >( numberM ) - 1 );
            indexMaxM = clamp( static_cast< int >( std::ceil( ( projection.getMaximumM( ) - minM ) / sizeM ) ) + 1,
                               indexMinM + 1,
                               static_cast< int >( numberM ) );
            // apply PIP to shadowing and update pixelation matrix (assign value of 1)
            arePointsInTriangle(
                    projection, gridCoordinatesL, gridCoordinatesM, pixelationMatrix, indexMinL, indexMaxL, indexMinM, indexMaxM );
        }
        fractionOnes = 0;
        fractionZeros = 0;
        int val;
        for( unsigned int i = 0; i < pixelationMatrix.size( ); i++ )
        {
            val = pixelationMatrix[ i ];
            fractionOnes += ( val == 1 );
            fractionZeros += ( val == 0 );
        }

        fractions[ i ] = clamp( 1 - fractionOnes / ( fractionOnes + fractionZeros ), 0.0, 1.0 );
    }
    return fractions;
}
// core of the SSH algorithm, inspiration from "SELF-SHADOWING OF A SPACECRAFT IN THE COMPUTATION OF SURFACE FORCES. AN EXAMPLE IN PLANETARY
// GEODESY" (Balmino et al., 2018) with many slight modification to take advantage of edge cases or hypotesis during macro-model design
// (e.g., no intersection of panels, only triangular panels) in depth explanation on the revised algorithm can be found in the MSc thesis
// "Non-conservative forces modelling for precise orbit determination" (Maistri, 2025)
void SelfShadowing::updateIlluminatedPanelFractions( const Eigen::Vector3d& incomingDirection )
{
    if( !isComputed_ )
    {
        int numberOfPanels = allPanels_.size( );
        // initialize sigma matrix for discrimination logic (-2 dummy value, -1 no-sh, 0 partial sh, 1 full sh)
        std::vector< std::vector< int > > sigmaMatrix( numberOfPanels, std::vector< int >( numberOfPanels, -2 ) );
        // first discrimination
        std::vector< int > firstDiscriminationIndexes;
        for( int i = 0; i < numberOfPanels; i++ )
        {
            if( !firstDiscrimination( allPanels_.at( i ), incomingDirection ) )
            {
                firstDiscriminationIndexes.push_back( i );
                for( int j = 0; j < numberOfPanels; j++ )
                {
                    sigmaMatrix.at( j ).at( i ) = -1;  // cannot shadow other plates
                }
            }
        }

        for( int i = 0; i < static_cast< int >( firstDiscriminationIndexes.size( ) ); i++ )
        {
            for( int j = 0; j < numberOfPanels; j++ )
            {
                sigmaMatrix.at( firstDiscriminationIndexes.at( i ) ).at( j ) = 1;  // is completely shadowed
            }
        }
        // exclude neighbouring plates
        for( int i = 0; i < numberOfPanels; i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                int neighbour = allPanels_.at( i )->getNeighboringSurfaces( )[ j ];
                if( sigmaMatrix.at( i ).at( neighbour ) == -2 )
                {
                    sigmaMatrix.at( i ).at( neighbour ) = -1;
                }
                if( sigmaMatrix.at( neighbour ).at( i ) == -2 )
                {
                    sigmaMatrix.at( neighbour ).at( i ) = -1;
                }
            }
        }
        // exclude the diagonal (a plate cannot shadow itself)
        for( int i = 0; i < numberOfPanels; i++ )
        {
            if( sigmaMatrix.at( i ).at( i ) == -2 )
            {
                sigmaMatrix.at( i ).at( i ) = -1;
            }
        }
        // general discrimination initialization
        std::vector< std::vector< ParallelProjection > > projections( numberOfPanels, std::vector< ParallelProjection >( numberOfPanels ) );
        std::vector< bool > shadowingInShadowed, shadowedInShadowing;
        // discrimination subcase initialization
        std::vector< bool > areLambdasActuallyPositive;
        std::vector< double > lambdas;
        system_models::Triangle2d triangle2dShadowing;
        Eigen::Vector2d vertexAShadowing, vertexBShadowing, vertexCShadowing;
        system_models::Triangle2d triangle2dShadowed;
        Eigen::Vector2d vertexAShadowed, vertexBShadowed, vertexCShadowed;
        std::vector< Eigen::Vector2d > edgesShadowed;
        std::vector< double > coordinatesL, coordinatesM;
        bool foundPoint;
        double t, lZero, mZero;
        Eigen::Vector2d startEdgeShadowing, endEdgeShadowing;
        const double EPSILON = 1e-12;
        ParallelProjection selfProjection;
        // exclusion logic algorithm
        for( int i = 0; i < numberOfPanels; i++ )
        {
            if( sigmaMatrix.at( i ).at( 0 ) == 1 )
            {
                continue;  // full shadow detected
            }
            // second discrimination
            for( int j = 0; j < numberOfPanels; j++ )
            {
                if( sigmaMatrix.at( i ).at( j ) != -2 )
                {
                    continue;
                }
                ParallelProjection projection(
                        allPanels_.at( i )->getBodyFixedTriangle3d( ), allPanels_.at( j )->getBodyFixedTriangle3d( ), incomingDirection );
                if( secondDiscrimination( projection ) )
                {
                    // testing for those panels that yield p-sh (this case)
                    // first test: min/max coordinates
                    selfProjection = allPanels_.at( i )->getSelfProjection( );
                    if( projection.getMinimumL( ) >= selfProjection.getMaximumL( ) - EPSILON ||
                        projection.getMaximumL( ) <= selfProjection.getMinimumL( ) + EPSILON ||
                        projection.getMinimumM( ) >= selfProjection.getMaximumM( ) - EPSILON ||
                        projection.getMaximumM( ) <= selfProjection.getMinimumM( ) + EPSILON )
                    {
                        sigmaMatrix.at( i ).at( j ) = -1;  // no-sh
                    }
                    else
                    {
                        // second test: use PIP algorithms
                        shadowingInShadowed = isTriangleInTriangle( selfProjection, projection );
                        shadowedInShadowing = isTriangleInTriangle( projection, selfProjection );
                        if( std::any_of( shadowingInShadowed.begin( ), shadowingInShadowed.end( ), []( bool x ) { return x; } ) )
                        {
                            // at least one vertix of the shadowing panel falls into the shadowed, partial sh
                            sigmaMatrix.at( i ).at( j ) = 0;
                            projections.at( i ).at( j ) = projection;
                            if( sigmaMatrix.at( j ).at( i ) == -2 )
                            {
                                sigmaMatrix.at( j ).at( i ) = -1;  // the opposite is a case of no-sh
                            }
                        }
                        if( std::all_of( shadowedInShadowing.begin( ), shadowedInShadowing.end( ), []( bool x ) { return x; } ) )
                        {
                            // all shadowed points are inside the shadowing, full sh
                            for( int u = 0; u < numberOfPanels; u++ )
                            {
                                sigmaMatrix.at( i ).at( u ) = 1;
                            }
                            for( int u = 0; u < numberOfPanels; u++ )
                            {
                                if( sigmaMatrix.at( u ).at( i ) == -2 )
                                {
                                    sigmaMatrix.at( u ).at( i ) = -1;
                                }
                            }
                        }
                        if( std::none_of( shadowingInShadowed.begin( ), shadowingInShadowed.end( ), []( bool x ) { return x; } ) )
                        {
                            // none of the shadowing vertices falls into the shadowed one
                            areLambdasActuallyPositive = projection.getAreLambdasActuallyPositive( );
                            lambdas = projection.getLambdas( );
                            triangle2dShadowing = projection.getTriangle2d( );
                            vertexAShadowing = triangle2dShadowing.getVertexA( );
                            vertexBShadowing = triangle2dShadowing.getVertexB( );
                            vertexCShadowing = triangle2dShadowing.getVertexC( );
                            triangle2dShadowed = allPanels_.at( i )->getSelfProjection( ).getTriangle2d( );
                            vertexAShadowed = triangle2dShadowed.getVertexA( );
                            vertexBShadowed = triangle2dShadowed.getVertexB( );
                            vertexCShadowed = triangle2dShadowed.getVertexC( );
                            edgesShadowed = { vertexAShadowed, vertexBShadowed, vertexCShadowed, vertexAShadowed };
                            coordinatesL = { vertexAShadowing[ 0 ], vertexBShadowing[ 0 ], vertexCShadowing[ 0 ] };
                            coordinatesM = { vertexAShadowing[ 1 ], vertexBShadowing[ 1 ], vertexCShadowing[ 1 ] };
                            foundPoint = false;
                            for( int ii = 0; ii < 3; ii++ )
                            {
                                if( !areLambdasActuallyPositive[ ii ] )
                                {
                                    continue;
                                }
                                for( int jj = 0; jj < 3; jj++ )
                                {
                                    if( ii == jj )
                                    {
                                        continue;
                                    }
                                    // found match, start edge PIP
                                    if( !areLambdasActuallyPositive[ jj ] )
                                    {
                                        t = lambdas[ ii ] / ( lambdas[ ii ] - lambdas[ jj ] );
                                        lZero = coordinatesL[ ii ] + ( coordinatesL[ jj ] - coordinatesL[ ii ] ) * t;
                                        mZero = coordinatesM[ ii ] + ( coordinatesM[ jj ] - coordinatesM[ ii ] ) * t;
                                        startEdgeShadowing = { coordinatesL[ ii ], coordinatesM[ ii ] };
                                        endEdgeShadowing = { lZero, mZero };
                                    }
                                    else
                                    {
                                        startEdgeShadowing = { coordinatesL[ ii ], coordinatesM[ ii ] };
                                        endEdgeShadowing = { coordinatesL[ jj ], coordinatesM[ jj ] };
                                    }

                                    for( int k = 0; k < 3; k++ )
                                    {
                                        foundPoint = doEdgesIntersect(
                                                startEdgeShadowing, endEdgeShadowing, edgesShadowed[ k ], edgesShadowed[ k + 1 ] );
                                        if( foundPoint )
                                        {
                                            break;
                                        }
                                    }
                                    if( foundPoint )
                                    {
                                        break;
                                    }
                                }
                                if( foundPoint )
                                {
                                    break;
                                }
                            }
                            if( foundPoint )
                            {
                                sigmaMatrix.at( i ).at( j ) = 0;
                                projections.at( i ).at( j ) = projection;
                                if( sigmaMatrix.at( j ).at( i ) == -2 )
                                {
                                    sigmaMatrix.at( j ).at( i ) = -1;  // the opposite is a case of no-sh
                                }
                            }
                            else
                            {
                                sigmaMatrix.at( i ).at( j ) = -1;
                            }
                        }
                    }
                }
                else
                {
                    // second discrimination yields no-sh
                    sigmaMatrix.at( i ).at( j ) = -1;
                }
            }
        }
        // pixelation
        std::vector< int > toBePixelated;
        for( int i = 0; i < numberOfPanels; i++ )
        {
            // exclude full shadowed surfaces
            if( sigmaMatrix.at( i ).at( 0 ) == 1 )
            {
                continue;
            }
            toBePixelated.push_back( i );
        }
        int numberOfPanelsToBePixelated = toBePixelated.size( );
        std::vector< double > illuminatedPanelFractions( allPanels_.size( ), 0.0 );
        if( multiThreading_ )
        {
            // option to toggle on multi-treading to split surfaces to be pixelated (currently not user-available)
            int numberOfThreads = std::thread::hardware_concurrency( ) / 2;
            std::vector< std::thread > threads;
            std::vector< std::vector< double > > results( numberOfThreads );
            int chunk = static_cast< int >( std::round( numberOfPanelsToBePixelated / numberOfThreads ) );
            for( int i = 0; i < numberOfThreads - 1; i++ )
            {
                std::vector< int > indexes( toBePixelated.begin( ) + i * chunk, toBePixelated.begin( ) + ( i + 1 ) * chunk );
                threads.emplace_back( [ &, i, indexes ]( ) {
                    results[ i ] = computeFractionWithPixelation( sigmaMatrix, allPanels_, maximumNumberOfPixels_, projections, indexes );
                } );
            }
            int i = numberOfThreads - 1;
            std::vector< int > indexes( toBePixelated.begin( ) + i * chunk, toBePixelated.end( ) );
            threads.emplace_back( [ &, i, indexes ]( ) {
                results[ numberOfThreads - 1 ] =
                        computeFractionWithPixelation( sigmaMatrix, allPanels_, maximumNumberOfPixels_, projections, indexes );
            } );
            for( auto& th: threads )
            {
                th.join( );
            }
            // merge results
            std::vector< double > finalResults;
            for( const auto& res: results )
            {
                finalResults.insert( finalResults.end( ), res.begin( ), res.end( ) );
            }
            for( int i = 0; i < numberOfPanelsToBePixelated; i++ )
            {
                illuminatedPanelFractions.at( toBePixelated.at( i ) ) = finalResults.at( i );
            }
        }
        else
        {
            // default option for single-threading
            std::vector< double > finalResults =
                    computeFractionWithPixelation( sigmaMatrix, allPanels_, maximumNumberOfPixels_, projections, toBePixelated );
            for( int i = 0; i < numberOfPanelsToBePixelated; i++ )
            {
                illuminatedPanelFractions.at( toBePixelated.at( i ) ) = finalResults.at( i );
            }
        }
        illuminatedPanelFractions_ = illuminatedPanelFractions;
    }
    isComputed_ = true;
}

}  // namespace system_models

}  // namespace tudat