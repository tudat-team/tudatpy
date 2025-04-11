/*    Copyright (c) 2010-2019, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/
#ifndef TUDAT_SELFSHADOWING_H
#define TUDAT_SELFSHADOWING_H

#include <map>
#include <iostream>
#include <algorithm>
#include <memory>

#include <thread>
#include "tudat/astro/system_models/geometryUtils.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"

namespace tudat
{

namespace system_models
{

inline bool firstDiscrimination( const std::shared_ptr< system_models::VehicleExteriorPanel > panel, 
                          const Eigen::Vector3d& v) {
    if ( panel->getFrameFixedSurfaceNormal( )( ).dot(v) < 0 ) 
    {
        return 1;
    }
    else {
        return 0;
    }
}

inline bool secondDiscrimination( const ParallelProjection& projection )
{
    std::vector< double > lambdas = projection.getLambdas( );
    if ( *std::max_element( lambdas.begin( ), lambdas.end( )) <= 0.0 )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

inline std::vector< bool > isTriangleInTriangle( const ParallelProjection& projection1, const ParallelProjection& projection2 )
{
    std::vector< Eigen::Vector2d > testPoints = { projection2.getTriangle2d( ).getVertexA( ), 
                                                  projection2.getTriangle2d( ).getVertexB( ),
                                                  projection2.getTriangle2d( ).getVertexC( ) };
    std::vector< bool > results( 3 );
    for ( int i=0; i<3; i++ )
    {
        int count = 0;
        Eigen::Vector2d testPoint = testPoints[i];
        Eigen::Vector2d pointA = projection1.getTriangle2d( ).getVertexA( );
        Eigen::Vector2d pointB = projection1.getTriangle2d( ).getVertexB( );
        Eigen::Vector2d pointC = projection1.getTriangle2d( ).getVertexC( );
        // edge from pointA to pointB
        {
            double ymin = (pointA( 1 ) < pointB( 1 )) ? pointA( 1 ) : pointB( 1 );
            double ymax = (pointA( 1 ) < pointB( 1 )) ? pointB( 1 ) : pointA( 1 );
            double xmax = (pointA( 0 ) > pointB( 0 )) ? pointA( 0 ) : pointB( 0 );
            if ((testPoint( 1 ) > ymin) && (testPoint( 1 ) <= ymax) && (testPoint( 0 ) <= xmax)) {
                double xIntersect = pointA( 0 ) + (testPoint( 1 ) - pointA( 1 )) * 
                (pointB( 0 ) - pointA( 0 )) / (pointB( 1 ) - pointA( 1 ));
                if (pointA( 0 ) == pointB( 0 ) || testPoint( 0 ) <= xIntersect)
                    count++;
            }
        }
        // edge from pointB to pointC
        {
            double ymin = (pointB( 1 ) < pointC( 1 )) ? pointB( 1 ) : pointC( 1 );
            double ymax = (pointB( 1 ) < pointC( 1 )) ? pointC( 1 ) : pointB( 1 );
            double xmax = (pointB( 0 ) > pointC( 0 )) ? pointB( 0 ) : pointC( 0 );
            if (ymax == ymin && testPoint( 1 ) == ymax) {
                count++;
            }
            if ((testPoint( 1 ) > ymin) && (testPoint( 1 ) <= ymax) && (testPoint( 0 ) <= xmax)) {
                double xIntersect = pointB( 0 ) + (testPoint( 1 ) - pointB( 1 )) * 
                (pointC( 0 ) - pointB( 0 )) / (pointC( 1 ) - pointB( 1 ));
                if (pointB( 0 ) == pointC( 0 ) || testPoint( 0 ) <= xIntersect)
                    count++;
            }
        }
        // edge from pointC to pointA
        {
            double ymin = (pointC( 1 ) < pointA( 1 )) ? pointC( 1 ) : pointA( 1 );
            double ymax = (pointC( 1 ) < pointA( 1 )) ? pointA( 1 ) : pointC( 1 );
            double xmax = (pointC( 0 ) > pointA( 0 )) ? pointC( 0 ) : pointA( 0 );
            if ((testPoint( 1 ) > ymin) && (testPoint( 1 ) <= ymax) && (testPoint( 0 ) <= xmax)) {
                double xIntersect = pointC( 0 ) + (testPoint( 1 ) - pointC( 1 )) * 
                (pointA( 0 ) - pointC( 0 )) / (pointA( 1 ) - pointC( 1 ));
                if (pointC( 0 ) == pointA( 0 ) || testPoint( 0 ) <= xIntersect)
                    count++;
            }
        }
        results[i] = ( count & 1 ) != 0; 
    }
    return results;
}

inline std::vector< std::vector< int > > arePointsInTriangle( const ParallelProjection& projection_, 
    const std::vector<double>& gridCoordinatesL_, const std::vector<double>& gridCoordinatesM_,
    const std::vector< std::vector < int > >& pixelationMatrix_,
    int indexMinL_ = 0, int indexMaxL_ = -1, int indexMinM_ = 0, int indexMaxM_ = -1 )
{   
    std::vector< std::vector < int > > pixelationMatrixUpdated = pixelationMatrix_;
    int value = 1;
    if ( indexMaxL_ == -1 )
    {
        indexMaxL_ = gridCoordinatesL_.size( );
        value = 0;
    }
    if ( indexMaxM_ == -1 )
    {
        indexMaxM_ = gridCoordinatesM_.size( );
    }
    
    for ( int i=indexMinM_; i<indexMaxM_; i++ )
    {
        // loop over all y-parallel rays
        const double coordM = gridCoordinatesM_[i];
        const Eigen::Vector2d pointA = projection_.getTriangle2d( ).getVertexA( );
        const Eigen::Vector2d pointB = projection_.getTriangle2d( ).getVertexB( );
        const Eigen::Vector2d pointC = projection_.getTriangle2d( ).getVertexC( );
        std::vector< double > intersections;

        // edge from pointA to pointB
        double minLab = (pointA( 1 ) < pointB( 1 )) ? pointA( 1 ) : pointB( 1 );
        double maxLab = (pointA( 1 ) < pointB( 1 )) ? pointB( 1 ) : pointA( 1 );
        if ((coordM > minLab) && (coordM <= maxLab)) {
            double xIntersect = pointA( 0 ) + (coordM - pointA( 1 )) * (pointB( 0 ) - pointA( 0 )) / (pointB( 1 ) - pointA( 1 ));
            intersections.push_back(xIntersect);
        }
        // edge 2
        double minLbc = (pointB( 1 ) < pointC( 1 )) ? pointB( 1 ) : pointC( 1 );
        double maxLbc = (pointB( 1 ) < pointC( 1 )) ? pointC( 1 ) : pointB( 1 );
        if ((coordM > minLbc) && (coordM <= maxLbc)) {
            double xIntersect = pointB( 0 ) + (coordM - pointB( 1 )) * (pointC( 0 ) - pointB( 0 )) / (pointC( 1 ) - pointB( 1 ));
            intersections.push_back(xIntersect);
        }
        // edge 3
        double minLca = (pointC( 1 ) < pointA( 1 )) ? pointC( 1 ) : pointA( 1 );
        double maxLca = (pointC( 1 ) < pointA( 1 )) ? pointA( 1 ) : pointC( 1 );
        if ((coordM > minLca) && (coordM <= maxLca)) {
            double xIntersect = pointC( 0 ) + (coordM - pointC( 1 )) * (pointA( 0 ) - pointC( 0 )) / (pointA( 1 ) - pointC( 1 ));
            intersections.push_back(xIntersect);
        }
        if (intersections.empty()) {
            continue;
        }
        double minIntersection = *std::min_element(intersections.begin(), intersections.end());
        double maxIntersection = *std::max_element(intersections.begin(), intersections.end());
        for (int j = indexMinL_; j<indexMaxL_; j++)
        {
            if ( gridCoordinatesL_[j] >= minIntersection && gridCoordinatesL_[j] <= maxIntersection && 
                value == 0)
            {
                pixelationMatrixUpdated[i][j] = value;
            }
            if ( gridCoordinatesL_[j] >= minIntersection && gridCoordinatesL_[j] <= maxIntersection && 
                value == 1 && pixelationMatrixUpdated[i][j] == 0)
            {
                pixelationMatrixUpdated[i][j] = value;
            }
        }
    }
    return pixelationMatrixUpdated;

}

template<typename T>
inline std::vector<double> linspace(T start_in, T end_in, int num_in)
{
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 

  return linspaced;
}

inline std::vector< double > computeFractioWithPixelation( const std::vector<std::vector<int>>& sigmaMatrix_, 
    const std::vector< std::shared_ptr< VehicleExteriorPanel > >& allPanels_, 
    const int maximumNumberOfPixels_, const std::vector< std::vector < ParallelProjection > >& projections_, 
    const std::vector< int > toBePixelated_ )
{
    std::vector< double > fractions( toBePixelated_.size( ) );
    for ( int i=0; i<static_cast< int >( toBePixelated_.size( ) ); i++ )
    {
        int index = toBePixelated_[i];
        double minL = allPanels_[index]->getSelfProjection( ).getMinL( );
        double maxL = allPanels_[index]->getSelfProjection( ).getMaxL( );
        double minM = allPanels_[index]->getSelfProjection( ).getMinM( );
        double maxM = allPanels_[index]->getSelfProjection( ).getMaxM( );

        double deltaL = maxL - minL;
        double deltaM = maxM - minM;

        // choosing number of pixels
        double numberM, numberL;
        if (deltaM>=deltaL) {
            numberM = static_cast< double >( maximumNumberOfPixels_ );
            numberL = static_cast< double >( std::max(2, static_cast< int >(std::round( numberM * deltaL/deltaM))));
        }
        else {
            numberL = static_cast< double >( maximumNumberOfPixels_ );
            numberM = static_cast< double >( std::max(2, static_cast< int >(std::round( numberL * deltaM/deltaL))));
        }
        double sizeL = deltaL / numberL;
        double sizeM = deltaM / numberM;
        // creating the grid
        std::vector< double > gridCoordinatesL = linspace(minL + sizeL/2, maxL - sizeL/2, numberL );
        std::vector< double > gridCoordinatesM = linspace(minM + sizeM/2, maxM - sizeM/2, numberM );
        // initialize pixelation matrix with -1 values (W in literature)
        std::vector< std::vector< int > > pixelationMatrix( gridCoordinatesM.size(), std::vector< int >(gridCoordinatesL.size(), -1));
        // apply PIP to shadowed plate to find its discretized surface (assign value of 0)
        pixelationMatrix = arePointsInTriangle( allPanels_[index]->getSelfProjection( ), 
            gridCoordinatesL, gridCoordinatesM, pixelationMatrix );

        for ( int j = 0; j< static_cast< int >( sigmaMatrix_.size( ) ); j++)
        {
            if (sigmaMatrix_[index][j] != 0)
            {
                continue;
            }
            // resize pixelation grid to accelerate the process
            int indexMinL = clamp(static_cast< int >( 
                std::floor((projections_[index][j].getMinL( )-minL)/sizeL)) - 1, 0, static_cast< int >( numberL ) - 1);
            int indexMaxL = clamp(static_cast< int >( 
                std::ceil((projections_[index][j].getMaxL( )-minL)/sizeL)) + 1, indexMinL + 1, static_cast< int >( numberL ));
            int indexMinM = clamp(static_cast< int >(
                std::floor((projections_[index][j].getMinM( )-minM)/sizeM)) - 1, 0, static_cast< int >( numberM ) - 1);
            int indexMaxM = clamp(static_cast< int >(
                std::ceil((projections_[index][j].getMaxM( )-minM)/sizeM)) + 1, indexMinM + 1, static_cast< int >( numberM ));
            // apply PIP to shadowing and update pixelation matrix (assign value of 1)  
            pixelationMatrix = arePointsInTriangle(projections_[index][j], 
                gridCoordinatesL, gridCoordinatesM, pixelationMatrix, indexMinL, indexMaxL, indexMinM, indexMaxM);
        }
        double fractionOnes = 0;
        for (int p = 0; p< static_cast< int >( gridCoordinatesM.size() ); p++) 
        {
            for (int q = 0; q< static_cast< int >( gridCoordinatesL.size() ); q++) 
            {
                if (pixelationMatrix[p][q] == 1) 
                {
                    fractionOnes++;
                }
            }
        }
        double fractionZeros = 0;
        for (int p = 0; p< static_cast< int >( gridCoordinatesM.size() ); p++) 
        {
            for (int q = 0; q< static_cast< int >( gridCoordinatesL.size() ); q++) 
            {
                if (pixelationMatrix[p][q] == 0) 
                {
                    fractionZeros++;
                }
            }
        }
        
        fractions[i] = clamp( 1 - fractionOnes / (fractionOnes + fractionZeros ), 0.0, 1.0 );
    }
    return fractions;
}

class SelfShadowing
{
private:
    std::vector< double > illuminatedPanelFractions_;
public:
    SelfShadowing( ) : illuminatedPanelFractions_( ) { };

    SelfShadowing( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                   const Eigen::Vector3d& incomingDirection, 
                   const int maximumNumberOfPixels,
                   const bool multiThreading = false )
    { 
        int numberOfPanels = allPanels.size( );
        // initialize sigma matrix for discrimination logic (-2 dummy value, -1 no-sh, 0 partial sh, 1 full sh)
        std::vector<std::vector< int >> sigmaMatrix(numberOfPanels, std::vector< int >(numberOfPanels, -2));
        // first discrimination
        std::vector< int > firstDiscriminationIndexes;
        for (int i = 0; i < numberOfPanels; i++) {
            if ( !firstDiscrimination( allPanels[i], incomingDirection ) ){
                firstDiscriminationIndexes.push_back( i );
                for (int j = 0; j < numberOfPanels; j++) {
                    sigmaMatrix[j][i] = -1; // cannot shadow other plates
                }
            }
        }
        
        for ( int i = 0; i < static_cast< int >( firstDiscriminationIndexes.size( ) ); i++ ) {
            for ( int j = 0; j < numberOfPanels; j++ ) {
                sigmaMatrix[firstDiscriminationIndexes[i]][j] = 1; // is completely shadowed
            }
        }
        // exclude neighbouring plates
        for ( int i=0; i<numberOfPanels; i++ ) {
            for ( int j=0; j<3; j++ ) {
                int neighbour = allPanels[i]->getNeighboringSurfaces( )[j];
                if ( sigmaMatrix[i][neighbour]==-2 ) 
                {
                    sigmaMatrix[i][neighbour] = -1;
                }
                if (sigmaMatrix[neighbour][i]==-2) 
                {
                    sigmaMatrix[neighbour][i] = -1;
                }
            }
        }
        // exclude the diagonal (a plate cannot shadow itself)
        for  (int i = 0; i < numberOfPanels; i++)
        {
            if (sigmaMatrix[i][i] == -2) 
            {
                sigmaMatrix[i][i] = -1;
            }
        }
        std::vector< std::vector < ParallelProjection > > projections( numberOfPanels, 
                                                                       std::vector < ParallelProjection >( numberOfPanels ) );
        // exclusion logic algorithm
        for ( int i = 0; i<numberOfPanels; i++ )
        {
            if ( sigmaMatrix[i][0] == 1 )
            {
                continue; // full shadow detected
            }
            // second discrimination
            for ( int j = 0; j<numberOfPanels; j++ )
            {
                if ( sigmaMatrix[i][j] != -2 )
                {
                    continue;
                }
                ParallelProjection projection( allPanels[i]->getTriangle3d( ), allPanels[j]->getTriangle3d( ), incomingDirection );
                if ( secondDiscrimination( projection ) )
                {
                    // testing for those panels that yield p-sh (this case)
                        // first test: min/max coordinates
                    if ( projection.getMinL( ) >= allPanels[i]->getSelfProjection( ).getMaxL( ) ||
                         projection.getMaxL( ) <= allPanels[i]->getSelfProjection( ).getMinL( ) || 
                         projection.getMinM( ) >= allPanels[i]->getSelfProjection( ).getMaxM( ) ||
                         projection.getMaxM( ) <= allPanels[i]->getSelfProjection( ).getMinM( ) )
                    {
                        sigmaMatrix[i][j] = -1; // no-sh
                    }
                    else
                    {
                        // second test: use PIP algorithms
                        std::vector< bool > shadowingInShadowed = isTriangleInTriangle( allPanels[i]->getSelfProjection( ), 
                                                                                        projection );
                        std::vector< bool > shadowedInShadowing = isTriangleInTriangle( projection,
                                                                                        allPanels[i]->getSelfProjection( ) );
                        if ( std::any_of( shadowingInShadowed.begin( ), shadowingInShadowed.end( ), [](bool x) { return x; }))
                        {
                            // at least one vertix of the shadowing panel falls into the shadowed, partial sh
                            sigmaMatrix[i][j] = 0; 
                            projections[i][j] = projection;
                            if ( sigmaMatrix[j][i] == -2 )
                            {
                                sigmaMatrix[j][i] = -1; // the opposite is a case of no-sh
                            }
                        }
                        if ( std::all_of( shadowedInShadowing.begin( ), shadowedInShadowing.end( ), [](bool x) { return x; }))
                        {   
                            // all shadowed points are inside the shadowing, full sh
                            for ( int u=0; u<numberOfPanels; u++ )
                            {
                                sigmaMatrix[i][u] = 1; 
                            }
                            for ( int u=0; u<numberOfPanels; u++ )
                            {
                                if ( sigmaMatrix[u][i] == -2 )
                                {
                                    sigmaMatrix[u][i] = -1; 
                                }
                            }
                        }
                        else
                        {
                            // otherwise cannot decide, needs to be pixelated when uncertain
                            sigmaMatrix[i][j] = 0; 
                            projections[i][j] = projection;                            
                        }
                    }
                }
                else
                {
                    // second discrimination yields no-sh
                    sigmaMatrix[i][j] = -1;
                }
            }
        } 
        // pixelation 
        std::vector<int> toBePixelated;
        for ( int i = 0; i<numberOfPanels; i++ ) 
        {
            // exclude full shadowed surfaces
            if (sigmaMatrix[i][0] == 1) {
                continue; 
            }
            toBePixelated.push_back(i);    
        } 
        int numberOfPanelsToBePixelated = toBePixelated.size( );
        std::vector< double > illuminatedPanelFractions( allPanels.size( ), 0.0 );
        if ( multiThreading ) 
        {
            // option to toggle on multi-treading to split surfaces to be pixelated (currently not user-available)
            int numberOfThreads = std::thread::hardware_concurrency( ) / 2;
            std::vector< std::thread > threads;
            std::vector< std::vector< double > > results( numberOfThreads );
            int chunk = static_cast< int >( std::round( numberOfPanelsToBePixelated / numberOfThreads ) );
            for (int i = 0; i<numberOfThreads-1; i++) {
                std::vector<int> indexes( toBePixelated.begin() + i*chunk, toBePixelated.begin() + (i+1)*chunk);
                threads.emplace_back([&, i, indexes]() {
                    results[i] = computeFractioWithPixelation( sigmaMatrix, allPanels, maximumNumberOfPixels, projections,
                                                        indexes );
                });
            }
            int i = numberOfThreads - 1;
            std::vector<int> indexes(toBePixelated.begin() + i*chunk, toBePixelated.end());
            threads.emplace_back([&, i, indexes]() {
                results[numberOfThreads-1] = computeFractioWithPixelation( sigmaMatrix, allPanels, maximumNumberOfPixels, projections,
                    indexes );
            });
            for (auto& th : threads)
            {
                th.join();
            }
            // merge results
            std::vector<double> finalResults;
            for (const auto& res : results) {
                finalResults.insert( finalResults.end(), res.begin(), res.end() );
            }
            for (int i = 0; i<numberOfPanelsToBePixelated; i++)
            {
                illuminatedPanelFractions[toBePixelated[i]] = finalResults[i];
            }

        }
        else
        {
            // default option for single-threading
            std::vector<double> finalResults = computeFractioWithPixelation( sigmaMatrix, allPanels, maximumNumberOfPixels, projections,
                toBePixelated );
                
            for (int i = 0; i<numberOfPanelsToBePixelated; i++)
            {
                illuminatedPanelFractions[toBePixelated[i]] = finalResults[i];
            }
        }      
        illuminatedPanelFractions_ = illuminatedPanelFractions;  
    }    
    
    std::vector< double > getIlluminatedPanelFractions( )
    {
        return illuminatedPanelFractions_;
    }
};


}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_SELFSHADOWING_H
