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
#include "tudat/astro/system_models/panelGeometryUtils.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"

namespace tudat
{

namespace system_models
{

bool firstDiscrimination( const std::shared_ptr< system_models::VehicleExteriorPanel > panel, const Eigen::Vector3d& v );

bool secondDiscrimination( ParallelProjection& projection );

inline std::vector< bool > isTriangleInTriangle( const ParallelProjection& projection1, const ParallelProjection& projection2 );

bool isPointInTriangle( const ParallelProjection& projection_, const double coordL, const double coordM );

bool doEdgesIntersect( const Eigen::Vector2d& edge1Start,
                       const Eigen::Vector2d& edge1End,
                       const Eigen::Vector2d& edge2Start,
                       const Eigen::Vector2d& edge2End );

void arePointsInTriangle( const ParallelProjection& projection_,
                          const std::vector< double >& gridCoordinatesL_,
                          const std::vector< double >& gridCoordinatesM_,
                          std::vector< int >& pixelationMatrix_,
                          int indexMinL_ = 0,
                          int indexMaxL_ = -1,
                          int indexMinM_ = 0,
                          int indexMaxM_ = -1 );

template< typename T >
inline std::vector< double >& linspace( T start_in, T end_in, int num_in, std::vector< double >& linspaced )
{
    double start = static_cast< double >( start_in );
    double end = static_cast< double >( end_in );
    // resize array
    linspaced.resize( num_in );

    if( num_in == 0 )
    {
        return linspaced;
    }
    if( num_in == 1 )
    {
        linspaced[ 0 ] = start;
        return linspaced;
    }

    double delta = ( end - start ) / ( static_cast< double >( num_in ) - 1 );

    for( int i = 0; i < num_in - 1; ++i )
    {
        linspaced[ i ] = ( start + delta * i );
    }
    linspaced[ num_in - 1 ] = end;

    return linspaced;
}

std::vector< double > computeFractioWithPixelation( const std::vector< std::vector< int > >& sigmaMatrix_,
                                                    const std::vector< std::shared_ptr< VehicleExteriorPanel > >& allPanels_,
                                                    const int maximumNumberOfPixels_,
                                                    const std::vector< std::vector< ParallelProjection > >& projections_,
                                                    const std::vector< int >& toBePixelated_ );

class SelfShadowing
{
private:
    const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels_;
    int maximumNumberOfPixels_;
    bool multiThreading_;
    std::vector< double > illuminatedPanelFractions_;
    bool isComputed_;

public:
    SelfShadowing( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                   const int maximumNumberOfPixels,
                   const bool multiThreading = false ):
        allPanels_( allPanels ), maximumNumberOfPixels_( maximumNumberOfPixels ), multiThreading_( multiThreading ),
        isComputed_( false ){ };

    void updateIlluminatedPanelFractions( const Eigen::Vector3d& incomingDirection );

    std::vector< double >& getIlluminatedPanelFractions( )
    {
        return illuminatedPanelFractions_;
    }

    int getMaximumNumberOfPixels( ) const
    {
        return maximumNumberOfPixels_;
    }

    void reset( )
    {
        isComputed_ = false;
    };
};

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_SELFSHADOWING_H
