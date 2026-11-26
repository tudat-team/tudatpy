/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Eigen/LU>

#include "tudat/astro/system_models/camera.h"
namespace tudat
{
namespace system_models
{

Eigen::Vector2d PsfCameraProjectionModel::pixelLineToFocalPlaneCoordinates( const Eigen::Vector2d& pixelLine ) const
{
    const Eigen::Vector2d target = pixelLine - principalPoint_;
    Eigen::Vector2d focalPlaneCoordinates = kMatrix_.block< 2, 2 >( 0, 0 ).fullPivLu( ).solve( target );

    for( unsigned int i = 0; i < 20; ++i )
    {
        const double x = focalPlaneCoordinates.x( );
        const double y = focalPlaneCoordinates.y( );
        const Eigen::Vector2d residual = focalPlaneCoordinatesToPixelLine( x, y ) - pixelLine;
        if( residual.norm( ) < 1.0E-13 )
        {
            return focalPlaneCoordinates;
        }

        Eigen::Matrix2d jacobian;
        jacobian( 0, 0 ) = kMatrix_( 0, 0 ) + kMatrix_( 0, 2 ) * y;
        jacobian( 0, 1 ) = kMatrix_( 0, 1 ) + kMatrix_( 0, 2 ) * x;
        jacobian( 1, 0 ) = kMatrix_( 1, 0 ) + kMatrix_( 1, 2 ) * y;
        jacobian( 1, 1 ) = kMatrix_( 1, 1 ) + kMatrix_( 1, 2 ) * x;

        focalPlaneCoordinates -= jacobian.fullPivLu( ).solve( residual );
    }

    throw std::runtime_error( "Error when converting PSF pixel/line coordinates to camera-frame unit vector: Newton solve failed." );
}

}  // namespace system_models
}  // namespace tudat
