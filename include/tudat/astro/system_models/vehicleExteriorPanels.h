/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VEHICLEEXTERIORPANELS_H
#define TUDAT_VEHICLEEXTERIORPANELS_H

#include <map>
#include <iostream>

#include <memory>

#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/system_models/engineModel.h"

#include "tudat/astro/system_models/panelGeometryUtils.h"

namespace tudat
{

namespace system_models
{

class VehicleExteriorPanel
{
public:
    VehicleExteriorPanel( const Eigen::Vector3d& frameFixedSurfaceNormal,
                          const Eigen::Vector3d& frameFixedPositionVector,
                          const double panelArea,
                          const double panelTemperature = 273.0,
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [ = ]( ) { return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [ = ]( ) { return frameFixedPositionVector; } ), panelArea_( panelArea ),
        panelTemperature_( panelTemperature ), trackedBody_( "" ), reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const Eigen::Vector3d& frameFixedSurfaceNormal,
                          const Eigen::Vector3d& frameFixedPositionVector,
                          const double panelArea,
                          const double panelTemperature = 273.0,
                          const std::string trackedBody = "",
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [ = ]( ) { return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [ = ]( ) { return frameFixedPositionVector; } ), panelArea_( panelArea ),
        panelTemperature_( panelTemperature ), trackedBody_( trackedBody ), reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal,
                          const std::function< Eigen::Vector3d( ) > frameFixedPositionVector,
                          const double panelArea,
                          const double panelTemperature = 273.0,
                          const std::string trackedBody = "",
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( frameFixedSurfaceNormal ), frameFixedPositionVector_( frameFixedPositionVector ), panelArea_( panelArea ),
        panelTemperature_( panelTemperature ), trackedBody_( trackedBody ), reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const double panelArea,
                          const Eigen::Vector3d& frameFixedSurfaceNormal,
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [ = ]( ) { return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [ = ]( ) { return Eigen::Vector3d::Constant( TUDAT_NAN ); } ), panelArea_( panelArea ),
        panelTemperature_( TUDAT_NAN ), trackedBody_( "" ), reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const Eigen::Vector3d& frameFixedSurfaceNormal,
                          const double panelArea,
                          const std::string trackedBody = "",
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr,
                          const Eigen::Vector3d& frameFixedPosition = Eigen::Vector3d::Constant( TUDAT_NAN ) ):
        frameFixedSurfaceNormal_( [ = ]( ) { return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [ = ]( ) { return frameFixedPosition; } ), panelArea_( panelArea ), panelTemperature_( TUDAT_NAN ),
        trackedBody_( trackedBody ), reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal,
                          const double panelArea,
                          const std::string trackedBody = "",
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( frameFixedSurfaceNormal ),
        frameFixedPositionVector_( [ = ]( ) { return Eigen::Vector3d::Constant( TUDAT_NAN ); } ),
        panelArea_( panelArea ), trackedBody_( trackedBody ),
        reflectionLaw_( reflectionLaw ), frameOrigin_( Eigen::Vector3d::Zero( ) ), geometryDefined_( false )
    { }

    VehicleExteriorPanel( const std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal,
                          const std::function< Eigen::Vector3d( ) > frameFixedPositionVector,
                          const double panelArea,
                          const double panelTemperature = 273.0,
                          const std::string trackedBody = "",
                          const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr,
                          const Triangle3d& triangle3d = Triangle3d( ),
                          const Eigen::Vector3d frameOrigin = Eigen::Vector3d::Zero( ),
                          const bool geometryDefined = false ):
        frameFixedSurfaceNormal_( frameFixedSurfaceNormal ), frameFixedPositionVector_( frameFixedPositionVector ), panelArea_( panelArea ),
        panelTemperature_( panelTemperature ), trackedBody_( trackedBody ), reflectionLaw_( reflectionLaw ), 
        triangle3d_( triangle3d ), frameOrigin_( frameOrigin ), geometryDefined_( geometryDefined )
    { 
        neighboringSurfaces_ = {-1, -1, -1};
        if ( geometryDefined_ )
        {
            selfProjection_ = ParallelProjection( triangle3d_, triangle3d_, frameFixedSurfaceNormal_( ) );
            std::vector< bool > dummyVector = { true, true, true };
            selfProjection_.setAreLambdasActuallyPositive( dummyVector );
        }
        else
        {
            selfProjection_ = ParallelProjection( );
        }
    }
    // copy-rotation constructor
    VehicleExteriorPanel( const std::shared_ptr< VehicleExteriorPanel > panelToBeRotated, 
        const Eigen::Quaterniond& currentOrientation_ )
    {
        frameFixedSurfaceNormal_ = [=] () { return currentOrientation_*panelToBeRotated->getFrameFixedSurfaceNormal( )( ); };
        frameFixedPositionVector_ = [=] () { 
            Eigen::Vector3d frameFixedPositionRotated = currentOrientation_*panelToBeRotated->getFrameFixedPositionVector( )( ) + panelToBeRotated->getFrameOrigin( );
            return frameFixedPositionRotated; 
        };
        panelArea_ = panelToBeRotated->getPanelArea( );
        panelTemperature_ = panelToBeRotated->getPanelTemperature( );
        trackedBody_ = panelToBeRotated->getTrackedBody( );
        reflectionLaw_ = panelToBeRotated->getReflectionLaw( );
        Eigen::Vector3d vertexA = currentOrientation_*panelToBeRotated->getTriangle3d( ).getVertexA( ) + 
            panelToBeRotated->getFrameOrigin( );
        Eigen::Vector3d vertexB = currentOrientation_*panelToBeRotated->getTriangle3d( ).getVertexB( ) + 
            panelToBeRotated->getFrameOrigin( );
        Eigen::Vector3d vertexC = currentOrientation_*panelToBeRotated->getTriangle3d( ).getVertexC( ) + 
            panelToBeRotated->getFrameOrigin( );
        frameOrigin_ = panelToBeRotated->getFrameOrigin( );
        triangle3d_ = Triangle3d( vertexA, vertexB, vertexC );
        selfProjection_ = panelToBeRotated->getSelfProjection( );
        neighboringSurfaces_ = panelToBeRotated->getNeighboringSurfaces( );
        geometryDefined_ = true;
    }


    void setReflectionLaw( const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw )
    {
        reflectionLaw_ = reflectionLaw;
    }

    std::shared_ptr< electromagnetism::ReflectionLaw > getReflectionLaw( ) const
    {
        return reflectionLaw_;
    }

    std::function< Eigen::Vector3d( ) > getFrameFixedSurfaceNormal( ) const
    {
        return frameFixedSurfaceNormal_;
    }

    std::function< Eigen::Vector3d( ) > getFrameFixedPositionVector( ) const
    {
        return frameFixedPositionVector_;
    }

    double getPanelArea( ) const
    {
        return panelArea_;
    }

    double getPanelTemperature( ) const
    {
        return panelTemperature_;
    }

    std::string getTrackedBody( )
    {
        return trackedBody_;
    }

    std::string getPanelTypeId( ) const
    {
        return panelTypeId_;
    }

    void setPanelTypeId( const std::string panelTypeId )
    {
        panelTypeId_ = panelTypeId;
    }

    // addition for self-shadowing (SSH) algorithm
    Triangle3d getTriangle3d( ) const
    {
        return triangle3d_;
    }
    std::vector<int> getNeighboringSurfaces( ) const
    {
        return neighboringSurfaces_;
    }
    void setNeighboringSurfaces( const std::vector<int>& neighboringSurfaces ) 
    {
        neighboringSurfaces_ = neighboringSurfaces;
    }
    bool isGeometryDefined( ) const
    {
        return geometryDefined_;
    }
    ParallelProjection getSelfProjection( ) const
    {
        return selfProjection_;
    }

    Eigen::Vector3d getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

protected:
    std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal_;

    std::function< Eigen::Vector3d( ) > frameFixedPositionVector_;

    double panelArea_;

    double panelTemperature_;

    std::string trackedBody_;

    std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw_;

    std::string panelTypeId_;

    // addition for self-shadowing (SSH) algorithm
    Triangle3d triangle3d_;

    std::vector< int > neighboringSurfaces_;

    ParallelProjection selfProjection_;

    Eigen::Vector3d frameOrigin_;

    bool geometryDefined_;

};

}  // namespace system_models

}  // namespace tudat

#endif  // TUDAT_VEHICLEEXTERIORPANELS_H
