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

namespace tudat
{

namespace system_models
{


class VehicleExteriorPanel
{
public:

    VehicleExteriorPanel(
        const Eigen::Vector3d& frameFixedSurfaceNormal,
        const Eigen::Vector3d& frameFixedPositionVector,
        const double panelArea,
        const double panelTemperature = 273.0,
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [=]( ){ return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [=]( ){ return frameFixedPositionVector; } ),
        panelArea_( panelArea ),
        panelTemperature_( panelTemperature ),
        trackedBody_( "" ),
        reflectionLaw_( reflectionLaw ){ }

    VehicleExteriorPanel(
        const Eigen::Vector3d& frameFixedSurfaceNormal,
        const Eigen::Vector3d& frameFixedPositionVector,
        const double panelArea,
        const double panelTemperature = 273.0,
        const std::string trackedBody = "",
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [=]( ){ return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [=]( ){ return frameFixedPositionVector; } ),
        panelArea_( panelArea ),
        panelTemperature_( panelTemperature ),
        trackedBody_( trackedBody ),
        reflectionLaw_( reflectionLaw ){ }

    VehicleExteriorPanel(
        const std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal,
        const std::function< Eigen::Vector3d( ) > frameFixedPositionVector,
        const double panelArea,
        const double panelTemperature = 273.0,
        const std::string trackedBody = "",
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( frameFixedSurfaceNormal ),
        frameFixedPositionVector_( frameFixedPositionVector ),
        panelArea_( panelArea ),
        panelTemperature_( panelTemperature ),
        trackedBody_( trackedBody ),
        reflectionLaw_( reflectionLaw ){ }

        VehicleExteriorPanel(
        const double panelArea,
        const Eigen::Vector3d& frameFixedSurfaceNormal,
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( [=]( ){ return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [=]( ){ return Eigen::Vector3d::Constant( TUDAT_NAN ); } ),
        panelArea_( panelArea ),
        panelTemperature_( TUDAT_NAN ),
        trackedBody_( "" ),
        reflectionLaw_( reflectionLaw ){ }

    VehicleExteriorPanel(
        const Eigen::Vector3d& frameFixedSurfaceNormal,
        const double panelArea,
        const std::string trackedBody = "",
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr,
        const Eigen::Vector3d& frameFixedPosition = Eigen::Vector3d::Constant( TUDAT_NAN ) ):
        frameFixedSurfaceNormal_( [=]( ){ return frameFixedSurfaceNormal; } ),
        frameFixedPositionVector_( [=]( ){ return frameFixedPosition; } ),
        panelArea_( panelArea ),
        panelTemperature_( TUDAT_NAN ),
        trackedBody_( trackedBody ),
        reflectionLaw_( reflectionLaw ){ }

    VehicleExteriorPanel(
        const std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal,
        const double panelArea,
        const std::string trackedBody = "",
        const std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw = nullptr ):
        frameFixedSurfaceNormal_( frameFixedSurfaceNormal ),
        panelArea_( panelArea ),
        trackedBody_( trackedBody ),
        reflectionLaw_( reflectionLaw ){ }

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

protected:

    std::function< Eigen::Vector3d( ) > frameFixedSurfaceNormal_;

    std::function< Eigen::Vector3d( ) > frameFixedPositionVector_;

    double panelArea_;

    double panelTemperature_;

    std::string trackedBody_;

    std::shared_ptr< electromagnetism::ReflectionLaw > reflectionLaw_;

    std::string panelTypeId_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_VEHICLEEXTERIORPANELS_H
