/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_METEOROLOGICALCONDITIONS_H
#define TUDAT_METEOROLOGICALCONDITIONS_H

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{

namespace ground_stations
{

// Function to compute saturation vapor pressure (Pa) using Tetens' formula
double computeSaturationWaterVaporPressure( const double temperature);

double computeDewPoint( const double waterVaporPartialPressure);

class StationMeteoData
{
public:
    StationMeteoData( ){ }

    virtual ~StationMeteoData( ){ }

    virtual double getTemperature( const double currentUtc ) = 0;

    virtual double getPressure( const double currentUtc ) = 0;

    virtual double getWaterVaporPartialPressure( const double currentUtc ) = 0;

    virtual double getRelativeHumidity( const double currentUtc )
    {
        return getWaterVaporPartialPressure( currentUtc ) /
        computeSaturationWaterVaporPressure( getTemperature( currentUtc ) );
    }

    virtual double getDewPointTemperature( const double currentUtc )
    {
        return computeDewPoint( getWaterVaporPartialPressure( currentUtc ) );
    }


};

class InterpolatedStationVmfMeteoData: public StationMeteoData
{
public:

    InterpolatedStationVmfMeteoData(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > meteoDataInterpolator ):
        meteoDataInterpolator_( meteoDataInterpolator ){ }

    ~InterpolatedStationVmfMeteoData( ){ }

    double getTemperature( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 1 );
    }

    double getPressure( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 0 );
    }

    double getWaterVaporPartialPressure( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 2 );
    }

private:

    void updateData( const double currentUtc )
    {
        if( !( currentUtc == currentUtc_ ) )
        {
            currentUtc_ = currentUtc;
            currentData_ = meteoDataInterpolator_->interpolate( currentUtc_ );
        }
    }

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector3d > > meteoDataInterpolator_;

    Eigen::VectorXd currentData_;

    double currentUtc_;
};


class InterpolatedStationDsnMeteoData: public StationMeteoData
{
public:

    InterpolatedStationDsnMeteoData(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > meteoDataInterpolator ):
        meteoDataInterpolator_( meteoDataInterpolator ), currentUtc_( TUDAT_NAN ), currentData_( Eigen::VectorXd::Zero( 5 ) ){ }

    ~InterpolatedStationDsnMeteoData( ){ }


    double getTemperature( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 1 );
    }

    double getPressure( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 2 );
    }

    double getWaterVaporPartialPressure( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 3 );
    }

    virtual double getRelativeHumidity( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 4 );
    }

    virtual double getDewPointTemperature( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( 0 );
    }

private:

    void updateData( const double currentUtc )
    {
        if( !( currentUtc == currentUtc_ ) )
        {
            currentUtc_ = currentUtc;
            currentData_ = meteoDataInterpolator_->interpolate( currentUtc_ );
        }
    }

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > meteoDataInterpolator_;

    double currentUtc_;

    Eigen::VectorXd currentData_;

};


class StationTroposphereData
{
public:

    virtual ~StationTroposphereData( ){ }

    virtual Eigen::Vector2d getZenithDelay( const double currentUtc ) = 0;

    virtual Eigen::Vector2d getMappingFunction( const double currentUtc ) = 0;

    virtual Eigen::Vector4d getGradient( const double currentUtc ) = 0;

};

class InterpolatedStationTroposphereData: public StationTroposphereData
{
public:

    InterpolatedStationTroposphereData(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > troposphereInterpolator,
        const bool dataHasMappingFunction,
        const bool dataHasGradient ): gradientStartIndex_( dataHasMappingFunction ? 4 : 2 ){ }

    ~InterpolatedStationTroposphereData( ){ }

    Eigen::Vector2d getZenithDelay( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_.segment( 0, 2 );
    }

    Eigen::Vector2d getMappingFunction( const double currentUtc )
    {
        if( dataHasMappingFunction_ )
        {
            updateData( currentUtc );
        }
        else
        {
            throw std::runtime_error( "Error when retrieving mapping function, no data set." );
        }
        return currentData_.segment( 2, 2 );
    }

    Eigen::Vector4d getGradient( const double currentUtc )
    {
        if( dataHasGradient_ )
        {
            updateData( currentUtc );
        }
        else
        {
            throw std::runtime_error( "Error when retrieving mapping function, no data set." );
        }
        return currentData_.segment( gradientStartIndex_, 2 );
    }

private:

    void updateData( const double currentUtc )
    {
        if( !( currentUtc == currentUtc_ ) )
        {
            currentUtc_ = currentUtc;
            currentData_ = troposphereInterpolator_->interpolate( currentUtc_ );
        }
    }


    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > troposphereInterpolator_;

    bool dataHasMappingFunction_;

    bool dataHasGradient_;


    int gradientStartIndex_;


    Eigen::VectorXd currentData_;

    double currentUtc_;
};

}  // namespace ground_stations

}  // namespace tudat

#endif  // TUDAT_METEOROLOGICALCONDITIONS_H
