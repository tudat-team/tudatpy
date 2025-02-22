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


/*! Calculate the partial vapor pressure.
 *
 * Calculate the partial vapor pressure according to the Bean and Dutton (1966) model, as described by Estefan and Sovers
 * (1994), Eq. 16.
 *
 * @param relativeHumidity Relative humidity, defined in [0,1]
 * @param temperature Temperature in Kelvin
 * @return Partial vapor pressure in Pa
 */
double calculateBeanAndDuttonWaterVaporPartialPressure( double relativeHumidity, double temperature );

/*!
 * Returns a function that computes the water vapor partial pressure as a function of time, according to the Bean and
 * Dutton (1966) model.
 *
 * @param relativeHumidity Relative humidity as a function of time.
 * @param temperature Temperature as a function of time.
 * @return Water vapor partial pressure as a function of time.
 */
std::function< double( const double time ) > getBeanAndDuttonWaterVaporPartialPressureFunction(
    std::function< double( const double time ) > relativeHumidity,
    std::function< double( const double time ) > temperature );

// Function to compute saturation vapor pressure (Pa) using Tetens' formula
double computeSaturationWaterVaporPressure( const double temperature);

double computeDewPoint( const double relativeHumidity, const double temperature );


enum MeteoDataEntries
{
    temperature_meteo_data,
    pressure_meteo_data,
    water_vapor_pressure_meteo_data,
    relative_humidity_meteo_data,
    dew_point_meteo_data
};

class StationMeteoData
{
public:
    StationMeteoData( const std::map< MeteoDataEntries, int > vectorEntries ):
        vectorEntries_( vectorEntries ),
        currentUtc_( TUDAT_NAN ),
        currentData_( Eigen::VectorXd::Zero( vectorEntries.size( ) ) )
    {
        setAndValidateInput( vectorEntries );
    }

    virtual ~StationMeteoData( ){ }

    double getTemperature( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( vectorEntries_.at( temperature_meteo_data ) );
    }

    double getPressure( const double currentUtc )
    {
        updateData( currentUtc );
        return currentData_( vectorEntries_.at( pressure_meteo_data ) );
    }

    double getWaterVaporPartialPressure( const double currentUtc )
    {
        if( hasVaporPressure_ )
        {
            updateData( currentUtc );
            return currentData_( vectorEntries_.at( water_vapor_pressure_meteo_data ));
        }
        else
        {
            return getRelativeHumidity( currentUtc ) *
                   computeSaturationWaterVaporPressure( getTemperature( currentUtc ));
        }
    }

    virtual double getRelativeHumidity( const double currentUtc )
    {
        if( hasRelativeHumidity_ )
        {
            updateData( currentUtc );
            return currentData_( vectorEntries_.at( relative_humidity_meteo_data ));
        }
        else
        {
            return getWaterVaporPartialPressure( currentUtc ) /
                   computeSaturationWaterVaporPressure( getTemperature( currentUtc ));
        }
    }

    virtual double getDewPointTemperature( const double currentUtc )
    {
        if( hasDewPoint_ )
        {
            updateData( currentUtc );
            return currentData_( vectorEntries_.at( dew_point_meteo_data ));
        }
        else
        {
            return computeDewPoint( getRelativeHumidity( currentUtc ), getTemperature( currentUtc ));
        }
    }

protected:

    virtual void updateData( const double currentUtc ) = 0;

    void setAndValidateInput( const std::map< MeteoDataEntries, int >& vectorEntries )
    {
        if( vectorEntries.count( pressure_meteo_data ) == 0 )
        {
            throw std::runtime_error( "Error, meteo data requires pressure" );
        }
        if( vectorEntries.count( temperature_meteo_data ) == 0 )
        {
            throw std::runtime_error( "Error, meteo data requires temperature" );
        }
        hasVaporPressure_ = ( vectorEntries.count( water_vapor_pressure_meteo_data ) == 0 ) ? false : true;
        hasRelativeHumidity_ = ( vectorEntries.count( relative_humidity_meteo_data ) == 0 ) ? false : true;
        hasDewPoint_ = ( vectorEntries.count( dew_point_meteo_data ) == 0 ) ? false : true;

        if( !hasVaporPressure_ && !hasRelativeHumidity_ )
        {
            throw std::runtime_error( "Error, meteo data requires vapor pressure or relative humidity." );
        }
    }

    std::map< MeteoDataEntries, int > vectorEntries_;

    double currentUtc_;

    Eigen::VectorXd currentData_;

    bool hasVaporPressure_;

    bool hasRelativeHumidity_;

    bool hasDewPoint_;


};


class ContinuousInterpolatedMeteoData: public StationMeteoData
{
public:

    ContinuousInterpolatedMeteoData(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > meteoDataInterpolator,
        const std::map< MeteoDataEntries, int > vectorEntries ):StationMeteoData( vectorEntries ),
        meteoDataInterpolator_( meteoDataInterpolator ){ }

    ~ContinuousInterpolatedMeteoData( ){ }

private:

    void updateData( const double currentUtc );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > meteoDataInterpolator_;

};


class PiecewiseInterpolatedMeteoData: public StationMeteoData
{
public:

    PiecewiseInterpolatedMeteoData(
        const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > > meteoDataInterpolators,
        const std::map< MeteoDataEntries, int > vectorEntries ):StationMeteoData( vectorEntries ),
        meteoDataInterpolators_( meteoDataInterpolators )
    {
        for( unsigned int i = 0; i < meteoDataInterpolators_.size( ); i++ )
        {
            startTimes_.push_back( meteoDataInterpolators.at( i )->getIndependentValues( ).front( ) );
            endTimes_.push_back( meteoDataInterpolators.at( i )->getIndependentValues( ).back( ) );
        }
        lookUpScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( startTimes_ );
    }

    ~PiecewiseInterpolatedMeteoData( ){ }

    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > > getMeteoDataInterpolators( )
    {
        return meteoDataInterpolators_;
    }

    std::vector< double > getStartTimes( )
    {
        return startTimes_;
    }

    std::vector< double > getEndTimes( )
    {
        return endTimes_;
    }

private:

    void updateData( const double currentUtc );

    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > > meteoDataInterpolators_;

    std::vector< double > startTimes_;

    std::vector< double > endTimes_;

    std::shared_ptr< interpolators::LookUpScheme< double > > lookUpScheme_;

    int currentInterpolator_;
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
