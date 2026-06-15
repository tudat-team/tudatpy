/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_MARSCLIMATEDATABASECLIMATEMODEL_H
#define TUDAT_MARSCLIMATEDATABASECLIMATEMODEL_H

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include "tudat/astro/basic_astro/climateModel.h"

#if defined( TUDAT_MCD_USE_LLVM_FLANG_SYMBOLS )
#define TUDAT_MCD_FORTRAN_CALL_MCD _QMmcdPcall_mcd
#else
#define TUDAT_MCD_FORTRAN_CALL_MCD __mcd_MOD_call_mcd
#endif

extern "C" {
void TUDAT_MCD_FORTRAN_CALL_MCD( int* zkey,
                                 float* xz,
                                 float* xlon,
                                 float* xlat,
                                 int* hireskey,
                                 int* datekey,
                                 double* xdate,
                                 float* localtime,
                                 const char* dset,
                                 int* scena,
                                 int* perturkey,
                                 float* seedin,
                                 float* gwlength,
                                 int* extvarkeys,
                                 float* pres,
                                 float* dens,
                                 float* temp,
                                 float* zonwind,
                                 float* merwind,
                                 float* meanvar,
                                 float* extvar,
                                 float* seedout,
                                 int* ier,
                                 std::size_t dsetLength );
}

namespace tudat
{

namespace mcd_interface
{

using McdCacheKey = std::tuple< int, double, double, double, double >;

enum class MeanVar {
    mean_atmospheric_pressure = 0,
    mean_atmospheric_density = 1,
    mean_atmospheric_temperature = 2,
    mean_zonal_wind = 3,
    mean_meridional_wind = 4
};

enum class ExtVar {
    radial_distance_from_planet_center = 0,
    altitude_above_areoid = 1,
    altitude_above_local_surface = 2,
    orographic_height = 3,
    orographic_height_from_gcm = 4,
    local_slope_inclination = 5,
    local_slope_orientation = 6,
    distance_sun_mars = 7,
    solar_longitude = 8,
    local_true_solar_time = 9,
    local_mean_time = 10,
    universal_solar_time = 11,
    solar_zenith_angle = 12,
    surface_temperature = 13,
    surface_pressure = 14,
    surface_pressure_from_gcm = 15,
    potential_temperature = 16,
    vertical_wind_component = 17,
    zonal_slope_wind_component = 18,
    meridional_slope_wind_component = 19,
    surface_pressure_variations_rms = 20,
    surface_temperature_variations_rms = 21,
    atmospheric_pressure_variations_rms = 22,
    density_variations_rms = 23,
    atmospheric_temperature_variations_rms = 24,
    zonal_wind_rms = 25,
    meridional_wind_rms = 26,
    vertical_wind_rms = 27,
    incident_solar_flux_on_toa = 28,
    reflected_solar_flux_to_space = 29,
    incident_solar_flux_on_horizontal_surface = 30,
    incident_solar_flux_on_local_slope = 31,
    reflected_solar_flux_on_horizontal_surface = 32,
    thermal_ir_flux_to_space = 33,
    thermal_ir_flux_on_surface = 34,
    surface_roughness_from_gcm = 35,
    surface_thermal_inertia_from_gcm = 36,
    surface_bare_ground_albedo_from_gcm = 37,
    monthly_mean_dust_column_visible_above_surface = 38,
    daily_mean_dust_column_visible_above_surface = 39,
    dust_mass_mixing_ratio = 40,
    dust_effective_radius = 41,
    daily_mean_dust_deposition_rate = 42,
    monthly_mean_surface_co2_ice_layer = 43,
    monthly_mean_surface_h2o_ice_layer = 44,
    perennial_surface_water_ice_from_gcm = 45,
    water_vapor_column = 46,
    water_vapor_volume_mixing_ratio = 47,
    water_ice_column = 48,
    water_ice_mixing_ratio = 49,
    water_ice_effective_ratio = 50,
    convective_planetary_boundary_layer_height = 51,
    maximum_upward_convective_wind_within_pbl = 52,
    maximum_downward_convective_wind_within_pbl = 53,
    convective_vertical_wind_variance = 54,
    convective_eddy_vertical_heat_flux = 55,
    surface_wind_stress = 56,
    surface_sensible_heat_flux = 57,
    air_specific_heat_capacity = 58,
    ratio_of_specific_heats = 59,
    reduced_molecular_gas_constant = 60,
    air_viscosity_estimation = 61,
    scale_height = 62,
    co2_volume_mixing_ratio = 63,
    n2_volume_mixing_ratio = 64,
    ar_volume_mixing_ratio = 65,
    co_volume_mixing_ratio = 66,
    o_volume_mixing_ratio = 67,
    o2_volume_mixing_ratio = 68,
    o3_volume_mixing_ratio = 69,
    h_volume_mixing_ratio = 70,
    h2_volume_mixing_ratio = 71,
    he_volume_mixing_ratio = 72,
    co2_column = 73,
    n2_column = 74,
    ar_column = 75,
    co_column = 76,
    o_column = 77,
    o2_column = 78,
    o3_column = 79,
    h_column = 80,
    h2_column = 81,
    he_column = 82,
    electron_number_density = 83,
    total_electric_content = 84,

};

struct McdCache {
    McdCache( ) = default;

    McdCache( double density, double pressure, double temperature, double zonalWind, double meridionalWind ):
        density_( density ), pressure_( pressure ), temperature_( temperature ), zonalWind_( zonalWind ),
        meridionalWind_( meridionalWind ) {};

    //! Atmospheric density (kg/m^3)
    double density_;

    //! Atmospheric pressure (Pa)
    double pressure_;

    //! Atmospheric temperature (K)
    double temperature_;

    //! Zonal wind (m/s)
    double zonalWind_;

    //! Meridional wind (m/s)
    double meridionalWind_;

    //! Mean variables
    std::vector< double > meanVariables_;

    //! Extra variables
    std::vector< double > extraVariables_;
};

class MarsClimateDatabaseClimateModel : public environment::ClimateModel
{
public:
    //! Constructor
    /*!
     * Constructor for MCD climate model.
     * \param mcdDataPath Path to MCD data files (default: "" = use compile-time default)
     * \param dustScenario Dust and solar EUV scenario (1-8 or 24-35, default: 1)
     * \param perturbationKey Perturbation type (0-5, default: 0 = none)
     * \param perturbationSeed Random seed or scaling factor (default: 0.0)
     * \param gravityWaveLength Gravity wave wavelength in meters (default: 0.0 = use MCD default)
     * \param highResolutionMode High resolution topography flag (0 or 1, default: 0)
     */
    MarsClimateDatabaseClimateModel( const std::string& mcdDataPath = "",
                                     const int dustScenario = 1,
                                     const int perturbationKey = 0,
                                     const double perturbationSeed = 0.0,
                                     const double gravityWaveLength = 0.0,
                                     const int highResolutionMode = 0 );

    //! Destructor
    ~MarsClimateDatabaseClimateModel( ) override = default;

    void updateCache( const double verticalCoordinate, const double longitude, const double latitude, const double time );

    int* getExtraVariableKeys( )
    {
        return extraVariableKeys_;
    }

    void setZkey( int zkey );

    std::shared_ptr< McdCache > getCache( const double verticalCoordinate,
                                          const double longitude,
                                          const double latitude,
                                          const double time )
    {
        updateCache( verticalCoordinate, longitude, latitude, time );
        return mcdCache_.at( getCacheKey( verticalCoordinate, longitude, latitude, time ) );
    }

    double getMeanVariable( const MeanVar variable,
                            const double verticalCoordinate,
                            const double longitude,
                            const double latitude,
                            const double time )
    {
        updateCache( verticalCoordinate, longitude, latitude, time );
        return mcdCache_.at( getCacheKey( verticalCoordinate, longitude, latitude, time ) )
                ->meanVariables_[ static_cast< std::size_t >( variable ) ];
    }

    double getExtraVariable( const ExtVar variable,
                             const double verticalCoordinate,
                             const double longitude,
                             const double latitude,
                             const double time )
    {
        updateCache( verticalCoordinate, longitude, latitude, time );
        return mcdCache_.at( getCacheKey( verticalCoordinate, longitude, latitude, time ) )
                ->extraVariables_[ static_cast< std::size_t >( variable ) ];
    }

    void addExtraVariableKeys( std::vector< mcd_interface::ExtVar > requiredExtraVariables );

private:
    McdCacheKey getCacheKey( const double verticalCoordinate, const double longitude, const double latitude, const double time ) const
    {
        return { zkey_, verticalCoordinate, longitude, latitude, time };
    }

    int zkey_;

    //! Path to MCD data files
    std::string mcdDataPath_;

    //! Dust and solar EUV scenario (1-8 or 24-35)
    int dustScenario_;

    //! Perturbation type
    int perturbationKey_;

    //! Perturbation seed
    double perturbationSeed_;

    //! Gravity wave wavelength
    double gravityWaveLength_;

    //! High resolution mode flag
    int highResolutionMode_;

    //! Atmospheric density (kg/m^3)
    float density_;

    //! Atmospheric pressure (Pa)
    float pressure_;

    //! Atmospheric temperature (K)
    float temperature_;

    //! Zonal wind (m/s)
    float zonalWind_;

    //! Meridional wind (m/s)
    float meridionalWind_;

    //! Keys to extract external variables
    int extraVariableKeys_[ 100 ] = { 0 };

    //! Mean variables
    float meanVariables_[ 5 ] = { 0 };

    //! Extra variables
    float extraVariables_[ 100 ] = { 0 };

    std::map< McdCacheKey, std::shared_ptr< McdCache > > mcdCache_;
};

}  // namespace mcd_interface

}  // namespace tudat

#endif
