/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEATMOSPHEREMODEL_H
#define TUDAT_CREATEATMOSPHEREMODEL_H

#include <string>
#include <map>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <algorithm>
#include <sstream>

#include <memory>
#include <boost/date_time/posix_time/time_period.hpp>

#include "body.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/basics/identityElements.h"
#include "tudat/io/comaModelInputOutput.h"
#include <boost/variant.hpp>


namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;


//  List of wind models available in simulations
/*
 *  List of wind models available in simulations. Wind models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum WindModelTypes { empty_wind_model, constant_wind_model, custom_wind_model, coma_wind_model };

//  Class for providing settings for wind model.
/*
 *  Class for providing settings for automatic wind model creation. This class is a
 *  functional (base) class for settings of wind models that require no information in
 *  addition to their type. Wind model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(WindModelSettings.__docstring__)
class WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windModelType Type of wind model that is to be created
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    WindModelSettings( const WindModelTypes windModelType,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                       const bool includeCorotation = true ):
        windModelType_( windModelType ), associatedFrame_( associatedFrame ), includeCorotation_( includeCorotation )
    {
    }

    //  Destructor
    virtual ~WindModelSettings( )
    {
    }

    //  Function to retrieve type of wind model that is to be created
    /*
     * Function to retrieve type of wind model that is to be created
     * \return Type of wind model that is to be created
     */
    WindModelTypes getWindModelType( )
    {
        return windModelType_;
    }

    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

    //  Function to retrieve whether atmospheric co-rotation should be included
    /*
     * Function to retrieve whether atmospheric co-rotation should be included in aerodynamic computations
     * \return Boolean indicating if atmospheric co-rotation is included
     */
    bool getIncludeCorotation( ) const
    {
        return includeCorotation_;
    }

    //  Function to set whether atmospheric co-rotation should be included
    /*
     * Function to set whether atmospheric co-rotation should be included in aerodynamic computations
     * \param includeCorotation Boolean indicating if atmospheric co-rotation should be included
     */
    void setIncludeCorotation( const bool includeCorotation )
    {
        includeCorotation_ = includeCorotation;
    }

protected:
    //  Type of wind model that is to be created
    WindModelTypes windModelType_;

    reference_frames::AerodynamicsReferenceFrames associatedFrame_;

    //  Boolean flag indicating whether atmospheric co-rotation should be included in aerodynamic computations
    bool includeCorotation_;
};

//  Class for empty wind model settings
/*
 * Class for empty wind model settings. This represents a wind model with no physical wind
 * (always returns zero velocity), but allows specification of co-rotation behavior.
 */
class EmptyWindModelSettings : public WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    EmptyWindModelSettings( const bool includeCorotation = true ):
        WindModelSettings( empty_wind_model, reference_frames::vertical_frame, includeCorotation )
    {
    }
};

class ConstantWindModelSettings : public WindModelSettings
{
public:
    ConstantWindModelSettings( const Eigen::Vector3d constantWindVelocity,
                               const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                               const bool includeCorotation = true ):
        WindModelSettings( constant_wind_model, associatedFrame, includeCorotation ), constantWindVelocity_( constantWindVelocity )
    {
    }

    Eigen::Vector3d getConstantWindVelocity( )
    {
        return constantWindVelocity_;
    }

private:
    Eigen::Vector3d constantWindVelocity_;
};

//  Class to define settings for a custom, user-defined, wind model
class CustomWindModelSettings : public WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     * \param associatedFrame Reference frame in which the wind is defined
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    CustomWindModelSettings( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                             const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
                             const bool includeCorotation = true ):
        WindModelSettings( custom_wind_model, associatedFrame, includeCorotation ), windFunction_( windFunction )
    {
    }

    //  Destructor
    ~CustomWindModelSettings( )
    {
    }

    //  Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
     * \return Function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > getWindFunction( )
    {
        return windFunction_;
    }

    //  Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
     * \param windFunction New function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    void setWindFunction( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction )
    {
        windFunction_ = windFunction;
    }

protected:
    //  Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

/**
 * \class ComaWindModelSettings
 * \brief Configuration settings for coma wind models
 *
 * This class handles three separate input datasets for wind velocity components in the modified vertical frame.
 * Each dataset can contain either polynomial coefficients or pre-computed Stokes coefficients.
 *
 * The wind velocity components in these input datasets are defined in a modified vertical frame where:
 *   - X-component: Meridional direction (North, in meridian plane)
 *   - Y-component: Zonal direction (West, completing the right-handed frame)
 *   - Z-component: Radial direction pointing OUTWARD from the nucleus (opposite to standard vertical frame)
 *
 * The constructed ComaWindModel converts the radial component to Tudat's standard vertical_frame
 * convention before returning it to the aerodynamic-angle pipeline.
 */
class ComaWindModelSettings : public WindModelSettings
{
public:
    // Type alias for cleaner code
    using DataVariant = boost::variant< ComaPolyDataset, ComaStokesDataset >;

    /**
     * \brief Constructor from ComaWindDatasetCollection
     * \param datasetCollection Collection containing wind component datasets in modified vertical frame (X=meridional/North, Y=zonal/West, Z=radial outward)
     * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
     * \param associatedFrame Reference frame for the returned wind model vector.
     * \param includeCorotation Boolean indicating whether atmospheric co-rotation should be included
     */
    explicit ComaWindModelSettings( const ComaWindDatasetCollection& datasetCollection,
                                    const int requestedDegree = -1,
                                    const int requestedOrder = -1,
                                    const reference_frames::AerodynamicsReferenceFrames associatedFrame =
                                            reference_frames::vertical_frame,
                                    const bool includeCorotation = true ) :
        WindModelSettings( coma_wind_model, associatedFrame, includeCorotation ),
        requestedDegree_( requestedDegree ),
        requestedOrder_( requestedOrder )
    {
        // Extract datasets from collection based on type
        if (datasetCollection.isPoly())
        {
            xData_ = datasetCollection.getXPolyDataset();
            yData_ = datasetCollection.getYPolyDataset();
            zData_ = datasetCollection.getZPolyDataset();
        }
        else if (datasetCollection.isStokes())
        {
            xData_ = datasetCollection.getXStokesDataset();
            yData_ = datasetCollection.getYStokesDataset();
            zData_ = datasetCollection.getZStokesDataset();
        }
        else
        {
            throw std::runtime_error("ComaWindModelSettings: Invalid dataset collection type");
        }

        validateAndSetDefaults( );
    }

    /**
     * \brief Get the x-component dataset
     */
    const DataVariant& getXData( ) const
    {
        return xData_;
    }

    /**
     * \brief Get the y-component dataset
     */
    const DataVariant& getYData( ) const
    {
        return yData_;
    }

    /**
     * \brief Get the z-component dataset
     */
    const DataVariant& getZData( ) const
    {
        return zData_;
    }

    /**
     * \brief Get requested maximum degree
     */
    int getRequestedDegree( ) const
    {
        return requestedDegree_;
    }

    /**
     * \brief Get requested maximum order
     */
    int getRequestedOrder( ) const
    {
        return requestedOrder_;
    }

    /**
     * \brief Get the effective maximum degree available in the data
     */
    int getAvailableMaxDegree( ) const
    {
        return availableMaxDegree_;
    }

    /**
     * \brief Get the effective maximum order available in the data
     */
    int getAvailableMaxOrder( ) const
    {
        return availableMaxOrder_;
    }

    /**
     * \brief Check if x-component data contains polynomial coefficients
     */
    bool xHasPolyData( ) const
    {
        return xData_.type( ) == typeid( ComaPolyDataset );
    }

    /**
     * \brief Check if y-component data contains polynomial coefficients
     */
    bool yHasPolyData( ) const
    {
        return yData_.type( ) == typeid( ComaPolyDataset );
    }

    /**
     * \brief Check if z-component data contains polynomial coefficients
     */
    bool zHasPolyData( ) const
    {
        return zData_.type( ) == typeid( ComaPolyDataset );
    }

    /**
     * \brief Get x-component polynomial dataset if available
     * \throws std::runtime_error if data is not polynomial type
     */
    const ComaPolyDataset& getXPolyDataset( ) const
    {
        if(auto* p = boost::get< ComaPolyDataset >( &xData_ )) return *p;
        throw std::runtime_error( "X-component data does not contain polynomial data" );
    }

    /**
     * \brief Get y-component polynomial dataset if available
     * \throws std::runtime_error if data is not polynomial type
     */
    const ComaPolyDataset& getYPolyDataset( ) const
    {
        if(auto* p = boost::get< ComaPolyDataset >( &yData_ )) return *p;
        throw std::runtime_error( "Y-component data does not contain polynomial data" );
    }

    /**
     * \brief Get z-component polynomial dataset if available
     * \throws std::runtime_error if data is not polynomial type
     */
    const ComaPolyDataset& getZPolyDataset( ) const
    {
        if(auto* p = boost::get< ComaPolyDataset >( &zData_ )) return *p;
        throw std::runtime_error( "Z-component data does not contain polynomial data" );
    }

    /**
     * \brief Get x-component Stokes dataset if available
     * \throws std::runtime_error if data is not Stokes type
     */
    const ComaStokesDataset& getXStokesDataset( ) const
    {
        if(auto* p = boost::get< ComaStokesDataset >( &xData_ )) return *p;
        throw std::runtime_error( "X-component data does not contain Stokes data" );
    }

    /**
     * \brief Get y-component Stokes dataset if available
     * \throws std::runtime_error if data is not Stokes type
     */
    const ComaStokesDataset& getYStokesDataset( ) const
    {
        if(auto* p = boost::get< ComaStokesDataset >( &yData_ )) return *p;
        throw std::runtime_error( "Y-component data does not contain Stokes data" );
    }

    /**
     * \brief Get z-component Stokes dataset if available
     * \throws std::runtime_error if data is not Stokes type
     */
    const ComaStokesDataset& getZStokesDataset( ) const
    {
        if(auto* p = boost::get< ComaStokesDataset >( &zData_ )) return *p;
        throw std::runtime_error( "Z-component data does not contain Stokes data" );
    }

    /**
     * \brief Check if all components contain polynomial data
     */
    bool hasPolyData( ) const
    {
        return xHasPolyData( ) && yHasPolyData( ) && zHasPolyData( );
    }

    /**
     * \brief Check if all components contain Stokes data
     */
    bool hasStokesData( ) const
    {
        return !xHasPolyData( ) && !yHasPolyData( ) && !zHasPolyData( );
    }

private:
    /**
     * \brief Validate settings and set defaults for degree/order
     */
    void validateAndSetDefaults( )
    {
        // Determine available maxima from all three datasets
        availableMaxDegree_ = std::min( {
            getMaxDegreeFromData( xData_ ),
            getMaxDegreeFromData( yData_ ),
            getMaxDegreeFromData( zData_ )
        } );

        availableMaxOrder_ = std::min( {
            getMaxOrderFromData( xData_ ),
            getMaxOrderFromData( yData_ ),
            getMaxOrderFromData( zData_ )
        } );

        // Set defaults if -1
        if(requestedDegree_ < 0)
        {
            requestedDegree_ = availableMaxDegree_;
        }
        if(requestedOrder_ < 0)
        {
            requestedOrder_ = availableMaxOrder_;
        }

        // Validate requested values don't exceed available
        if(requestedDegree_ > availableMaxDegree_)
        {
            throw std::invalid_argument(
                    "Requested degree " + std::to_string( requestedDegree_ ) +
                    " exceeds available maximum " + std::to_string( availableMaxDegree_ ) );
        }
        if(requestedOrder_ > availableMaxOrder_)
        {
            throw std::invalid_argument(
                    "Requested order " + std::to_string( requestedOrder_ ) +
                    " exceeds available maximum " + std::to_string( availableMaxOrder_ ) );
        }
    }

    /**
     * \brief Get maximum degree from a data variant
     */
    static int getMaxDegreeFromData( const DataVariant& data )
    {
        if(data.type( ) == typeid( ComaPolyDataset ))
        {
            const auto& poly = boost::get< ComaPolyDataset >( data );
            int maxDeg = 0;
            for(std::size_t f = 0; f < poly.getNumFiles( ); ++f)
            {
                maxDeg = std::max( maxDeg, poly.getMaxDegreeSH( f ) );
            }
            return maxDeg;
        }
        else if(data.type( ) == typeid( ComaStokesDataset ))
        {
            const auto& stokes = boost::get< ComaStokesDataset >( data );
            return stokes.nmax( );
        }
        return 0;
    }

    /**
     * \brief Get maximum order from a data variant
     */
    static int getMaxOrderFromData( const DataVariant& data )
    {
        if(data.type( ) == typeid( ComaPolyDataset ))
        {
            const auto& poly = boost::get< ComaPolyDataset >( data );
            int maxOrd = 0;
            for(std::size_t f = 0; f < poly.getNumFiles( ); ++f)
            {
                const auto& indices = poly.getSHDegreeAndOrderIndices( f );
                int fileMaxOrd = indices.row( 1 ).abs( ).maxCoeff( );
                maxOrd = std::max( maxOrd, fileMaxOrd );
            }
            return maxOrd;
        }
        else if(data.type( ) == typeid( ComaStokesDataset ))
        {
            const auto& stokes = boost::get< ComaStokesDataset >( data );
            // For Stokes data, order equals degree in the triangular storage
            return stokes.nmax( );
        }
        return 0;
    }

    // Data members
    DataVariant xData_; // Dataset for x-component wind
    DataVariant yData_; // Dataset for y-component wind
    DataVariant zData_; // Dataset for z-component wind
    int requestedDegree_; // User-requested max degree
    int requestedOrder_; // User-requested max order
    int availableMaxDegree_{ 0 }; // Maximum available in data
    int availableMaxOrder_{ 0 }; // Maximum available in data
};

//  List of atmosphere models available in simulations
/*
 *  List of atmosphere models available in simulations. Atmosphere models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum AtmosphereTypes
{
    exponential_atmosphere,
    custom_constant_temperature_atmosphere,
    tabulated_atmosphere,
    nrlmsise00,
    mars_dtm_atmosphere,
    scaled_atmosphere,
    coma_model,
    mcd_atmosphere,
};

//  Class for providing settings for atmosphere model.
/*
 *  Class for providing settings for automatic atmosphere model creation. This class is a
 *  functional (base) class for settings of atmosphere models that require no information in
 *  addition to their type. Atmosphere model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(AtmosphereSettings.__docstring__)
class AtmosphereSettings
{
public:
    //  Constructor, sets type of atmosphere model.
    /*
     *  Constructor, sets type of atmosphere model. Settings for atmosphere models requiring
     *  additional information should be defined in a derived class.
     *  \param atmosphereType Type of atmosphere model that is to be created.
     */
    AtmosphereSettings( const AtmosphereTypes atmosphereType ):
        atmosphereType_( atmosphereType ),
        windSettings_( std::make_shared< EmptyWindModelSettings >( true ) )
    {
    }

    //  Destructor
    virtual ~AtmosphereSettings( )
    {
    }

    //  Function to return type of atmosphere model that is to be created.
    /*
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( )
    {
        return atmosphereType_;
    }

    //  Function to return settings for the atmosphere's wind model.
    /*
     *  Function to return settings for the atmosphere's wind model.
     *  \return Settings for the atmosphere's wind model.
     */
    std::shared_ptr< WindModelSettings > getWindSettings( )
    {
        return windSettings_;
    }

    //  Function to (re)set settings for the atmosphere's wind model.
    /*
     *  Function to (re)set settings for the atmosphere's wind model.
     *  \param windSettings Settings for the atmosphere's wind model.
     */
    void setWindSettings( const std::shared_ptr< WindModelSettings > windSettings )
    {
        windSettings_ = windSettings;
    }

private:
    //   Type of atmosphere model that is to be created.
    AtmosphereTypes atmosphereType_;

    //  Settings for the atmosphere's wind model.
    std::shared_ptr< WindModelSettings > windSettings_;
};

//  AtmosphereSettings for defining an exponential atmosphere.
//! @get_docstring(ExponentialAtmosphereSettings.__docstring__)
class ExponentialAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor, taking full atmosphere model parameters.
     *  \param densityScaleHeight Scale height for density profile of atmosphere.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at ground level.
     *  \param specificGasConstant Specific gas constant for (constant) atmospheric chemical
     *  composition.
     *  \param ratioOfSpecificHeats Ratio of specific heats for (constant) atmospheric chemical
     *  composition.
     */
    ExponentialAtmosphereSettings( const double densityScaleHeight,
                                   const double constantTemperature,
                                   const double densityAtZeroAltitude,
                                   const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                   const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( exponential_atmosphere ), densityScaleHeight_( densityScaleHeight ),
        constantTemperature_( constantTemperature ), densityAtZeroAltitude_( densityAtZeroAltitude ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        bodyWithPredefinedExponentialAtmosphere_( undefined_body )
    {
    }

    //  Default constructor.
    /*
     *  Default constructor, taking only the name of the body for which to load the predefined
     *  exponential atmosphere model parameters.
     *  \param bodyWithPredefinedExponentialAtmosphere Enumeration denoting the name of the body for which the
     *  predefined atmosphere model is to be loaded.
     */
    ExponentialAtmosphereSettings( const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere ):
        AtmosphereSettings( exponential_atmosphere ), bodyWithPredefinedExponentialAtmosphere_( bodyWithPredefinedExponentialAtmosphere )
    {
        // Check that the body name inserted is available
        switch(bodyWithPredefinedExponentialAtmosphere)
        {
            case earth:
            case mars:
                // all is good
                break;
            default:
                throw std::runtime_error(
                        "Error while creating exponential atmosphere. The body name provided "
                        "does not match any predefined atmosphere model. Available models for: "
                        "Earth, Mars." );
        }
    }

    //  Function to return scale heigh for density profile of atmosphere.
    /*
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( )
    {
        return densityScaleHeight_;
    }

    //  Function to return constant atmospheric temperature.
    /*
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Function to return atmospheric density at ground level.
    /*
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( )
    {
        return densityAtZeroAltitude_;
    }

    //  Function to return specific gas constant for (constant) atmospheric chemical
    /*
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats for (constant) atmospheric chemical
    /*
     *  Function to return ratio of specific heats for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return the name of the body for which to load the predefined
    //  atmosphere model parameters.
    BodiesWithPredefinedExponentialAtmospheres getBodyName( )
    {
        return bodyWithPredefinedExponentialAtmosphere_;
    }

private:
    //  Scale heigh for density profile of atmosphere.
    double densityScaleHeight_;

    //  Constant atmospheric temperature.
    double constantTemperature_;

    //  Atmospheric density at ground level.
    double densityAtZeroAltitude_;

    //  Specific gas constant for (constant) atmospheric chemical
    double specificGasConstant_;

    //  Ratio of specific heats for (constant) atmospheric chemical
    double ratioOfSpecificHeats_;

    //  Enumeration denoting the name of the body for which to load the predefined
    //  atmosphere model.
    BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere_;
};

//  AtmosphereSettings for defining custom constant temperature atmosphere.
class CustomConstantTemperatureAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Typedef for density function.
    typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;

    //  Default constructor.
    /*
     *  Default constructor setting all parameters manually.
     *  \param densityFunction Function to retireve the density at the current altitude,
     *      longitude, latitude and time.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    CustomConstantTemperatureAtmosphereSettings( const DensityFunction& densityFunction,
                                                 const double constantTemperature,
                                                 const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                                 const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunction_( densityFunction ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    {
    }

    //  Constructor.
    /*
     *  Constructor which uses one of the built-in density functions as input.
     *  \param densityFunctionType Enumeration denoting which density function to implement.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param modelSpecificParameters Vector of parameters to be used to set up the density
     *      function. Both meaning and number of parameters depends on the model.
     */
    CustomConstantTemperatureAtmosphereSettings( const AvailableConstantTemperatureAtmosphereModels densityFunctionType,
                                                 const double constantTemperature,
                                                 const double specificGasConstant,
                                                 const double ratioOfSpecificHeats,
                                                 const std::vector< double >& modelSpecificParameters ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunctionType_( densityFunctionType ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ), modelSpecificParameters_( modelSpecificParameters )
    {
    }

    //  Get the function to compute the density at the current conditions.
    /*
     *  Get the function to compute the density at the current conditions.
     *  \return Function to compute the density at the current conditions.
     */
    DensityFunction getDensityFunction( )
    {
        return densityFunction_;
    }

    //  Get the type of function to compute the density at the current conditions.
    /*
     *  Get the type of function to compute the density at the current conditions.
     *  \return Type of function to compute the density at the current conditions.
     */
    AvailableConstantTemperatureAtmosphereModels getDensityFunctionType( )
    {
        return densityFunctionType_;
    }

    //  Get constant temperature.
    /*
     *  Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     *  Kelvin.
     *  \return constantTemperature Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Get specific gas constant.
    /*
     *  Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Get ratio of specific heats.
    /*
     *  Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Get model specific parameters.
    /*
     *  Get model specific parameters.
     *  \return Vector of parameters to be used to set up the density function.
     */
    std::vector< double > getModelSpecificParameters( )
    {
        return modelSpecificParameters_;
    }

private:
    //  Function to compute the density at the current conditions.
    /*
     *  Function to compute the density at the current conditions. Note that the independent
     *  variables need to be altitude, longitude, latitude and time, in this precise order.
     */
    DensityFunction densityFunction_;

    //  Enumeration denoting which density function to implement.
    /*
     *  Enumeration denoting which density function to implement
     */
    AvailableConstantTemperatureAtmosphereModels densityFunctionType_;

    //  Constant temperature.
    /*
     *  The atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    const double constantTemperature_;

    //  Specific gas constant.
    /*
     *  Specific gas constant of the atmosphere, its value is assumed constant, due to the
     *  assumption of constant atmospheric composition.
     */
    const double specificGasConstant_;

    //  Ratio of specific heats at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    const double ratioOfSpecificHeats_;

    //  Vector of parameters to be used to set up the density function.
    /*
     *  Vector of parameters to be used to set up the density function. Both meaning and number of parameters depends on the model.
     */
    std::vector< double > modelSpecificParameters_;
};

//  AtmosphereSettings for defining an NRLMSISE00 atmosphere reading space weather data from a text file.
class NRLMSISE00AtmosphereSettings : public AtmosphereSettings
{
public:
    //  Constructor.
    /*
     *  Constructor.
     *  \param spaceWeatherFile File containing space weather data, as in
     *  https://celestrak.com/SpaceData/sw19571001.txt
     */
    NRLMSISE00AtmosphereSettings( const std::string& spaceWeatherFile,
                                  const bool useStormConditions = false,
                                  const bool useAnomalousOxygen = true ):
        AtmosphereSettings( nrlmsise00 ), spaceWeatherFile_( spaceWeatherFile ), useStormConditions_( useStormConditions ),
        useAnomalousOxygen_( useAnomalousOxygen )
    {
    }

    //  Function to return file containing space weather data.
    /*
     *  Function to return file containing space weather data.
     *  \return Filename containing space weather data.
     */
    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

    //  Function to return geomagnetic activity setting.
    /*
     *  Function to return geomagnetic activity setting.
     *  \return Geomagnetic activity value.
     */
    bool getUseStormConditions( )
    {
        return useStormConditions_;
    }

    bool getUseAnomalousOxygen( )
    {
        return useAnomalousOxygen_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;

    //  Boolean denoting how to deal with geomagnetic activity setting.
    /*
     *  Controls the geomagnetic activity behavior.
     *  If true, it uses full vector of Ap values (switch for geomagnetic activity in nrlmsise to -1).
     *  If false, it only uses daily value of Ap (switch for geomagnetic activity in nrlmsise to 0).
     */
    bool useStormConditions_;

    bool useAnomalousOxygen_;
};

class MarsDtmAtmosphereSettings : public AtmosphereSettings
{
public:
    MarsDtmAtmosphereSettings( const std::string& spaceWeatherFile = "" ):
        AtmosphereSettings( mars_dtm_atmosphere ), spaceWeatherFile_( spaceWeatherFile )
    {
    }

    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;
};

#if TUDAT_BUILD_WITH_MCD_INTERFACE

class McdAtmosphereSettings: public AtmosphereSettings
{
public:

    McdAtmosphereSettings( ):
        AtmosphereSettings( mcd_atmosphere ){ }

};

#endif

//  AtmosphereSettings for defining an atmosphere with tabulated data from file.
// //! @get_docstring(TabulatedAtmosphereSettings.__docstring__)
class TabulatedAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::map< int, std::string >& atmosphereTableFile,
            const std::vector< AtmosphereIndependentVariables >& independentVariablesNames = { altitude_dependent_atmosphere },
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = {},
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ), boundaryHandling_( boundaryHandling ),
        defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with single boundary handling parameters.
    /*
     *  Constructor with single boundary handling parameters. The specifier is assumed to be the same for
     *  each (in)dependent variable.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const double specificGasConstant,
                                 const double ratioOfSpecificHeats,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                specificGasConstant,
                ratioOfSpecificHeats,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor compatible with old version.
    /*
     *  Constructor compatible with old version.
     *  \param atmosphereTableFile File containing atmospheric properties. The file name of the atmosphere table. The file
     *      should contain four columns of data, containing altitude (first column), and the associated density, pressure and
     *      density values in the second, third and fourth columns.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::string& atmosphereTableFile,
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
            const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings( { { 0, atmosphereTableFile } },
                                     { altitude_dependent_atmosphere },
                                     dependentVariablesNames,
                                     specificGasConstant,
                                     ratioOfSpecificHeats,
                                     { boundaryHandling },
                                     std::vector< std::vector< std::pair< double, double > > >(
                                             dependentVariablesNames.size( ),
                                             std::vector< std::pair< double, double > >(
                                                     1,
                                                     std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector).
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< double >& defaultExtrapolationValue ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling )
    {
        // Assign default values
        defaultExtrapolationValue_.resize( dependentVariablesNames.size( ) );
        for(unsigned int i = 0; i < dependentVariablesNames.size( ); i++)
        {
            for(unsigned int j = 0; j < independentVariablesNames.size( ); j++)
            {
                if(boundaryHandling_.at( j ) == interpolators::use_default_value ||
                    boundaryHandling_.at( j ) == interpolators::use_default_value_with_warning)
                {
                    defaultExtrapolationValue_.at( i ).push_back(
                            std::make_pair( defaultExtrapolationValue.at( i ), defaultExtrapolationValue.at( i ) ) );
                }
                else
                {
                    defaultExtrapolationValue_.at( i ).push_back( std::make_pair( IdentityElement::getAdditionIdentity< double >( ),
                                                                                  IdentityElement::getAdditionIdentity< double >( ) ) );
                }
            }
        }
    }

    //  Constructor with no specific gas constant nor ratio of specific heats, and with
    //  single boundary handling parameters.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector). Only one boundary handling parameter is specified, which is then repeated for
     *  dimension.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                1.4,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Map of filenames containing atmospheric properties.
     */
    std::map< int, std::string > getAtmosphereFile( )
    {
        return atmosphereFile_;
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Filename containing atmospheric properties.
     */
    std::string getAtmosphereFile( const unsigned int fileIndex )
    {
        return atmosphereFile_.at( fileIndex );
    }

    //  Function to return independent variables names.
    /*
     *  Function to return independent variables names.
     *  \return Independent variables.
     */
    std::vector< AtmosphereIndependentVariables > getIndependentVariables( )
    {
        return independentVariables_;
    }

    //  Function to return dependent variables names.
    /*
     *  Function to return dependent variables names.
     *  \return Dependent variables.
     */
    std::vector< AtmosphereDependentVariables > getDependentVariables( )
    {
        return dependentVariables_;
    }

    //  Function to return specific gas constant of the atmosphere.
    /*
     *  Function to return specific gas constant of the atmosphere.
     *  \return Specific gas constant of the atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats of the atmosphere.
    /*
     *  Function to return ratio of specific heats of the atmosphere.
     *  \return Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return boundary handling method.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< interpolators::BoundaryInterpolationType > getBoundaryHandling( )
    {
        return boundaryHandling_;
    }

    //  Function to return default extrapolation value.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< std::vector< std::pair< double, double > > > getDefaultExtrapolationValue( )
    {
        return defaultExtrapolationValue_;
    }

private:
    //  File containing atmospheric properties.
    /*
     *  File containing atmospheric properties, file should contain
     *  columns of atmospheric data with at least density, pressure and temperature,
     *  (whose order is specified in dependentVariables), and with at least one
     *  indendent variables.
     */
    std::map< int, std::string > atmosphereFile_;

    //  A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //  A vector of strings containing the names of the variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //  Specific gas constant of the atmosphere.
    /*
     * Specific gas constant of the atmosphere.
     */
    double specificGasConstant_;

    //  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double ratioOfSpecificHeats_;

    //  Behavior of interpolator when independent variable is outside range.
    /*
     *  Behavior of interpolator when independent variable is outside range.
     */
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling_;

    //  Default value to be used for extrapolation.
    /*
     *  Default value to be used for extrapolation.
     */
    std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValue_;
};

class ScaledAtmosphereSettings : public AtmosphereSettings
{
public:
    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                              const double scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( [ = ]( const double ) {
            return scaling;
        } ),
        isScalingAbsolute_( isScalingAbsolute )
    {
    }

    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings >& baseSettings,
                              const std::function< double( const double ) > scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( scaling ), isScalingAbsolute_( isScalingAbsolute )
    {
    }

    std::shared_ptr< AtmosphereSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< double( const double ) > getScaling( )
    {
        return scaling_;
    }

    bool getIsScalingAbsolute( ) const
    {
        return isScalingAbsolute_;
    }

protected:
    std::shared_ptr< AtmosphereSettings > baseSettings_;

    std::function< double( const double ) > scaling_;

    bool isScalingAbsolute_;
};


/**
 * \class ComaSettings
 * \brief Configuration settings for coma atmosphere models
 *
 * This class can be initialized with either polynomial coefficients or
 * pre-computed Stokes coefficients. It provides a unified interface for
 * passing data to the ComaModel while maintaining flexibility in input types.
 */
class ComaSettings final : public AtmosphereSettings
{
public:
    // Type alias for cleaner code
    using DataVariant = boost::variant< ComaPolyDataset, ComaStokesDataset >;

    /**
     * \brief Constructor with polynomial coefficient data
     * \param polyData Pre-loaded polynomial coefficient dataset
     * \param molecularWeight Molecular weight of the gas species [kg/mol]
     * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
     * in aerodynamic computations (default true for backward compatibility).
     */
    explicit ComaSettings( const ComaPolyDataset& polyData,
                           const double molecularWeight,
                           const int requestedDegree = -1,
                           const int requestedOrder = -1,
                           const bool isLog2Data = true ) :
        AtmosphereSettings( coma_model ),
        data_( polyData ),
        molecularWeight_( molecularWeight ),
        requestedDegree_( requestedDegree ),
        requestedOrder_( requestedOrder ),
        isLog2Data_( isLog2Data )
    {
        validateAndSetDefaults( );
    }

    /**
     * \brief Constructor with Stokes coefficient data
     * \param stokesData Pre-computed Stokes coefficient dataset
     * \param molecularWeight Molecular weight of the gas species [kg/mol]
     * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
     * in aerodynamic computations (default true for backward compatibility).
     */
    explicit ComaSettings( const ComaStokesDataset& stokesData,
                           const double molecularWeight,
                           const int requestedDegree = -1,
                           const int requestedOrder = -1,
                           const bool isLog2Data = true ) :
        AtmosphereSettings( coma_model),
        data_( stokesData ),
        molecularWeight_( molecularWeight ),
        requestedDegree_( requestedDegree ),
        requestedOrder_( requestedOrder ),
        isLog2Data_( isLog2Data )
    {
        validateAndSetDefaults( );
    }

    /**
     * \brief Get the underlying data (poly or Stokes coefficients)
     * \return Variant containing either ComaPolyDataset or ComaStokesDataset
     */
    const DataVariant& getData( ) const
    {
        return data_;
    }

    /**
     * \brief Check if settings contain polynomial coefficient data
     */
    bool hasPolyData( ) const
    {
        return data_.type( ) == typeid( ComaPolyDataset );
    }

    /**
     * \brief Check if settings contain Stokes coefficient data
     */
    bool hasStokesData( ) const
    {
        return data_.type( ) == typeid( ComaStokesDataset );
    }

    /**
     * \brief Get polynomial dataset if available
     * \throws std::runtime_error if data is not polynomial type
     */
    const ComaPolyDataset& getPolyDataset( ) const
    {
        if(auto* p = boost::get< ComaPolyDataset >( &data_ )) return *p;
        throw std::runtime_error( "ComaSettings does not contain polynomial data" );
    }

    /**
     * \brief Get Stokes dataset if available
     * \throws std::runtime_error if data is not Stokes type
     */
    const ComaStokesDataset& getStokesDataset( ) const
    {
        if(auto* p = boost::get< ComaStokesDataset >( &data_ )) return *p;
        throw std::runtime_error( "ComaSettings does not contain Stokes data" );
    }

    /**
     * \brief Get requested maximum degree
     */
    int getRequestedDegree( ) const
    {
        return requestedDegree_;
    }

    /**
     * \brief Get requested maximum order
     */
    int getRequestedOrder( ) const
    {
        return requestedOrder_;
    }

    /**
     * \brief Get molecular weight
     */
    double getMolecularWeight( ) const
    {
        return molecularWeight_;
    }

    /**
     * \brief Get the effective maximum degree available in the data
     */
    int getAvailableMaxDegree( ) const
    {
        return availableMaxDegree_;
    }

    /**
     * \brief Get the effective maximum order available in the data
     */
    int getAvailableMaxOrder( ) const
    {
        return availableMaxOrder_;
    }

    /**
     * \brief Add temperature model from polynomial coefficient data
     * \param temperaturePolyData Pre-loaded polynomial coefficient dataset for temperature
     * \param maxDegree Maximum spherical harmonic degree for temperature (-1 for auto)
     * \param maxOrder Maximum spherical harmonic order for temperature (-1 for auto)
     * \param gamma Heat capacity ratio (default 1.33 for water vapor)
     */
    void addTemperatureModel( const ComaPolyDataset& temperaturePolyData,
                              const int maxDegree = -1,
                              const int maxOrder = -1,
                              const double gamma = 1.33 )
    {
        temperatureData_ = temperaturePolyData;
        temperatureMaxDegree_ = maxDegree;
        temperatureMaxOrder_ = maxOrder;
        heatCapacityRatio_ = gamma;
        hasTemperatureModel_ = true;
    }

    /**
     * \brief Add temperature model from Stokes coefficient data
     * \param temperatureStokesData Pre-computed Stokes coefficient dataset for temperature
     * \param maxDegree Maximum spherical harmonic degree for temperature (-1 for auto)
     * \param maxOrder Maximum spherical harmonic order for temperature (-1 for auto)
     * \param gamma Heat capacity ratio (default 1.33 for water vapor)
     */
    void addTemperatureModel( const ComaStokesDataset& temperatureStokesData,
                              const int maxDegree = -1,
                              const int maxOrder = -1,
                              const double gamma = 1.33 )
    {
        temperatureData_ = temperatureStokesData;
        temperatureMaxDegree_ = maxDegree;
        temperatureMaxOrder_ = maxOrder;
        heatCapacityRatio_ = gamma;
        hasTemperatureModel_ = true;
    }

    /**
     * \brief Check if temperature model has been added
     */
    bool hasTemperatureModel( ) const
    {
        return hasTemperatureModel_;
    }

    /**
     * \brief Get temperature data variant
     */
    const DataVariant& getTemperatureData( ) const
    {
        return temperatureData_;
    }

    /**
     * \brief Check if temperature data is polynomial type
     */
    bool hasTemperaturePolyData( ) const
    {
        return hasTemperatureModel_ && temperatureData_.type( ) == typeid( ComaPolyDataset );
    }

    /**
     * \brief Check if temperature data is Stokes type
     */
    bool hasTemperatureStokesData( ) const
    {
        return hasTemperatureModel_ && temperatureData_.type( ) == typeid( ComaStokesDataset );
    }

    /**
     * \brief Get temperature polynomial dataset if available
     * \throws std::runtime_error if temperature data is not polynomial type
     */
    const ComaPolyDataset& getTemperaturePolyDataset( ) const
    {
        if( auto* p = boost::get< ComaPolyDataset >( &temperatureData_ ) ) return *p;
        throw std::runtime_error( "ComaSettings does not contain temperature polynomial data" );
    }

    /**
     * \brief Get temperature Stokes dataset if available
     * \throws std::runtime_error if temperature data is not Stokes type
     */
    const ComaStokesDataset& getTemperatureStokesDataset( ) const
    {
        if( auto* p = boost::get< ComaStokesDataset >( &temperatureData_ ) ) return *p;
        throw std::runtime_error( "ComaSettings does not contain temperature Stokes data" );
    }

    /**
     * \brief Get temperature maximum degree
     */
    int getTemperatureMaxDegree( ) const
    {
        return temperatureMaxDegree_;
    }

    /**
     * \brief Get temperature maximum order
     */
    int getTemperatureMaxOrder( ) const
    {
        return temperatureMaxOrder_;
    }

    /**
     * \brief Get heat capacity ratio
     */
    double getHeatCapacityRatio( ) const
    {
        return heatCapacityRatio_;
    }

    /**
     * \brief Check if coefficient data represents log2-transformed number density
     */
    bool getIsLog2Data( ) const
    {
        return isLog2Data_;
    }

private:
    /**
     * \brief Validate settings and set defaults for degree/order
     */
    void validateAndSetDefaults( )
    {
        // Determine available maxima from data
        if(hasPolyData( ))
        {
            const auto& poly = getPolyDataset( );
            availableMaxDegree_ = determineMaxDegreeFromPoly( poly );
            availableMaxOrder_ = determineMaxOrderFromPoly( poly );
        }
        else if(hasStokesData( ))
        {
            const auto& stokes = getStokesDataset( );
            availableMaxDegree_ = stokes.nmax( );
            // For Stokes data, order equals degree in the triangular storage
            availableMaxOrder_ = stokes.nmax( );
        }

        // Set defaults if -1
        if(requestedDegree_ < 0)
        {
            requestedDegree_ = availableMaxDegree_;
        }
        if(requestedOrder_ < 0)
        {
            requestedOrder_ = availableMaxOrder_;
        }

        // Validate requested values don't exceed available
        if(requestedDegree_ > availableMaxDegree_)
        {
            throw std::invalid_argument(
                    "Requested degree " + std::to_string( requestedDegree_ ) +
                    " exceeds available maximum " + std::to_string( availableMaxDegree_ ) );
        }
        if(requestedOrder_ > availableMaxOrder_)
        {
            throw std::invalid_argument(
                    "Requested order " + std::to_string( requestedOrder_ ) +
                    " exceeds available maximum " + std::to_string( availableMaxOrder_ ) );
        }
    }

    /**
     * \brief Determine maximum degree from polynomial dataset
     */
    static int determineMaxDegreeFromPoly( const ComaPolyDataset& poly )
    {
        int maxDeg = 0;
        for(std::size_t f = 0; f < poly.getNumFiles( ); ++f)
        {
            maxDeg = std::max( maxDeg, poly.getMaxDegreeSH( f ) );
        }
        return maxDeg;
    }

    /**
     * \brief Determine maximum order from polynomial dataset
     */
    static int determineMaxOrderFromPoly( const ComaPolyDataset& poly )
    {
        int maxOrd = 0;
        for(std::size_t f = 0; f < poly.getNumFiles( ); ++f)
        {
            const auto& indices = poly.getSHDegreeAndOrderIndices( f );
            int fileMaxOrd = indices.row( 1 ).abs( ).maxCoeff( );
            maxOrd = std::max( maxOrd, fileMaxOrd );
        }
        return maxOrd;
    }

    // Data members
    DataVariant data_; // Holds either poly or Stokes data
    double molecularWeight_; // Molecular weight of gas species
    int requestedDegree_; // User-requested max degree
    int requestedOrder_; // User-requested max order
    int availableMaxDegree_{ 0 }; // Maximum available in data
    int availableMaxOrder_{ 0 }; // Maximum available in data

    // Temperature model data members
    bool hasTemperatureModel_{ false }; // Flag indicating if temperature model is added
    DataVariant temperatureData_; // Holds either poly or Stokes data for temperature
    int temperatureMaxDegree_{ -1 }; // Maximum degree for temperature model
    int temperatureMaxOrder_{ -1 }; // Maximum order for temperature model
    double heatCapacityRatio_{ 1.33 }; // Heat capacity ratio (gamma)
    bool isLog2Data_{ true }; // Whether coefficients represent log2-transformed number density
};


//! @get_docstring(exponentialAtmosphereSettings,2)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const double densityScaleHeight,
        const double densityAtZeroAltitude,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4)
{
    return std::make_shared< ExponentialAtmosphereSettings >(
            densityScaleHeight,
            constantTemperature,
            densityAtZeroAltitude,
            specificGasConstant,
            ratioOfSpecificHeats);
}

//! @get_docstring(exponentialAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const double densityScaleHeight,
                                                                            const double densityAtZeroAltitude )
{
    return std::make_shared< ExponentialAtmosphereSettings >( densityScaleHeight, TUDAT_NAN, densityAtZeroAltitude, TUDAT_NAN, TUDAT_NAN );
}

//! @get_docstring(exponentialAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const std::string& bodyName )
{
    BodiesWithPredefinedExponentialAtmospheres bodyId;
    if(bodyName == "Earth")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::earth;
    }
    else if(bodyName == "Mars")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::mars;
    }
    else
    {
        throw std::runtime_error(
                "Error while creating exponential atmosphere. The body name provided "
                "does not match any predefined atmosphere model. Available models for: "
                "Earth, Mars." );
    }
    return std::make_shared< ExponentialAtmosphereSettings >( bodyId );
}

//! @get_docstring(nrlmsise00AtmosphereSettings)
inline std::shared_ptr< AtmosphereSettings > nrlmsise00AtmosphereSettings( const std::string dataFile = paths::getSpaceWeatherDataPath( ) +
                                                                                   "/sw19571001.txt",
                                                                           const bool useStormConditions = true,
                                                                           const bool useAnomalousOxygen = true )
{
    return std::make_shared< NRLMSISE00AtmosphereSettings >( dataFile, useStormConditions, useAnomalousOxygen );
}

inline std::shared_ptr< AtmosphereSettings > marsDtmAtmosphereSettings( )
{
    return std::make_shared< MarsDtmAtmosphereSettings >( "" );
}

inline std::shared_ptr< AtmosphereSettings > mcdAtmosphereSettings( )
{
    return std::make_shared< McdAtmosphereSettings >( );
}


typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;
//! @get_docstring(customConstantTemperatureAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const std::function< double( const double ) > densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    DensityFunction fullDensityFunction = [ = ]( const double altitude, const double, const double, const double ) {
        return densityFunction( altitude );
    };
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            fullDensityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(customConstantTemperatureAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const DensityFunction densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            densityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(scaledAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const std::function< double( const double ) > scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const double scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

inline std::shared_ptr< AtmosphereSettings > tabulatedAtmosphereSettings(
        const std::string& atmosphereTableFile,
        const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                       pressure_dependent_atmosphere,
                                                                                       temperature_dependent_atmosphere },
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile,
                                                            dependentVariablesNames,
                                                            specificGasConstant,
                                                            ratioOfSpecificHeats,
                                                            interpolators::throw_exception_at_boundary,
                                                            IdentityElement::getAdditionIdentity< double >( ) );
}


//@get_docstring(ComaSettings,0)
/*!
 * \brief Create coma atmosphere settings from polynomial coefficient data
 * \param polyData Pre-loaded polynomial coefficient dataset
 * \param molecularWeight Molecular weight of the gas species [kg/mol]
 * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
 * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
 * \return Shared pointer to AtmosphereSettings configured for coma model
 */
inline std::shared_ptr< AtmosphereSettings > comaSettings(
        const ComaPolyDataset& polyData,
        const double molecularWeight,
        const int requestedDegree = -1,
        const int requestedOrder = -1,
        const bool isLog2Data = true )
{
    return std::make_shared< ComaSettings >( polyData, molecularWeight, requestedDegree, requestedOrder, isLog2Data );
}

//@get_docstring(ComaSettings,1)
/*!
 * \brief Create coma atmosphere settings from Stokes coefficient data
 * \param stokesData Pre-computed Stokes coefficient dataset
 * \param molecularWeight Molecular weight of the gas species [kg/mol]
 * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
 * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
 * \return Shared pointer to AtmosphereSettings configured for coma model
 */
inline std::shared_ptr< AtmosphereSettings > comaSettings(
        const ComaStokesDataset& stokesData,
        const double molecularWeight,
        const int requestedDegree = -1,
        const int requestedOrder = -1,
        const bool isLog2Data = true )
{
    return std::make_shared< ComaSettings >( stokesData, molecularWeight, requestedDegree, requestedOrder, isLog2Data );
}


//! @get_docstring(emptyWindModelSettings)
inline std::shared_ptr< WindModelSettings > emptyWindModelSettings(
        const bool includeCorotation = true )
{
    return std::make_shared< EmptyWindModelSettings >( includeCorotation );
}

//! @get_docstring(customWindModelSettings)
inline std::shared_ptr< WindModelSettings > customWindModelSettings(
        const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
        const bool includeCorotation = true )
{
    return std::make_shared< CustomWindModelSettings >( windFunction, associatedFrame, includeCorotation );
}

//! @get_docstring(constantWindModelSettings)
inline std::shared_ptr< WindModelSettings > constantWindModelSettings(
        const Eigen::Vector3d constantWindVelocity,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
        const bool includeCorotation = true )
{
    return std::make_shared< ConstantWindModelSettings >( constantWindVelocity, associatedFrame, includeCorotation );
}

//! @get_docstring(comaWindModelSettings)
inline std::shared_ptr< WindModelSettings > comaWindModelSettings(
        const ComaWindDatasetCollection& datasetCollection,
        const int requestedDegree = -1,
        const int requestedOrder = -1,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame,
        const bool includeCorotation = true )
{
    return std::make_shared< ComaWindModelSettings >( datasetCollection, requestedDegree, requestedOrder, associatedFrame, includeCorotation );
}

//  Function to create a wind model.
/*
 *  Function to create a wind model based on model-specific settings for the wind model.
 *  \param windSettings Settings for the wind model that is to be created, defined
 *  a pointer to an object of class (derived from) WindModelSettings.
 *  \param body Name of the body for which the wind model is to be created.
 *  \return Wind model created according to settings in windSettings.
 */
std::shared_ptr< aerodynamics::WindModel > createWindModel( const std::shared_ptr< WindModelSettings > windSettings,
                                                            const std::string& body,
                                                            const std::shared_ptr< aerodynamics::AtmosphereModel >& atmosphereModel,
                                                            const SystemOfBodies& bodies );

//  Function to create an atmosphere model.
/*
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel( const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
                                                                        const std::string& body,
                                                                        const SystemOfBodies& bodies = SystemOfBodies() );
} // namespace simulation_setup
} // namespace tudat

#endif  // TUDAT_CREATEATMOSPHEREMODEL_H
