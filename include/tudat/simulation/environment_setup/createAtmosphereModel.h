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

#include <memory>

#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/basics/identityElements.h"

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
enum WindModelTypes
{
    constant_wind_model,
    custom_wind_model
};

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
     */
    WindModelSettings(
            const WindModelTypes windModelType,
            const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        windModelType_( windModelType ), associatedFrame_( associatedFrame ){ }

    //  Destructor
    virtual ~WindModelSettings( ){ }

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


protected:

    //  Type of wind model that is to be created
    WindModelTypes windModelType_;

    reference_frames::AerodynamicsReferenceFrames associatedFrame_;

};

class ConstantWindModelSettings: public WindModelSettings
{
public:

    ConstantWindModelSettings(
            const Eigen::Vector3d constantWindVelocity,
            const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModelSettings( constant_wind_model, associatedFrame ),
        constantWindVelocity_( constantWindVelocity ){ }

    Eigen::Vector3d getConstantWindVelocity( )
    {
        return constantWindVelocity_;
    }

private:

    Eigen::Vector3d constantWindVelocity_;

};

//  Class to define settings for a custom, user-defined, wind model
class CustomWindModelSettings: public WindModelSettings
{
public:

    //  Constructor
    /* 
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     */
    CustomWindModelSettings(
            const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
            const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame  ):
        WindModelSettings( custom_wind_model, associatedFrame ), windFunction_( windFunction )
    { }

    //  Destructor
    ~CustomWindModelSettings( ){ }

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
    void setWindFunction(
            const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction )
    {
        windFunction_ = windFunction;
    }

protected:

    //  Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
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
    scaled_atmosphere
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
        atmosphereType_( atmosphereType ){ }

    //  Destructor
    virtual ~AtmosphereSettings( ){ }

    //  Function to return type of atmosphere model that is to be created.
    /* 
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( ){ return atmosphereType_; }

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
class ExponentialAtmosphereSettings: public AtmosphereSettings
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
    ExponentialAtmosphereSettings(
            const double densityScaleHeight,
            const double constantTemperature,
            const double densityAtZeroAltitude,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( exponential_atmosphere ),
        densityScaleHeight_( densityScaleHeight ), constantTemperature_( constantTemperature ),
        densityAtZeroAltitude_( densityAtZeroAltitude ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ), bodyWithPredefinedExponentialAtmosphere_( undefined_body )
    { }

    //  Default constructor.
    /* 
     *  Default constructor, taking only the name of the body for which to load the predefined
     *  exponential atmosphere model parameters.
     *  \param bodyWithPredefinedExponentialAtmosphere Enumeration denoting the name of the body for which the
     *  predefined atmosphere model is to be loaded.
     */
    ExponentialAtmosphereSettings(
            const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere ):
        AtmosphereSettings( exponential_atmosphere ),
        bodyWithPredefinedExponentialAtmosphere_( bodyWithPredefinedExponentialAtmosphere )
    {
        // Check that the body name inserted is available
        switch ( bodyWithPredefinedExponentialAtmosphere )
        {
        case earth:
        case mars:
            // all is good
            break;
        default:
            throw std::runtime_error( "Error while creating exponential atmosphere. The body name provided "
                                      "does not match any predefined atmosphere model. Available models for: "
                                      "Earth, Mars." );
        }
    }

    //  Function to return scale heigh for density profile of atmosphere.
    /* 
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( ){ return densityScaleHeight_; }

    //  Function to return constant atmospheric temperature.
    /* 
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( ){ return constantTemperature_; }

    //  Function to return atmospheric density at ground level.
    /* 
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( ){ return densityAtZeroAltitude_; }

    //  Function to return specific gas constant for (constant) atmospheric chemical
    /* 
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( ){ return specificGasConstant_; }

    //  Function to return ratio of specific heats for (constant) atmospheric chemical
    /* 
     *  Function to return ratio of specific heats for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getRatioOfSpecificHeats( ){ return ratioOfSpecificHeats_; }

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
class CustomConstantTemperatureAtmosphereSettings: public AtmosphereSettings
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
    CustomConstantTemperatureAtmosphereSettings(
            const DensityFunction& densityFunction,
            const double constantTemperature,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 ) :
        AtmosphereSettings( custom_constant_temperature_atmosphere ),
        densityFunction_( densityFunction ),
        constantTemperature_( constantTemperature ),
        specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    { }

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
    CustomConstantTemperatureAtmosphereSettings(
            const AvailableConstantTemperatureAtmosphereModels densityFunctionType,
            const double constantTemperature,
            const double specificGasConstant,
            const double ratioOfSpecificHeats,
            const std::vector< double >& modelSpecificParameters ) :
        AtmosphereSettings( custom_constant_temperature_atmosphere ),
        densityFunctionType_( densityFunctionType ),
        constantTemperature_( constantTemperature ),
        specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        modelSpecificParameters_( modelSpecificParameters )
    { }

    //  Get the function to compute the density at the current conditions.
    /* 
     *  Get the function to compute the density at the current conditions.
     *  \return Function to compute the density at the current conditions.
     */
    DensityFunction getDensityFunction( ) { return densityFunction_; }

    //  Get the type of function to compute the density at the current conditions.
    /* 
     *  Get the type of function to compute the density at the current conditions.
     *  \return Type of function to compute the density at the current conditions.
     */
    AvailableConstantTemperatureAtmosphereModels getDensityFunctionType( ) { return densityFunctionType_; }

    //  Get constant temperature.
    /* 
     *  Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     *  Kelvin.
     *  \return constantTemperature Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( ) { return constantTemperature_; }

    //  Get specific gas constant.
    /* 
     *  Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( ) { return specificGasConstant_; }

    //  Get ratio of specific heats.
    /* 
     *  Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( ) { return ratioOfSpecificHeats_; }

    //  Get model specific parameters.
    /* 
     *  Get model specific parameters.
     *  \return Vector of parameters to be used to set up the density function.
     */
    std::vector< double > getModelSpecificParameters( ) { return modelSpecificParameters_; }

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
class NRLMSISE00AtmosphereSettings: public AtmosphereSettings
{
public:

    //  Constructor.
    /* 
     *  Constructor.
     *  \param spaceWeatherFile File containing space weather data, as in
     *  https://celestrak.com/SpaceData/sw19571001.txt
     */
    NRLMSISE00AtmosphereSettings( const std::string& spaceWeatherFile ):
        AtmosphereSettings( nrlmsise00 ), spaceWeatherFile_( spaceWeatherFile ){ }

    //  Function to return file containing space weather data.
    /* 
     *  Function to return file containing space weather data.
     *  \return Filename containing space weather data.
     */
    std::string getSpaceWeatherFile( ){ return spaceWeatherFile_; }

private:

    //  File containing space weather data.
    /* 
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;
};


//  AtmosphereSettings for defining an atmosphere with tabulated data from file.
// //! @get_docstring(TabulatedAtmosphereSettings.__docstring__)
class TabulatedAtmosphereSettings: public AtmosphereSettings
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
            pressure_dependent_atmosphere, temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = { },
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = { } ) :
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    { }

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
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ) :
        TabulatedAtmosphereSettings( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                                     specificGasConstant, ratioOfSpecificHeats,
                                     std::vector< interpolators::BoundaryInterpolationType >(
                                         independentVariablesNames.size( ), boundaryHandling ),
                                     std::vector< std::vector< std::pair< double, double > > >(
                                         dependentVariablesNames.size( ), std::vector< std::pair< double, double > >(
                                             independentVariablesNames.size( ), std::make_pair(
                                                 defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    { }

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
            pressure_dependent_atmosphere, temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
            const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ) :
        TabulatedAtmosphereSettings( { { 0, atmosphereTableFile } }, { altitude_dependent_atmosphere },
                                     dependentVariablesNames, specificGasConstant,
                                     ratioOfSpecificHeats, { boundaryHandling },
                                     std::vector< std::vector< std::pair< double, double > > >(
                                         dependentVariablesNames.size( ), std::vector< std::pair< double, double > >(
                                             1, std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    { }

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
    TabulatedAtmosphereSettings(
            const std::map< int, std::string >& atmosphereTableFile,
            const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = { } ) :
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    { }

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
                                 const std::vector< double >& defaultExtrapolationValue ) :
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling )
    {
        // Assign default values
        defaultExtrapolationValue_.resize( dependentVariablesNames.size( ) );
        for ( unsigned int i = 0; i < dependentVariablesNames.size( ); i++ )
        {
            for ( unsigned int j = 0; j < independentVariablesNames.size( ); j++ )
            {
                if ( boundaryHandling_.at( j ) == interpolators::use_default_value ||
                     boundaryHandling_.at( j ) == interpolators::use_default_value_with_warning )
                {
                    defaultExtrapolationValue_.at( i ).push_back( std::make_pair( defaultExtrapolationValue.at( i ),
                                                                                  defaultExtrapolationValue.at( i ) ) );
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
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ) :
        TabulatedAtmosphereSettings( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                                     physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                                     std::vector< interpolators::BoundaryInterpolationType >(
                                         independentVariablesNames.size( ), boundaryHandling ),
                                     std::vector< std::vector< std::pair< double, double > > >(
                                         dependentVariablesNames.size( ), std::vector< std::pair< double, double > >(
                                             independentVariablesNames.size( ), std::make_pair(
                                                 defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    { }

    //  Function to return file containing atmospheric properties.
    /* 
     *  Function to return file containing atmospheric properties.
     *  \return Map of filenames containing atmospheric properties.
     */
    std::map< int, std::string > getAtmosphereFile( ){ return atmosphereFile_; }

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
    std::vector< AtmosphereIndependentVariables > getIndependentVariables( ){ return independentVariables_; }

    //  Function to return dependent variables names.
    /* 
     *  Function to return dependent variables names.
     *  \return Dependent variables.
     */
    std::vector< AtmosphereDependentVariables > getDependentVariables( ){ return dependentVariables_; }

    //  Function to return specific gas constant of the atmosphere.
    /* 
     *  Function to return specific gas constant of the atmosphere.
     *  \return Specific gas constant of the atmosphere.
     */
    double getSpecificGasConstant( ){ return specificGasConstant_; }

    //  Function to return ratio of specific heats of the atmosphere.
    /* 
     *  Function to return ratio of specific heats of the atmosphere.
     *  \return Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double getRatioOfSpecificHeats( ){ return ratioOfSpecificHeats_; }

    //  Function to return boundary handling method.
    /* 
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< interpolators::BoundaryInterpolationType > getBoundaryHandling( ){ return boundaryHandling_; }

    //  Function to return default extrapolation value.
    /* 
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< std::vector< std::pair< double, double > > > getDefaultExtrapolationValue( ){ return defaultExtrapolationValue_; }

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



class ScaledAtmosphereSettings: public AtmosphereSettings
{
public:

    ScaledAtmosphereSettings(
            const std::shared_ptr< AtmosphereSettings > baseSettings,
            const double scaling,
            const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ),
        baseSettings_( baseSettings ), scaling_( [=]( const double ){ return scaling; } ), isScalingAbsolute_( isScalingAbsolute ){ }

    ScaledAtmosphereSettings(
            const std::shared_ptr< AtmosphereSettings > baseSettings,
            const std::function< double( const double ) > scaling,
            const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ),
    baseSettings_( baseSettings ), scaling_( scaling ), isScalingAbsolute_( isScalingAbsolute ){ }

    std::shared_ptr< AtmosphereSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< double( const double ) > getScaling( )
    {
        return scaling_;
    }

    bool getIsScalingAbsolute( )
    {
        return isScalingAbsolute_;
    }

protected:

    std::shared_ptr< AtmosphereSettings > baseSettings_;

    std::function< double( const double ) > scaling_;

    bool isScalingAbsolute_;
};


//! @get_docstring(exponentialAtmosphereSettings,2)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const double densityScaleHeight,
        const double densityAtZeroAltitude,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< ExponentialAtmosphereSettings >(
                densityScaleHeight, constantTemperature, densityAtZeroAltitude, specificGasConstant,
                ratioOfSpecificHeats );
}

//! @get_docstring(exponentialAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const double densityScaleHeight,
        const double densityAtZeroAltitude )
{
    return std::make_shared< ExponentialAtmosphereSettings >(
                densityScaleHeight, TUDAT_NAN, densityAtZeroAltitude, TUDAT_NAN, TUDAT_NAN );
}

//! @get_docstring(exponentialAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const std::string& bodyName )
{
    BodiesWithPredefinedExponentialAtmospheres bodyId;
    if( bodyName == "Earth" )
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::earth;
    }
    else if( bodyName == "Mars" )
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::mars;
    }
    else
    {
        throw std::runtime_error( "Error while creating exponential atmosphere. The body name provided "
                                  "does not match any predefined atmosphere model. Available models for: "
                                  "Earth, Mars." );
    }
    return std::make_shared< ExponentialAtmosphereSettings >(
                bodyId );
}

//! @get_docstring(nrlmsise00AtmosphereSettings)
inline std::shared_ptr< AtmosphereSettings > nrlmsise00AtmosphereSettings(
        const std::string dataFile = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt" )
{
    return std::make_shared< NRLMSISE00AtmosphereSettings >( dataFile );
}

typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;
//! @get_docstring(customConstantTemperatureAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const std::function< double( const double ) > densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    DensityFunction fullDensityFunction = [=](const double altitude, const double, const double, const double ){
        return densityFunction( altitude ); };
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
                fullDensityFunction, constantTemperature, specificGasConstant, ratioOfSpecificHeats );
}

//! @get_docstring(customConstantTemperatureAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const DensityFunction densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
                densityFunction, constantTemperature, specificGasConstant, ratioOfSpecificHeats );
}

//! @get_docstring(scaledAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings(
        const std::shared_ptr< AtmosphereSettings > baseSettings,
        const std::function< double( const double ) > scaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings(
        const std::shared_ptr< AtmosphereSettings > baseSettings,
        const double scaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}


//! @get_docstring(customWindModelSettings)
inline std::shared_ptr< WindModelSettings > customWindModelSettings(
        const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< CustomWindModelSettings >( windFunction, associatedFrame );
}

//! @get_docstring(constantWindModelSettings)
inline std::shared_ptr< WindModelSettings > constantWindModelSettings(
        const Eigen::Vector3d constantWindVelocity,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< ConstantWindModelSettings >( constantWindVelocity, associatedFrame );
}

//  Function to create a wind model.
/* 
 *  Function to create a wind model based on model-specific settings for the wind model.
 *  \param windSettings Settings for the wind model that is to be created, defined
 *  a pointer to an object of class (derived from) WindModelSettings.
 *  \param body Name of the body for which the wind model is to be created.
 *  \return Wind model created according to settings in windSettings.
 */
std::shared_ptr< aerodynamics::WindModel > createWindModel(
        const std::shared_ptr< WindModelSettings > windSettings,
        const std::string& body);

//  Function to create an atmosphere model.
/* 
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel(
        const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
        const std::string& body );

} // namespace simulation_setup

} // namespace tudat


#endif // TUDAT_CREATEATMOSPHEREMODEL_H
