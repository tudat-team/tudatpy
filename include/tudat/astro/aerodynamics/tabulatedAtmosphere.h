/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The provided USSA1976 table file, generated with the pascal file, has a small error which
 *      can be observed at the pressure at sea level. This in 101320 in the file but should be
 *      101325. If this error is not acceptable, another table file should be used.
 *
 */

#ifndef TUDAT_TABULATED_ATMOSPHERE_H
#define TUDAT_TABULATED_ATMOSPHERE_H

#include <string>

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/utilityMacros.h"

#include "tudat/astro/aerodynamics/standardAtmosphere.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/io/tabulatedAtmosphereReader.h"

namespace tudat
{

namespace aerodynamics
{

//! Check uniqueness of input.
/*!
 *  Function to check uniqueness of input for (in)dependent variables. The function works by checking that there
 *  are no duplicates, after having sorted the vector. For this reason, a local copy of the input vector is taken.
 *  \tparam VariableType Type belonging to the input vector.
 *  \param variables Vector of variables which needs to be checked.
 */
template< typename VariableType >
void checkVariableUniqueness( std::vector< VariableType > variables );

//! Tabulated atmosphere class.
/*!
 *  Tabulated atmospheres class, for example US1976. The default path from which the files are
 *  obtained is: /external/AtmosphereTables
 */
class TabulatedAtmosphere : public AtmosphereModel
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
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
    TabulatedAtmosphere(
            const std::map< int, std::string >& atmosphereTableFile,
            const std::vector< AtmosphereIndependentVariables >& independentVariablesNames = { altitude_dependent_atmosphere },
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
            pressure_dependent_atmosphere, temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = { },
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = { } ) :
        atmosphereTableFile_( atmosphereTableFile ), independentVariables_( independentVariablesNames ),
        dependentVariables_( dependentVariablesNames ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ), boundaryHandling_( boundaryHandling ),
        defaultExtrapolationValue_( defaultExtrapolationValue )
    {
        // Set default dependent variables
        dependentVariablesDependency_ = std::vector< bool >( 6, false ); // only 6 dependent variables supported
        dependentVariableIndices_ = std::vector< unsigned int >( 6, 0 ); // only 6 dependent variables supported

        // Initialize atmosphere
        createAtmosphereInterpolators( );
        independentVariableData_.resize( numberOfIndependentVariables_ );
    }

    //! Constructor with default gas constant and specific heat ratio.
    /*!
     *  Constructor with default gas constant and specific heat ratio.
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
    TabulatedAtmosphere( const std::map< int, std::string >& atmosphereTableFile,
                         const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                         const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                         const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                         const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue ) :
        TabulatedAtmosphere( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                             physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4, boundaryHandling, defaultExtrapolationValue )
    { }

    //! Constructor compatible with old version.
    /*!
     *  Constructor compatible with old version.
     *  \param atmosphereTableFile File containing atmospheric properties.
     *      The file name of the atmosphere table. The file should contain four columns of data,
     *      containing altitude (first column), and the associated density, pressure and density values
     *      in the second, third and fourth columns.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphere(
            const std::string& atmosphereTableFile,
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
            pressure_dependent_atmosphere, temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
            const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ) :
        TabulatedAtmosphere( { { 0, atmosphereTableFile } }, { altitude_dependent_atmosphere },
                             dependentVariablesNames, specificGasConstant, ratioOfSpecificHeats, { boundaryHandling },
                             std::vector< std::vector< std::pair< double, double > > >(
                                 dependentVariablesNames.size( ), { { defaultExtrapolationValue, defaultExtrapolationValue } } ) )
    { }

    //! Constructor.
    /*!
     *  Constructor.
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
    TabulatedAtmosphere( const std::map< int, std::string >& atmosphereTableFile,
                         const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                         const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                         const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                         const std::vector< double >& defaultExtrapolationValue ) :
        atmosphereTableFile_( atmosphereTableFile ), independentVariables_( independentVariablesNames ),
        dependentVariables_( dependentVariablesNames ), specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ),
        ratioOfSpecificHeats_( 1.4 ), boundaryHandling_( boundaryHandling )
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

        // Set default dependent variables
        dependentVariablesDependency_ = std::vector< bool >( 6, false ); // only 6 dependent variables supported
        dependentVariableIndices_ = std::vector< unsigned int >( 6, 0 ); // only 6 dependent variables supported

        // Initialize atmosphere
        createAtmosphereInterpolators( );
        independentVariableData_.resize( numberOfIndependentVariables_ );
    }

    //! Destructor
    ~TabulatedAtmosphere( ){ }

    //! Get atmosphere table file name.
    /*!
     *  Returns atmosphere table file name.
     *  \return The atmosphere table file.
     */
    std::map< int, std::string > getAtmosphereTableFile( ) { return atmosphereTableFile_; }

    //! Get local density.
    /*!
     *  Returns the local density parameter of the atmosphere in kg per meter^3, at the specified conditions.
     *  \param altitude Altitude at which density is to be computed.
     *  \param longitude Longitude at which density is to be computed.
     *  \param latitude Latitude at which density is to be computed.
     *  \param time Time at which density is to be computed.
     *  \return Atmospheric density at specified conditions.
     */
    double getDensity( const double altitude, const double longitude = 0.0,
                       const double latitude = 0.0, const double time = 0.0 )
    {
        // Get list of independent variables
        for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
        {
            switch ( independentVariables_.at( i ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData_[ i ] = altitude;
                break;
            case longitude_dependent_atmosphere:
                independentVariableData_[ i ] = longitude;
                break;
            case latitude_dependent_atmosphere:
                independentVariableData_[ i ] = latitude;
                break;
            case time_dependent_atmosphere:
                independentVariableData_[ i ] = time;
                break;
            }
        }

        // Give output
        return interpolatorForDensity_->interpolate( independentVariableData_ );
    }

    //! Get local pressure.
    /*!
     *  Returns the local pressure of the atmosphere in Newton per meter^2, at the specified conditions.
     *  \param altitude Altitude  at which pressure is to be computed.
     *  \param longitude Longitude at which pressure is to be computed.
     *  \param latitude Latitude at which pressure is to be computed.
     *  \param time Time at which pressure is to be computed.
     *  \return Atmospheric pressure at specified conditions.
     */
    double getPressure( const double altitude, const double longitude = 0.0,
                        const double latitude = 0.0, const double time = 0.0 )
    {
        // Get list of independent variables
        std::vector< double > independentVariableData;
        for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
        {
            switch ( independentVariables_.at( i ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData.push_back( altitude );
                break;
            case longitude_dependent_atmosphere:
                independentVariableData.push_back( longitude );
                break;
            case latitude_dependent_atmosphere:
                independentVariableData.push_back( latitude );
                break;
            case time_dependent_atmosphere:
                independentVariableData.push_back( time );
                break;
            }
        }

        // Give output
        return interpolatorForPressure_->interpolate( independentVariableData );
    }

    //! Get local temperature.
    /*!
     *  Returns the local temperature of the atmosphere in Kelvin, at the specified conditions.
     *  \param altitude Altitude at which temperature is to be computed
     *  \param longitude Longitude at which temperature is to be computed.
     *  \param latitude Latitude at which temperature is to be computed.
     *  \param time Time at which temperature is to be computed.
     *  \return constantTemperature Atmospheric temperature at specified conditions.
     */
    double getTemperature( const double altitude, const double longitude = 0.0,
                           const double latitude = 0.0, const double time = 0.0 )
    {
        // Get list of independent variables
        std::vector< double > independentVariableData;
        for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
        {
            switch ( independentVariables_.at( i ) )
            {
            case altitude_dependent_atmosphere:
                independentVariableData.push_back( altitude );
                break;
            case longitude_dependent_atmosphere:
                independentVariableData.push_back( longitude );
                break;
            case latitude_dependent_atmosphere:
                independentVariableData.push_back( latitude );
                break;
            case time_dependent_atmosphere:
                independentVariableData.push_back( time );
                break;
            }
        }

        // Give output
        return interpolatorForTemperature_->interpolate( independentVariableData );
    }

    //! Get specific gas constant.
    /*!
     *  Returns the specific gas constant of the atmosphere in J/(kg K), at the specified conditions.
     *  \param altitude Altitude at which specific gas constant is to be computed.
     *  \param longitude Longitude at which specific gas constant is to be computed.
     *  \param latitude Latitude at which specific gas constant is to be computed.
     *  \param time Time at which specific gas constant is to be computed.
     *  \return specificGasConstant Specific gas constant at specified conditions.
     */
    double getSpecificGasConstant( const double altitude, const double longitude = 0.0,
                                   const double latitude = 0.0, const double time = 0.0 )
    {
        if ( dependentVariablesDependency_.at( gas_constant_dependent_atmosphere ) )
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return interpolatorForGasConstant_->interpolate( independentVariableData );
        }
        else
        {
            return specificGasConstant_;
        }
    }

    //! Get ratio of specific heats.
    /*!
     *  Returns the ratio of specific heats of the atmosphere at the specified conditions.
     *  \param altitude Altitude at which ratio of specific heats is to be computed
     *  \param longitude Longitude at which ratio of specific heats is to be computed.
     *  \param latitude Latitude at which ratio of specific heats is to be computed.
     *  \param time Time at which ratio of specific heats is to be computed.
     *  \return Ratio of specific heats at specified conditions.
     */
    double getRatioOfSpecificHeats( const double altitude, const double longitude = 0.0,
                                    const double latitude = 0.0, const double time = 0.0 )
    {
        if ( dependentVariablesDependency_.at( specific_heat_ratio_dependent_atmosphere ) )
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return interpolatorForSpecificHeatRatio_->interpolate( independentVariableData );
        }
        else
        {
            return ratioOfSpecificHeats_;
        }
    }

    //! Get molar mass.
    /*!
     *  Returns the molar mass of the atmosphere in kilograms per mole, at the specified conditions.
     *  \param altitude Altitude at which molar mass is to be computed
     *  \param longitude Longitude at which molar mass is to be computed.
     *  \param latitude Latitude at which molar mass is to be computed.
     *  \param time Time at which molar mass is to be computed.
     *  \return Molar mass at specified conditions.
     */
    double getMolarMass( const double altitude, const double longitude = 0.0,
                         const double latitude = 0.0, const double time = 0.0 )
    {
        if ( dependentVariablesDependency_.at( molar_mass_dependent_atmosphere ) )
        {
            // Get list of independent variables
            std::vector< double > independentVariableData;
            for ( unsigned int i = 0; i < numberOfIndependentVariables_; i++ )
            {
                switch ( independentVariables_.at( i ) )
                {
                case altitude_dependent_atmosphere:
                    independentVariableData.push_back( altitude );
                    break;
                case longitude_dependent_atmosphere:
                    independentVariableData.push_back( longitude );
                    break;
                case latitude_dependent_atmosphere:
                    independentVariableData.push_back( latitude );
                    break;
                case time_dependent_atmosphere:
                    independentVariableData.push_back( time );
                    break;
                }
            }

            // Give output
            return interpolatorForMolarMass_->interpolate( independentVariableData );
        }
        else
        {
            throw std::runtime_error( "Error in tabulated atmosphere. The molar mass needs to be specified in the atmosphere "
                                      "table files." );
        }
    }

    //! Get local speed of sound in the atmosphere.
    /*!
     *  Returns the speed of sound in the atmosphere in m/s.
     *  \param altitude Altitude at which speed of sound is to be computed.
     *  \param longitude Longitude at which speed of sound is to be computed.
     *  \param latitude Latitude at which speed of sound is to be computed.
     *  \param time Time at which speed of sound is to be computed.
     *  \return Atmospheric speed of sound at specified conditions.
     */
    double getSpeedOfSound( const double altitude, const double longitude = 0.0,
                            const double latitude = 0.0, const double time = 0.0 )
    {
        return computeSpeedOfSound( getTemperature( altitude, longitude, latitude, time ),
                                    getSpecificGasConstant( altitude, longitude, latitude, time ),
                                    getRatioOfSpecificHeats( altitude, longitude, latitude, time ) );
    }

protected:

private:

    //! Function to create the interpolators based on the tabulated atmosphere files.
    /*!
     *  Function to create the interpolators based on the tabulated atmosphere files, and the provided interpolation settings. This
     *  function also checks the compatibility of (in)dependent variables and determines which (in)dependent variables are in use.
     */
    void createAtmosphereInterpolators( );

    //! Create interpolators for specified dependent variables, taking into consideration the number
    //! of independent variables (which is greater than one).
    /*!
     *  Create interpolators for specified dependent variables, taking into consideration the variable
     *  size of independent variables (which is greater than one).
     *  \tparam Number of independent variables to be used by the interpolator.
     */
    template< unsigned int NumberOfIndependentVariables >
    void createMultiDimensionalAtmosphereInterpolators( );

    //! The file name of the atmosphere table.
    /*!
     *  The file name of the atmosphere table. The file should contain four columns of data,
     *  containing altitude (first column), and the associated density, pressure and density values
     *  in the second, third and fourth columns.
     */
    std::map< int, std::string > atmosphereTableFile_;

    //! A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*!
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //! Vector containing the independent variables.
    std::vector< std::vector< double > > independentVariablesData_;

    //! Integer specifying number of independent variables.
    unsigned int numberOfIndependentVariables_;

    //! A vector of strings containing the names of the variables contained in the atmosphere file.
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //! Vector of booleans that determines if the atmosphere file contains dentity, pressure, temperature,
    //! gas constant and/or ratio of specific heats.
    std::vector< bool > dependentVariablesDependency_;

    //! Vector of integers that specifies the order of dentity, pressure, temperature, gas constant and
    //! ratio of specific heats are located.
    std::vector< unsigned int > dependentVariableIndices_;

    //! Specific gas constant of the atmosphere.
    double specificGasConstant_;

    //! Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    double ratioOfSpecificHeats_;

    //! Interpolation for density. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForDensity_;

    //! Interpolation for pressure. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForPressure_;

    //! Interpolation for temperature. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForTemperature_;

    //! Interpolation for specific gas constant. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForGasConstant_;

    //! Interpolation for ratio of specific heats. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForSpecificHeatRatio_;

    //! Interpolation for molar mass. Note that type of interpolator depends on number of independent variables specified.
    std::shared_ptr< interpolators::Interpolator< double, double > > interpolatorForMolarMass_;

    //! Behavior of interpolator when independent variable is outside range.
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling_;

    //! Default values to be used for extrapolation.
    /*!
     *  Default values to be used for extrapolation. The structure of the vector is as follows:
     *      - outer vector: one element for each dependent variable (density, pressure, temperature, ect.)
     *      - inner vector: one element for each independent variable (latitude, longitude and altitude)
     *      - pair: first element in case independent variable is below its lower definition limit, second element
     *          in case it is above the upper definition limit
     */
    std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValue_;

    std::vector< double > independentVariableData_;

};

//! Typedef for shared-pointer to TabulatedAtmosphere object.
typedef std::shared_ptr< TabulatedAtmosphere > TabulatedAtmospherePointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_TABULATED_ATMOSPHERE_H
