/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#include "tudat/astro/basic/unitConversions.h"
#include "tudat/io/basicInputOutput.h"

//! Execute examples on tabulated atmosphere.
int main( )
{
    // Using statements to make code more neat
    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::interpolators;
    using namespace tudat::unit_conversions;

    // Example 1: one-dimensional tabulated atmosphere
    std::cout << "Example 1. ----------------------------------------------------------------------- " << std::endl;
    {
        // Create vector of dependent variables
        std::vector< AtmosphereDependentVariables > dependentVariables;
        dependentVariables.push_back( density_dependent_atmosphere );
        dependentVariables.push_back( pressure_dependent_atmosphere );
        dependentVariables.push_back( temperature_dependent_atmosphere );
        dependentVariables.push_back( gas_constant_dependent_atmosphere );
        dependentVariables.push_back( specific_heat_ratio_dependent_atmosphere );
        dependentVariables.push_back( molar_mass_dependent_atmosphere );

        // Set boundary handling to default value
        BoundaryInterpolationType boundaryHandling = use_boundary_value_with_warning;

        // Create a tabulated atmosphere object.
        std::string tabulatedAtmosphereFile = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphere.dat";
        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFile, // path to file
                                                 dependentVariables, // list of dependent variables
                                                 0.0, // gas constant (value does not matter, since it is specified as dependent variable)
                                                 0.0, // specific heat ratio (see above)
                                                 boundaryHandling ); // choice of boundary handling method in case of out-of-range altitude

        // Set altitudes at which density and speed of sound need to be retireved
        // Note that the atmosphere is only defined between 50 and 10000 kilometers altitude
        std::vector< double > vectorOfAltitudes;
        vectorOfAltitudes.push_back( 0.0 ); // 0 km altitude, will return lower boundary (same as output of conditions below)
        vectorOfAltitudes.push_back( 50.0e3 ); // 50 km altitude
        vectorOfAltitudes.push_back( 10.0e6 ); // 10 000 km altitude
        vectorOfAltitudes.push_back( 100.0e6 ); // 100 000 km altitude, will return upper boundary (same as output of conditions above)

        // Retrieve values of atmosphere at some specified conditions
        for ( unsigned int i = 0; i < vectorOfAltitudes.size( ); i++ )
        {
            std::cout << "Altitude: " << vectorOfAltitudes.at( i ) / 1.0e3 << " km. Density: "
                      << tabulatedAtmosphere.getDensity( vectorOfAltitudes.at( i ) ) << " kg/m^3. Speed of sound: "
                      << tabulatedAtmosphere.getSpeedOfSound( vectorOfAltitudes.at( i ) ) << " m/s." << std::endl;
        }
    }

    // Example 2: multi-dimensional tabulated atmosphere with default extrapolation value
    std::cout << std::endl << "Example 2. ----------------------------------------------------------------------- " << std::endl;
    {
        // Create vector of independent variables
        std::vector< AtmosphereIndependentVariables > independentVariables;
        independentVariables.push_back( longitude_dependent_atmosphere );
        independentVariables.push_back( latitude_dependent_atmosphere );
        independentVariables.push_back( altitude_dependent_atmosphere );

        // Create vector of dependent variables
        std::vector< AtmosphereDependentVariables > dependentVariables;
        dependentVariables.push_back( density_dependent_atmosphere );
        dependentVariables.push_back( pressure_dependent_atmosphere );
        dependentVariables.push_back( temperature_dependent_atmosphere );
        dependentVariables.push_back( gas_constant_dependent_atmosphere );
        dependentVariables.push_back( specific_heat_ratio_dependent_atmosphere );

        // Give names of files in same order as dependent variables
        std::map< int, std::string > tabulatedAtmosphereFiles;
        tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
        tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
        tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
        tabulatedAtmosphereFiles[ 3 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
        tabulatedAtmosphereFiles[ 4 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";

        // Set boundary handling to default value for each independent variable
        std::vector< BoundaryInterpolationType > boundaryHandling = std::vector< BoundaryInterpolationType >( 3, use_default_value );

        // Set default value for each independent variable
        std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValues =
                            std::vector< std::vector< std::pair< double, double > > >(
                                5, std::vector< std::pair< double, double > >( 3, std::make_pair( 0.0, 0.0 ) ) );
        // In Tudat, longitude and latitude are always wrapped between 0 and 2 PI, thus there is no need to specify bounds for
        // these two variables. However, if you want to specify the default values for altitude, you also need to do so for these
        // other two (although you can set their boundaryHandling to something like use_boundary_value and whatever you input here
        // will be safely ignored).

        // Overwrite values corresponding to latitude with NaNs
        defaultExtrapolationValues.at( 0 ).at( 1 ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ); // for density
        defaultExtrapolationValues.at( 1 ).at( 1 ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ); // for pressure
        defaultExtrapolationValues.at( 2 ).at( 1 ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ); // for temperature
        defaultExtrapolationValues.at( 3 ).at( 1 ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ); // for gas constant
        defaultExtrapolationValues.at( 4 ).at( 1 ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ); // for specific heat ratio

        // Overwrite values corresponding to altitude with desired numbers
        defaultExtrapolationValues.at( 0 ).at( 2 ) = std::make_pair( 0.02, 0.0 );
        // for density, give 0.02 kg/m^3 if altitude goes below the lower bound (this roughly corresponds to the density at
        // zero altitude), and give 0.0 kg/m^3 if above the largest value
        defaultExtrapolationValues.at( 1 ).at( 2 ) = std::make_pair( 2.35, 0.0 ); // for pressure
        defaultExtrapolationValues.at( 2 ).at( 2 ) = std::make_pair( 161.0, 186.813 ); // for temperature
        defaultExtrapolationValues.at( 3 ).at( 2 ) = std::make_pair( 190.7, 8183.0 ); // for gas constant
        defaultExtrapolationValues.at( 4 ).at( 2 ) = std::make_pair( 1.377, 1.667 ); // for specific heat ratio

        // Create a tabulated atmosphere object.
        TabulatedAtmosphere tabulatedAtmosphere( tabulatedAtmosphereFiles, // path to file
                                                 independentVariables,  // list of independent variables
                                                 dependentVariables, // list of dependent variables
                                                 boundaryHandling, // choise of boundary handling method in case of out-of-range altitude
                                                 defaultExtrapolationValues ); // value to be output in case of out-of-range altitude
        // Here, the constructor that does not require the input of of constant values for the gas constant and specific heat ratio
        // is used, since they are defined as dependent variables

        // Set longitudes at which density and speed of sound need to be retireved
        std::vector< double > vectorOfLongitudes;
        vectorOfLongitudes.push_back( convertDegreesToRadians( 360.0 ) ); // will return 0.0
        vectorOfLongitudes.push_back( convertDegreesToRadians( 135.0 ) ); // within bounds

        // Set longitudes at which density and speed of sound need to be retireved
        std::vector< double > vectorOfLatitudes;
        vectorOfLatitudes.push_back( convertDegreesToRadians( -180.0 ) ); // will return NaN
        vectorOfLatitudes.push_back( convertDegreesToRadians( 45.0 ) ); // within bounds

        // Set altitudes at which density and speed of sound need to be retireved
        // Note that the atmosphere is only defined between 50 and 10000 kilometers altitude
        std::vector< double > vectorOfAltitudes;
        vectorOfAltitudes.push_back( 0.0 ); // will return 0.02
        vectorOfAltitudes.push_back( 100.0e3 ); // within bounds

        // Retrieve values of atmosphere at some specified conditions
        for ( unsigned int i = 0; i < vectorOfLongitudes.size( ); i++ )
        {
            for ( unsigned int j = 0; j < vectorOfLatitudes.size( ); j++ )
            {
                for ( unsigned int k = 0; k < vectorOfAltitudes.size( ); k++ )
                {
                    std::cout << "Longitude: " << convertRadiansToDegrees( vectorOfLongitudes.at( i ) ) << " deg. "
                              << "Latitude: " << convertRadiansToDegrees( vectorOfLatitudes.at( j ) ) << " deg. "
                              << "Altitude: " << vectorOfAltitudes.at( k ) / 1.0e3 << " km. " << std::endl << "    "
                              << "Density: " << tabulatedAtmosphere.getDensity( vectorOfAltitudes.at( k ),
                                                                                vectorOfLongitudes.at( i ),
                                                                                vectorOfLatitudes.at( j ) ) << " kg/m^3. "
                              << "Speed of sound: " << tabulatedAtmosphere.getSpeedOfSound( vectorOfAltitudes.at( k ),
                                                                                            vectorOfLongitudes.at( i ),
                                                                                            vectorOfLatitudes.at( j ) ) << " m/s."
                              << std::endl;
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

