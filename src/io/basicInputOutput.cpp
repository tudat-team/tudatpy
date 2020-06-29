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
 *      The function printStandardScientificNotation() has been implemented to cope with
 *      cross-platform incompatibilities in the printed output of floating-point numbers in
 *      scientific notation.
 *
 */
 
#include <sstream>
#include <string>

#include "tudat/io/basicInputOutput.h"

namespace tudat
{

namespace paths {

    std::string get_resource_path() {
        // TODO: Improve robustness of resource path. Downloading of resources if they are not found should be
        //  added. Users should be able to purge (delete) the resources after arduous use of the software, starting
        //  with no resources.
        return std::string(get_resources_path()).c_str();
    }

    std::string get_tudat_data_path() {
        return std::string(get_resource_path()).c_str();
    }

    std::string get_tudat_path() {
        return std::string(get_resource_path()).c_str();
    }

    std::string get_default_output_path() {
        return std::string(get_resource_path()).c_str();
    }

    std::string getTudatTestDataPath() {
        return std::string(TEST_DATA_FOLDER).c_str();
    }

    std::string getEphemerisDataFilesPath() {
        return std::string(get_ephemeris_path()).c_str();
    }

    std::string getEarthOrientationDataFilesPath() {
        return std::string(get_earth_orientation_path()).c_str();
    }

    std::string getQuadratureDataPath() {
        return std::string(get_quadrature_path()).c_str();
    }

    std::string getSpiceKernelPath() {
        return std::string(get_spice_kernels_path()).c_str();
    }

    std::string getAtmosphereTablesPath() {
        return std::string(get_atmosphere_tables_path()).c_str();
    }

    std::string getGravityModelsPath() {
        return std::string(get_gravity_models_path()).c_str();
    }

    std::string getSpaceWeatherDataPath() {
        return std::string(get_space_weather_path()).c_str();
    }

}
namespace input_output
{

//! Print floating-point number in formatted scientific notation.
std::string printInFormattedScientificNotation( const double floatingPointNumber,
                                                const int basePrecision,
                                                const int exponentWidth )
{
    // Declare string representation of the floating-point number.
    std::string floatingPointNumberString;

    // Extract the decimal part (base) of the floating-point number with specified precision.
    // Default precision is numeric_limits< double >::digits10.
    std::stringstream buffer;
    buffer << std::scientific
           << std::setprecision( basePrecision )
           << std::uppercase
           << floatingPointNumber
           << std::endl;

    // Write base to string.
    buffer >> floatingPointNumberString;

    // Determine the location where the exponent part of the number starts.
    std::string::size_type exponentLocation = floatingPointNumberString.find( 'E' );

    // Set the exponent to zero.
    int exponent = 0;

    // Set flag indicating if exponent is positive or negative to positive.
    bool isExponentPositive = true;

    // Check if the exponent is present.
    if ( exponentLocation != std::string::npos )
    {
        // Extract the exponent to the buffer.
        buffer << floatingPointNumberString.substr( exponentLocation + 1,
                                                    floatingPointNumberString.length( )
                                                    - exponentLocation - 1 )
               << std::endl;

        // Write the exponent part from the buffer.
        buffer >> exponent;

        // Check if the exponent is negative.
        if ( exponent < 0 )
        {
            // If it is negative, set the flag to false, and switch signs of the exponent.
            isExponentPositive = false;
            exponent = -exponent;
        }
    }

    // Extract the exponent with a fixed field width.
    // Default width is 2.
    std::stringstream outputStream;
    outputStream << floatingPointNumberString.substr( 0, exponentLocation + 1 )
                 << ( isExponentPositive ? '+' : '-' )
                 << std::setw( exponentWidth )
                 << std::setfill( '0' )
                 << exponent;

    // Declare output string and write bu ffer output to the string.
    std::string outputString;
    outputStream >> outputString;

    // Return the formatted string representation of the floating-point number.
    return outputString;
}

//! Lists all files in directory.
std::vector< boost::filesystem::path > listAllFilesInDirectory(
    const boost::filesystem::path& directory, const bool isRecurseIntoSubdirectories )
{
    // Declare local variables.
    std::vector < boost::filesystem::path > listOfFileNamesWithPath_;

    if ( boost::filesystem::exists( directory ) )
    {
        boost::filesystem::directory_iterator iteratorPastEndOfDirectory_;

        for ( boost::filesystem::directory_iterator directoryIterator_( directory );
              directoryIterator_ != iteratorPastEndOfDirectory_ ; ++directoryIterator_ )
        {
            if ( boost::filesystem::is_directory( *directoryIterator_ ) )
            {
                if ( isRecurseIntoSubdirectories )
                {
                    listAllFilesInDirectory( *directoryIterator_ );
                }
            }

            else
            {
                listOfFileNamesWithPath_.push_back( directoryIterator_->path( ).filename( ) );
            }
        }
    }

    // Return container of filenames.
    return listOfFileNamesWithPath_;
}

} // namespace input_output
} // namespace tudat
