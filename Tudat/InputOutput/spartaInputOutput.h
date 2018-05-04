/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPARTA_INPUT_OUTPUT_H
#define TUDAT_SPARTA_INPUT_OUTPUT_H

#include <string>
#include <vector>

#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace input_output
{

//! Get SPARTA data files path.
/*!
 * Returns path in which SPARTA data files are located.
 * \return Path containing SPARTA data files.
 */
static inline std::string getSpartaDataPath( )
{
    return getTudatRootPath( ) + "External/SPARTA/";
}

//! Get SPARTA output directory name.
/*!
 * Returns name of directory with SPARTA results.
 * \return Name of directory with SPARTA results.
 */
static inline std::string getSpartaOutputDirectory( )
{
    return "results";
}

//! Get SPARTA output directory path.
/*!
 * Returns path to directory with SPARTA results.
 * \return path to directory with SPARTA results.
 */
static inline std::string getSpartaOutputPath( )
{
    return getSpartaDataPath( ) + getSpartaOutputDirectory( );
}

//! Get path to SPARTA input file.
/*!
 * Returns path to SPARTA input file.
 * \return Path to SPARTA input file.
 */
static inline std::string getSpartaInputFile( )
{
    return getSpartaDataPath( ) + "in.sparta";
}

//! Get path to SPARTA input file template.
/*!
 * Returns path to SPARTA input file template.
 * \return Path to SPARTA input file template.
 */
static inline std::string getSpartaInputFileTemplate( )
{
    return getSpartaDataPath( ) + "SPARTAInputTemplate.txt";
}

//! Get path to SPARTA geometry file for internal use.
/*!
 * Returns path to SPARTA geometry file for internal use.
 * \return Path to SPARTA geometry file for internal use.
 */
static inline std::string getSpartaInternalGeometryFile( )
{
    return getSpartaDataPath( ) + "data.shape";
}

} // namespace input_output

} // namespace tudat

#endif // TUDAT_SPARTA_INPUT_OUTPUT_H
