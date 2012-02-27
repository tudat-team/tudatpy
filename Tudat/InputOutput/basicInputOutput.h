/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100902    K. Kumar          File header and footer added.
 *      100916    D. Dirkx          Added minor comments and placeholder tag during checking.
 *      100928    K. Kumar          Added reference and adjusted include statements.
 *      100929    K. Kumar          Changed EigenRoutines.h include statement
 *                                  to linearAlgebra.h and removed placeholder comment.
 *                                  Added small comment modifications.
 *      110202    K. Kumar          Added overload for map with State* for
 *                                  computeNearestLeftNeighborUsingBinarySearch( ).
 *      110803    J. Leloux         Added convertStringToTemplate.
 *      110805    J. Leloux         Added outputCurrentRunningTime.
 *      110810    J. Leloux         Minor comment modifications.
 *      110913    K. Kumar          Implemented automatic root-path functions based on
 *                                  suggestions by M. Persson.
 *      111117    K. Kumar          Added listAllFilesInDirectory( ) function.
 *      120127    K. Kumar          Adapted for Tudat Core.
 *
 *    References
 *
 */

#ifndef TUDAT_BASIC_INPUT_OUTPUT_H
#define TUDAT_BASIC_INPUT_OUTPUT_H

#include <string>
#include <boost/filesystem.hpp>

namespace tudat
{
namespace input_output
{

//! Get root-path for Tudat library.
/*!
 * Returns root-path corresponding with root-directory of Tudat library as a string with
 * trailing slash included.
 * \return Tudat root-path.
 */
static inline std::string getTudatRootPath( )
{
#ifdef TUDAT_CUSTOM_ROOT_PATH
    return std::string( TUDAT_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path in the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "InputOutput/basicInputOutput.h" ).length( ) );
#endif
}

//! List all files in directory.
/*!
 * Lists all files in a given directory. There is a recursion option to allow
 * all files in subdirectories to be listed as well.
 * \param directory Absolute directory path.
 * \param isRecurseIntoSubdirectories Flag to set if algorithm should recurse through
 *          subdirectories. Set to false by default.
 * \return Container of filenames in directory, stored as Boost path variables.
 */
std::vector< boost::filesystem3::path > listAllFilesInDirectory(
    const boost::filesystem3::path& directory, bool isRecurseIntoSubdirectories = false );

} // namespace input_output
} // namespace tudat

#endif // TUDAT_BASIC_INPUT_OUTPUT_H
