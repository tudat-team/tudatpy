/*    Copyright (c) 2010 Delft University of Technology.
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
 
#include <iostream>
#include <vector>
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

//! Lists all files in directory.
std::vector< boost::filesystem3::path > listAllFilesInDirectory(
    const boost::filesystem3::path& directory, bool isRecurseIntoSubdirectories )
{
    // Declare local variables.
    std::vector < boost::filesystem3::path > listOfFileNamesWithPath_;

    if ( boost::filesystem3::exists( directory ) )
    {
        boost::filesystem3::directory_iterator iteratorPastEndOfDirectory_;

        for ( boost::filesystem3::directory_iterator directoryIterator_( directory );
              directoryIterator_ != iteratorPastEndOfDirectory_ ; ++directoryIterator_ )
        {
            if ( boost::filesystem3::is_directory( *directoryIterator_ ) )
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
