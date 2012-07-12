/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      120712    K. Kumar          Updated use of filesystem3 boost-namespace to filesystem.
 *
 *    References
 *
 */
 
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

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
