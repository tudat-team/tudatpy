/*! \file basicFunctions.cpp
 *    Source file that defines the basicFunctions namespace, containing all
 *    basic functions contained in Tudat.
 *
 *    Path              : /Basics/
 *    Version           : 11
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : S. Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Date created      : 1 September, 2010
 *    Last modified     : 17 November, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *      Need to add bounds checking with respect to targetValue for
 *      computeNearestLeftNeighborUsingBinarySearch( ).
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      100916    D. Dirkx          Added minor comments during checking.
 *      100928    K. Kumar          Small comment modifications.
 *      100929    K. Kumar          Small comment modifications.
 *      110202    K. Kumar          Added overload for map with State* for
 *                                  computeNearestLeftNeighborUsingBinarySearch().
 *      110803    J. Leloux         Added convertStringToTemplate.
 *      110805    J. Leloux         Added outputCurrentRunningTime().
 *      110807    K. Kumar          Minor comment modifications.
 *      110810    J. Leloux         Minor comment modifications.
 *      110913    K. Kumar          Implemented automatic root-path functions based on
 *                                  suggestions by M. Persson.
 *      111117    K. Kumar          Added listAllFilesInDirectory() function.
 */

// Include statements.
#include <iterator>
#include "Basics/basicFunctions.h"

//! Tudat library namespace.
namespace tudat
{

//! Basic functions namespace.
namespace basic_functions
{

//! Get root-path for Tudat library.
string getRootPath( )
{
#ifdef TUDAT_CUSTOM_ROOT_PATH

    return string( TUDAT_CUSTOM_ROOT_PATH );
#else
    // Declare file path string assigned to filePath.
    string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( )
                             - string( "/Basics/basicFunctions.cpp" ).length( ) ) + "/";
#endif
}

//! Nearest left neighbor binary search.
int computeNearestLeftNeighborUsingBinarySearch(
        VectorXd& vectorOfSortedData,
        double& targetValueInVectorOfSortedData )
{
    // Declare local variables.
    // Declare bounds of vector of sorted data and current position.
    int leftLimitOfVectorOfSortedData = 0;
    int rightLimitOfVectorOfSortedData = vectorOfSortedData.rows( ) - 1;
    int currentPositionInVectorOfSortedData;

    // Check if data is sorted in ascending order.
    // ( true if ascending, else false ).
    bool isVectorOfSortedDataAscending
            = ( vectorOfSortedData[ rightLimitOfVectorOfSortedData ]
                >= vectorOfSortedData[ leftLimitOfVectorOfSortedData ] );

    // Loop through vector of sorted data until left and right limits
    // are neighbours.
    while ( rightLimitOfVectorOfSortedData
            - leftLimitOfVectorOfSortedData > 1 )
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInVectorOfSortedData
                = ( rightLimitOfVectorOfSortedData
                    + leftLimitOfVectorOfSortedData ) >> 1;

        // Check which limit to replace ( if ascending and target datum
        // is in right half, replace left limit ).
        if ( targetValueInVectorOfSortedData
             >= vectorOfSortedData[ currentPositionInVectorOfSortedData ]
             && isVectorOfSortedDataAscending )
        {
            // Set left limit to current position in vector of sorted data.
            leftLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }

        else
        {
            // Set right limit to current position in vector of sorted data.
            rightLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }
    }

    // Set current position to left limit.
    currentPositionInVectorOfSortedData = leftLimitOfVectorOfSortedData;

    // Return current position in vector.
    return currentPositionInVectorOfSortedData;
}

//! Nearest left neighbor binary search.
int computeNearestLeftNeighborUsingBinarySearch(
        std::map < double, VectorXd >& sortedIndepedentAndDependentVariables,
        double& targetValueInMapOfData )
{
    // Declare local variables.
    // Declare bounds of key of map of data and current position.
    int leftLimitOfKeyOfMapOfData = 0;
    int rightLimitOfKeyOfMapOfData = sortedIndepedentAndDependentVariables
                                     .size( ) - 1;
    int currentPositionInKeyOfMapOfData;

    // Declare map iterator
    std::map < double, VectorXd >::iterator mapIterator;

    // Loop through vector of sorted data until left and right limits
    // are neighbours
    while ( rightLimitOfKeyOfMapOfData - leftLimitOfKeyOfMapOfData > 1 )
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInKeyOfMapOfData
                = ( rightLimitOfKeyOfMapOfData
                   + leftLimitOfKeyOfMapOfData ) >> 1;

        // Set map iterator to begin begin of map of sorted independent and
        // dependent variables.
        mapIterator = sortedIndepedentAndDependentVariables.begin( );

        // Advance iterator to location of current position in key of map of
        // data.
        advance( mapIterator, currentPositionInKeyOfMapOfData );

        // Check that target value lies to the right of lower bound.
        if ( targetValueInMapOfData
             >= mapIterator->first )
        {
            // Set left limit to current position in map of data.
            leftLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }

        else
        {
            // Set right limit to current position in map of data.
            rightLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }
    }

    // Set current position to left limit.
    currentPositionInKeyOfMapOfData = leftLimitOfKeyOfMapOfData;

    // Return current position in map of data.
    return currentPositionInKeyOfMapOfData;
}

//! Nearest left neighbor binary search.
int computeNearestLeftNeighborUsingBinarySearch(
        std::map < double, State* >& sortedIndepedentAndDependentVariables,
        double& targetValueInMapOfData )
{
    // Declare local variables.
    // Declare bounds of key of map of data and current position.
    int leftLimitOfKeyOfMapOfData = 0;
    int rightLimitOfKeyOfMapOfData = sortedIndepedentAndDependentVariables
                                     .size( ) - 1;
    int currentPositionInKeyOfMapOfData;

    // Declare map iterator
    std::map < double, State* >::iterator mapIterator;

    // Loop through State objects of sorted data until left and right limits
    // are neighbours.
    while ( rightLimitOfKeyOfMapOfData - leftLimitOfKeyOfMapOfData > 1 )
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInKeyOfMapOfData
                = ( rightLimitOfKeyOfMapOfData
                   + leftLimitOfKeyOfMapOfData ) >> 1;

        // Set map iterator to begin begin of map of sorted independent and
        // dependent variables.
        mapIterator = sortedIndepedentAndDependentVariables.begin( );

        // Advance iterator to location of current position in key of map of
        // data.
        advance( mapIterator, currentPositionInKeyOfMapOfData );

        // Check that target value lies to the right of lower bound.
        if ( targetValueInMapOfData
             >= mapIterator->first )
        {
            // Set left limit to current position in map of data.
            leftLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }

        else
        {
            // Set right limit to current position in map of data.
            rightLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }
    }

    // Set current position to left limit.
    currentPositionInKeyOfMapOfData = leftLimitOfKeyOfMapOfData;

    // Return current position in map of data.
    return currentPositionInKeyOfMapOfData;
}

//! Write the current running time and status to vector.
vector< string > outputCurrentRunningTime( clock_t start_clock, const string& status )
{
    // Declare local variables.
    // Declare vector of strings to store current running time and status of executed application.
    vector< string > runningTimeAndStatusContainer_;

    // Declare stringstream.
    std::stringstream currentRunningTimeStatement_;

    // Set current and running clock times.
    clock_t current_clock_ = clock( );
    clock_t running_clocks_ = current_clock_ - start_clock;

    // Compute current running time.
    double runningTime_ = running_clocks_;

    // Create output string.
    currentRunningTimeStatement_ << "Current running time in seconds is: "
                                 << runningTime_ / 1000;

    // Store data in vector of strings.
    runningTimeAndStatusContainer_.push_back( status );
    runningTimeAndStatusContainer_.push_back( currentRunningTimeStatement_.str( ) );

    // Return string container.
    return runningTimeAndStatusContainer_;
}

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

}

}

// End of file.
