/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantability or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author              Comment
 *      100914    K. Kumar            File created.
 *      100928    K. Kumar            Completed missing comments, changed
 *                                    writeIntegrationHistoryToFile( ) to
 *                                    writePropagationHistoryToFile( ).
 *      100929    B. Romgens          Spelling mistakes corrected and output to
 *                                    file corrected.
 *      100929    K. Kumar            Added checked code written by D. Dirkx.
 *      110202    K. Kumar            Updated writePropagationHistoryToFile( )
 *                                    to work with State*.
 *      120207    D. Dirkx            Split writingOutputToFile to separate free functions
 *
 *    References
 *
 */

#include <fstream>
#include "Tudat/Astrodynamics/States/writeStateHistoryToFile.h"

namespace tudat
{
namespace output
{

//! Write propagation history to file.
void writeStateHistoryToFile(
    std::map< double, State >& propagationHistory, const std::string& outputFilename )
{
    std::ofstream outputFile_;

    // Declare local variables.
    // Declare iterator for propagation history.
    std::map< double, State >::iterator iteratorPropagationHistory_;

    // Open output file.
    outputFile_.open( outputFilename.c_str( ) );

    // Loop over map of propagation history.
    for ( iteratorPropagationHistory_ = propagationHistory.begin( );
          iteratorPropagationHistory_ != propagationHistory.end( );
          iteratorPropagationHistory_++ )
    {
        // Print map key to output file.
        outputFile_ << iteratorPropagationHistory_->first;

        // Loop over map data.
        for ( int i = 0;
              i < iteratorPropagationHistory_->second.state.rows( ); i++ )
        {
            // Print map data to file.
            outputFile_.precision( 10 );
            outputFile_ << ", " << iteratorPropagationHistory_->second.state( i );
        }

        // End line of output file.
        outputFile_ << std::endl;
    }

    // Close output file.
    outputFile_.close( );
}

} // namespace output
} // namespace tudat
