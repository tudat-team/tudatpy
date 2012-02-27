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

#ifndef TUDAT_WRITE_STATE_HISTORY_TO_FILE_H
#define TUDAT_WRITE_STATE_HISTORY_TO_FILE_H

#include "Tudat/Astrodynamics/States/state.h"
#include <map>
#include <string>

namespace tudat
{
namespace output
{

//! Write state history to file.
/*!
 * Writes state history stored in a map to file. Map can be for instance of States and times.
 * \param propagationHistory Map of propagation history.
 * \param outputFilename Output filename.
 */
void writeStateHistoryToFile( std::map< double, State >& propagationHistory,
                                    const std::string& outputFilename );

} // namespace output
} // namespace tudat

#endif // TUDAT_WRITE_STATE_HISTORY_TO_FILE_H
