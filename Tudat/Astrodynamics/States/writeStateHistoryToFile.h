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

#include <map>
#include <string>

#include "Tudat/Astrodynamics/States/state.h"

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
void writeStateHistoryToFile( std::map< double, astrodynamics::states::State >& propagationHistory,
                              const std::string& outputFilename );

} // namespace output
} // namespace tudat

#endif // TUDAT_WRITE_STATE_HISTORY_TO_FILE_H
