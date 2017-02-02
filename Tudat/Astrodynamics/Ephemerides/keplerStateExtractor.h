/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      111103    S. Billemont      File created.
 *      111205    T. Secretin       Added functionalities to previous shell code.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_KEPLER_STATE_EXTRACTOR_H
#define TUDAT_KEPLER_STATE_EXTRACTOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/extractor.h"

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Keplerian State extractor class.
/*!
 * This class contains the functionality of extracting the Keplerian elements from a parsed data
 * line map.
 */
class KeplerStateExtractor : public input_output::Extractor< basic_mathematics::Vector6d >
{
public:

    //! Extract the Keplerian Elements.
    /*!
     * Returns a KeplerianElements object containing the orbital parameters found in the input data
     * line map. If one of the elements is not present, the function throws an exception.
     * If no true anomaly field type is found, the function searches for the mean
     * anomaly field type and performs the necessary conversions. The if-statements can be easily
     * extended to deal with other non-standard field types.
     *
     * \param dataLineMap Data map corresponding to a parsed line.
     * \return A KeplerianElements object containing the orbital parameters found in the data map.
     */
    boost::shared_ptr< basic_mathematics::Vector6d > extract( ParsedDataLineMapPtr dataLineMap );

protected:

private:
};

//! Typedef for shared-pointer to KeplerStateExtractor object.
typedef boost::shared_ptr< KeplerStateExtractor > KeplerStateExtractorPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_KEPLER_STATE_EXTRACTOR_H
