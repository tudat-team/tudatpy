/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      120521    T. Secretin       File created.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_CARTESIAN_STATE_EXTRACTOR_H
#define TUDAT_CARTESIAN_STATE_EXTRACTOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/extractor.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace ephemerides
{

//! Cartesian State extractor class.
/*!
 * This class contains the functionality of extracting the Cartesian elements from a parsed data
 * line map.
 */
class CartesianStateExtractor : public input_output::Extractor< basic_mathematics::Vector6d >
{
public:

    //! Extract the Cartesian elements.
    /*!
     * Returns a CartesianElements object containing the cartesian elements found in the input data
     * line map. If one of the elements is not present, the function throws an exception.
     *
     * \param dataLineMap Data map corresponding to a parsed line.
     * \return A CartesianElements object containing the orbital parameters found in the data map.
     */
    boost::shared_ptr< basic_mathematics::Vector6d > extract( ParsedDataLineMapPtr dataLineMap );

protected:

private:
};

//! Typedef for shared-pointer to CartesianStateExtractor object.
typedef boost::shared_ptr< CartesianStateExtractor > CartesianStateExtractorPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_CARTESIAN_STATE_EXTRACTOR_H
