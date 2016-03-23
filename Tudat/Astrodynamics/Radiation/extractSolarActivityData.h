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
 *      120607    A. Ronse          Creation of code.
 *
 *    References
 *
 */

#ifndef EXTRACTSOLARACTIVITY_H
#define EXTRACTSOLARACTIVITY_H

#include "Tudat/Astrodynamics/Radiation/solarActivityData.h"
#include "Tudat/InputOutput/extractor.h"

#include <boost/shared_ptr.hpp>

namespace tudat
{
namespace radiation
{
namespace solar_activity
{

//! Solar activity extractor class.
/*!
 * This class extracts the numeric information from a ParsedDataLineMapPtr containing parsed
 * solar activity data designed and placec them in a container of class SolarActivityData
 */
class ExtractSolarActivityData : public tudat::input_output::Extractor<
            tudat::radiation::solar_activity::SolarActivityData >
{

public:
    //! Extracts the solar activity data to a SolarActivityData container.
    /*!
     * Extracts the solar activity data from a "ParsedDataLineMap" object and saves it in a
     * "SolarActivityData" contatiner.
     */
    boost::shared_ptr<tudat::radiation::solar_activity::SolarActivityData> extract(
                tudat::input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr data );

protected:

private:

};

}   // namespace solar_activity
}   // namespace radiation
}   // namespace tudat

#endif  // EXTRACTSOLARACTIVITY_H
