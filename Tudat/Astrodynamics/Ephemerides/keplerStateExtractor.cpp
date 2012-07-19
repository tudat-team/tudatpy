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
 *      111103    S. Billemont      First creation of code.
 *      111205    T. Secretin       Added functionalities to previous shell code.
 *
 *    References
 *
 *    Notes
 *
 */

#include <stdexcept>

#include <boost/make_shared.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace ephemerides
{

//! Extract the Keplerian Elements.
boost::shared_ptr< KeplerianElements > KeplerStateExtractor::extract(
        ParsedDataLineMapPtr dataLineMap )
{    
    // Short-hand notation.
    namespace parsed_data_vector_utilities = tudat::input_output::parsed_data_vector_utilities;

    // Create a new KeplerianElements object.
    boost::shared_ptr< KeplerianElements > keplerianElements =
            boost::make_shared< KeplerianElements >( );

    // Find and set semi-major axis.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 tudat::input_output::field_types::state::semiMajorAxis ) )
    {
        keplerianElements->setSemiMajorAxis(
                    parsed_data_vector_utilities::getField< double >(
                        dataLineMap, tudat::input_output::field_types::state::semiMajorAxis ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error( "No semi-major axis entry found." ) ) );
    }

    // Find and set eccentricity.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 tudat::input_output::field_types::state::eccentricity ) )
    {
        keplerianElements->setEccentricity(
                    parsed_data_vector_utilities::getField< double >( dataLineMap,
                                       tudat::input_output::field_types::state::eccentricity ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error( "No eccentricity entry found." ) ) );
    }

    // Find and set inclination.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 tudat::input_output::field_types::state::inclination ) )
    {
        keplerianElements->setInclination(
                    parsed_data_vector_utilities::getField< double >( dataLineMap,
                                        tudat::input_output::field_types::state::inclination ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error( "No inclination entry found." ) ) );
    }

    // Find and set longitude of ascending node.
    if ( checkOptionalFieldType( dataLineMap, 1,
                              tudat::input_output::field_types::state::longitudeOfAscendingNode ) )
    {
        keplerianElements->setLongitudeOfAscendingNode(
                    parsed_data_vector_utilities::getField< double >( dataLineMap,
                            tudat::input_output::field_types::state::longitudeOfAscendingNode ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                            std::runtime_error(
                                        "No longitude of ascending node entry found." ) ) );
    }

    // Find and set argument of periapsis.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 tudat::input_output::field_types::state::argumentOfPeriapsis ) )
    {
        keplerianElements->setArgumentOfPeriapsis(
                    parsed_data_vector_utilities::getField< double >( dataLineMap,
                                tudat::input_output::field_types::state::argumentOfPeriapsis ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                                  std::runtime_error(
                                        "No argument of periapsis entry found." ) ) );
    }

    // Find and set true anomaly.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 tudat::input_output::field_types::state::trueAnomaly ) )
    {
        keplerianElements->setTrueAnomaly(
                    parsed_data_vector_utilities::getField< double >( dataLineMap,
                                        tudat::input_output::field_types::state::trueAnomaly ) );
    }

    // If the true anomaly is not present, check for mean anomaly.
    else if ( checkOptionalFieldType( dataLineMap, 1,
                                      tudat::input_output::field_types::state::meanAnomaly ) )
    {
        // Store mean anomaly.
        double meanAnomaly = parsed_data_vector_utilities::getField< double >( dataLineMap,
                                            tudat::input_output::field_types::state::meanAnomaly );

        // Retrieve eccentricity.
        double eccentricity = keplerianElements->getEccentricity( );

        // Create Newton-Raphson object.
        boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson =
                boost::make_shared< tudat::NewtonRaphson >( );

        // Declare mean to eccentric anomaly conversion class.
        orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
                convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly, newtonRaphson );

        // Convert to eccentric anomaly.
        double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

        // Convert eccentric anomaly to true anomaly and set the latter.
        keplerianElements->setTrueAnomaly(
                    orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                        eccentricAnomaly, eccentricity ) );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                        "No true anomaly or mean anomaly entries found." ) ) );
    }

    return keplerianElements;
}

} // namespace ephemerides
} // namespace tudat
