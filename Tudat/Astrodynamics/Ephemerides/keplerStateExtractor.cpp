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
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *      130317    K. Kumar          Updated allocation of Vector6d to fix problem with Eigen types.
 *
 *    References
 *      Eigen. http://eigen.tuxfamily.org/dox/TopicUnalignedArrayAssert.html, last accessed: 18th
 *          March, 2013.
 *
 *    Notes
 *      This implementation uses the KeplerianElements class, which is marked for deprecation. The
 *      code will be updated to use a simple "Vector6d" object from the Eigen library instead.
 *
 */

#include <stdexcept>

#include <boost/exception/all.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace ephemerides
{

//! Extract the Keplerian Elements.
boost::shared_ptr< basic_mathematics::Vector6d > KeplerStateExtractor::extract(
        ParsedDataLineMapPtr dataLineMap )
{
    // Short-hand notation.
    namespace parsed_data_vector_utilities = input_output::parsed_data_vector_utilities;
    using basic_mathematics::Vector6d;

    // Create a new KeplerianElements object.
    boost::shared_ptr< Vector6d > keplerianElements
            = boost::allocate_shared< Vector6d >( Eigen::aligned_allocator< Vector6d >( ) );

    // Find and set semi-major axis.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::semiMajorAxis ) )
    {
        ( *keplerianElements )( orbital_element_conversions::semiMajorAxisIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::semiMajorAxis );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( "No semi-major axis entry found." ) ) );
    }

    // Find and set eccentricity.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::eccentricity ) )
    {
        ( *keplerianElements )( orbital_element_conversions::eccentricityIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::eccentricity );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( "No eccentricity entry found." ) ) );
    }

    // Find and set inclination.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::inclination ) )
    {
        ( *keplerianElements )( orbital_element_conversions::inclinationIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::inclination );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( "No inclination entry found." ) ) );
    }

    // Find and set longitude of ascending node.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::longitudeOfAscendingNode ) )
    {
        ( *keplerianElements )( orbital_element_conversions::longitudeOfAscendingNodeIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::longitudeOfAscendingNode );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error(
                                        "No longitude of ascending node entry found." ) ) );
    }

    // Find and set argument of periapsis.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::argumentOfPeriapsis ) )
    {
        ( *keplerianElements )( orbital_element_conversions::argumentOfPeriapsisIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::argumentOfPeriapsis );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error(
                                        "No argument of periapsis entry found." ) ) );
    }

    // Find and set true anomaly.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::trueAnomaly ) )
    {
        ( *keplerianElements )( orbital_element_conversions::trueAnomalyIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::trueAnomaly );
    }

    // If the true anomaly is not present, check for mean anomaly.
    else if ( checkOptionalFieldType( dataLineMap, 1,
                                      input_output::field_types::state::meanAnomaly ) )
    {
        // Store mean anomaly.
        const double meanAnomaly = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::meanAnomaly );

        // Retrieve eccentricity.
        const double eccentricity
                = ( *keplerianElements )( orbital_element_conversions::eccentricityIndex );

        // Declare mean to eccentric anomaly conversion class.
        orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
                convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly );

        // Convert to eccentric anomaly.
        const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

        // Convert eccentric anomaly to true anomaly and set the latter.
        ( *keplerianElements )( orbital_element_conversions::trueAnomalyIndex )
                = orbital_element_conversions::
                convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error(
                                        "No true anomaly or mean anomaly entries found." ) ) );
    }

    return keplerianElements;
}

} // namespace ephemerides
} // namespace tudat
