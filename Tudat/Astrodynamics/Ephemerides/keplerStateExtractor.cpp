/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerStateExtractor.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace ephemerides
{

//! Extract the Keplerian Elements.
std::shared_ptr< Eigen::Vector6d > KeplerStateExtractor::extract(
        ParsedDataLineMapPtr dataLineMap )
{
    // Short-hand notation.
    namespace parsed_data_vector_utilities = input_output::parsed_data_vector_utilities;
    using Eigen::Vector6d;

    // Create a new KeplerianElements object.
    std::shared_ptr< Vector6d > keplerianElements
            = std::allocate_shared< Vector6d >( Eigen::aligned_allocator< Vector6d >( ) );

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
        throw std::runtime_error( "No semi-major axis entry found." );
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
        throw std::runtime_error( "No eccentricity entry found." );
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
        throw std::runtime_error( "No inclination entry found." );
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
        throw std::runtime_error( "No longitude of ascending node entry found." );
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
        throw std::runtime_error( "No argument of periapsis entry found." );
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

        // Convert to eccentric anomaly.
        const double eccentricAnomaly = orbital_element_conversions::
                convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly );

        // Convert eccentric anomaly to true anomaly and set the latter.
        ( *keplerianElements )( orbital_element_conversions::trueAnomalyIndex )
                = orbital_element_conversions::
                convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );
    }

    else
    {
        throw std::runtime_error( "No true anomaly or mean anomaly entries found." );
    }

    return keplerianElements;
}

} // namespace ephemerides
} // namespace tudat
