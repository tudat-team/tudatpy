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
 *      130318    K. Kumar          Updated allocation of Vector6d to fix problem with Eigen types.
 *
 *    References
 *      Eigen. http://eigen.tuxfamily.org/dox/TopicUnalignedArrayAssert.html, last accessed: 18th
 *          March, 2013.
 *
 *    Notes
 *
 */

#include <stdexcept>

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Ephemerides/cartesianStateExtractor.h"

namespace tudat
{
namespace ephemerides
{

boost::shared_ptr< basic_mathematics::Vector6d > CartesianStateExtractor::extract(
        ParsedDataLineMapPtr dataLineMap )
{
    // Short-hand notation.
    namespace parsed_data_vector_utilities = input_output::parsed_data_vector_utilities;
    using basic_mathematics::Vector6d;

    // Create a new CartesianElements object.
    boost::shared_ptr< Vector6d > cartesianElements
            = boost::allocate_shared< Vector6d >( Eigen::aligned_allocator< Vector6d >( ) );

    // Find and set Cartesian x coordinate.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianXCoordinate ) )
    {
        ( *cartesianElements )( basic_astrodynamics::xCartesianPositionIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianXCoordinate );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                 std::runtime_error(
                                        "No Cartesian x coordinate entry found." ) ) );
    }

    // Find and set Cartesian y coordinate.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianYCoordinate ) )
    {
        ( *cartesianElements )( basic_astrodynamics::yCartesianPositionIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianYCoordinate );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                 std::runtime_error( "No Cartesian y coordinate entry found." ) ) );
    }

    // Find and set Cartesian z coordinate.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianZCoordinate ) )
    {
        ( *cartesianElements )( basic_astrodynamics::zCartesianPositionIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianZCoordinate );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                 std::runtime_error(
                                        "No Cartesian z coordinate entry found." ) ) );
    }

    // Find and set Cartesian x velocity.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianXVelocity ) )
    {
        ( *cartesianElements )( basic_astrodynamics::xCartesianVelocityIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianXVelocity );
    }
    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error(
                                        "No Cartesian x velocity entry found." ) ) );
    }

    // Find and set Cartesian y velocity.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianYVelocity ) )
    {
        ( *cartesianElements )( basic_astrodynamics::yCartesianVelocityIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianYVelocity );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error(
                                        "No Cartesian y velocity entry found." ) ) );
    }

    // Find and set Cartesian z velocity.
    if ( checkOptionalFieldType( dataLineMap, 1,
                                 input_output::field_types::state::cartesianZVelocity ) )
    {
        ( *cartesianElements )( basic_astrodynamics::zCartesianVelocityIndex )
                = parsed_data_vector_utilities::getField< double >(
                    dataLineMap, input_output::field_types::state::cartesianZVelocity );
    }

    else
    {
        boost::throw_exception( boost::enable_error_info(
                                   std::runtime_error(
                                        "No Cartesian z velocity entry found." ) ) );
    }

    return cartesianElements;
}

} // namespace ephemerides
} // namespace tudat
