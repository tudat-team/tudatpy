/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TEST_FILE_UTILITIES_H
#define TUDAT_TEST_FILE_UTILITIES_H

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>

namespace tudat
{
namespace unit_tests
{

inline std::string createTemporaryFilePath( const std::string& prefix, const std::string& suffix )
{
    static std::atomic< std::uint64_t > counter{ 0 };
    const auto timestamp = std::chrono::steady_clock::now( ).time_since_epoch( ).count( );
    const auto sequence = counter.fetch_add( 1, std::memory_order_relaxed );
    return ( std::filesystem::temp_directory_path( ) /
             ( prefix + "-" + std::to_string( timestamp ) + "-" + std::to_string( sequence ) + suffix ) )
            .string( );
}

}  // namespace unit_tests
}  // namespace tudat

#endif  // TUDAT_TEST_FILE_UTILITIES_H
