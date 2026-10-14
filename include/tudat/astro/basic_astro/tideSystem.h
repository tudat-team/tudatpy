/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_TIDE_SYSTEM_H
#define TUDAT_TIDE_SYSTEM_H
namespace tudat
{
namespace basic_astrodynamics
{
//! Permanent-tide convention explicitly declared by a gravity-field provider.
enum class TideSystem { unknown, tideFree, zeroTide, meanTide };
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
