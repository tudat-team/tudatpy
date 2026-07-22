/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H

/**
 * @file registrations_propagation.h
 * @brief Retains and initializes propagation polymorphic registrations.
 */

#if TUDAT_BUILD_WITH_SERIALIZATION
#include <cereal/types/polymorphic.hpp>

#include <typeinfo>
#include <utility>
#include <vector>

namespace tudat
{
namespace serialization
{
namespace diagnostics
{

using RegistrationTypeInfo = std::pair< const char*, const std::type_info* >;

//! Temporary CI diagnostic: retrieve RTTI objects as seen by the registration translation unit.
std::vector< RegistrationTypeInfo > getPropagationRegistrationTypeInfo( );

}  // namespace diagnostics
}  // namespace serialization
}  // namespace tudat

CEREAL_FORCE_DYNAMIC_INIT( tudat_serialization_propagation )
#endif

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
