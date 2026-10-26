/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_CORE_H
#define TUDAT_SERIALIZATION_CORE_H

/**
 * @file serialization/core.h
 * @brief Minimal archive-independent Cereal support used by serialization-aware model headers.
 *
 * Keep archive implementations, streams, provenance handling, pybind11 and Tudat polymorphic
 * registrations out of this header. Public model headers may include this file when their
 * intrusive save/load declarations require Cereal, without making every consumer parse the
 * binary and JSON archive implementations.
 */

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#endif  // TUDAT_SERIALIZATION_CORE_H
