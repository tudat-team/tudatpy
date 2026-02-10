/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_H
#define TUDAT_SERIALIZATION_H

/**
 * @file serialization.h
 * @brief Master header for Tudat Boost serialization infrastructure.
 *
 * Provides core serialization helpers (Eigen, archives, serialize/deserialize
 * utilities). Does NOT include BOOST_CLASS_EXPORT_IMPLEMENT headers, since
 * those must appear in exactly one translation unit.
 *
 * Sub-headers (include directly as needed):
 *   - serialization/base.h          Core infrastructure (Eigen, helpers, archives)
 *   - serialization/observations.h  BOOST_CLASS_EXPORT_IMPLEMENT for observation types
 *                                    (include in exactly ONE .cpp file)
 *   - serialization/propagation.h   BOOST_CLASS_EXPORT_IMPLEMENT for propagation/estimation types
 *                                    (include in exactly ONE .cpp file)
 */

#include "tudat/io/serialization/base.h"

#endif  // TUDAT_SERIALIZATION_H
