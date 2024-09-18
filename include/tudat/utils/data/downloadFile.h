/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

namespace tudat {
namespace utils {
namespace data {

std::string download_file(const char *remote_url,
                          const char *cache = "true",
                          int verbosity = 1,
                          bool try_unzip = true,
                          const char *prefix = nullptr);

} // namespace data
} // namespace utils
} // namespace tudat
