/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/mapTextFileReader.h"

namespace tudat
{

namespace input_output
{

template std::map< double, std::vector< double > > readStlVectorMapFromFile< double, double >( const std::string& relativePath,
                                                                                               const std::string& separators,
                                                                                               const std::string& skipLinesCharacter );

template std::map< double, double > readFloatingPointMapFromFile< double, double >( const std::string& relativePath,
                                                                                    const std::string& separators,
                                                                                    const std::string& skipLinesCharacter );

}  // namespace input_output

}  // namespace tudat
