/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/preProcessFdetsFile.h"

namespace tudat
{

namespace input_output
{

data::PlainLinkDefinition getFdetsLinkEnds( const std::string& spacecraftName,
                                            const std::string& earthName,
                                            const std::string& transmittingStationName,
                                            const std::string& receivingStationName )
{
    return { { { earthName, transmittingStationName }, "transmitter" },
             { { spacecraftName, "" }, "reflector_1" },
             { { earthName, receivingStationName }, "receiver" } };
}

}  // namespace input_output

}  // namespace tudat
