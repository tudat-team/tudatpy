/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"

namespace tudat
{
namespace simulation_setup
{

template void setRelativisticTimeConverter< double, double >(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > >& conversionSettings,
        const SystemOfBodies& bodyMap );

template void setRelativisticTimeConverters< double, double >(
        const SystemOfBodies& bodyMap,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > > >& converterSettings );

}  // namespace simulation_setup
}  // namespace tudat
