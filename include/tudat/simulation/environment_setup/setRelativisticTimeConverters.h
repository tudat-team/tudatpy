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

#ifndef TUDAT_SET_RELATIVISTIC_TIME_CONVERTERS_H
#define TUDAT_SET_RELATIVISTIC_TIME_CONVERTERS_H

#include <map>
#include <memory>
#include <string>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/relativisticTimeConverterSettings.h"

namespace tudat
{
namespace simulation_setup
{

//! Setup a single relativistic time converter (barycentric -> bodycentric -> topocentric).
/*!
 *  \param conversionSettings Settings defining the converter and all conversion legs.
 *  \param bodies System of bodies in which the converter is created and stored.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverter(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > >& conversionSettings,
        const SystemOfBodies& bodies );

//! Setup multiple relativistic time converters (one per associated body).
/*!
 *  \param bodies System of bodies in which converters are created and stored.
 *  \param converterSettings Map from body name to converter settings.
 */
template< typename StateScalarType = double, typename TimeType = double >
void setRelativisticTimeConverters(
        const SystemOfBodies& bodies,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< StateScalarType, TimeType > > >&
                converterSettings );

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template void setRelativisticTimeConverter< double, double >(
        const std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > >& conversionSettings,
        const SystemOfBodies& bodies );

extern template void setRelativisticTimeConverters< double, double >(
        const SystemOfBodies& bodies,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< double, double > > >& converterSettings );
#endif

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_SET_RELATIVISTIC_TIME_CONVERTERS_H
