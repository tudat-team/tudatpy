/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_COPYSYSTEMOFBODIES_H
#define TUDAT_COPYSYSTEMOFBODIES_H

#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to copy the environment models of a system of bodies, without processing the frame definitions of the copy
/*!
 * Function to copy the environment models of a system of bodies, without processing the frame definitions of the copy. See
 * createSystemOfBodiesCopy for details.
 * \param bodies System of bodies that is to be copied
 * \param printWarnings Boolean denoting whether warnings for environment models that are shared between original and copy
 * are to be printed
 * \return Copy of system of bodies, without frame definitions
 */
SystemOfBodies copySystemOfBodiesWithoutFrameProcessing( const SystemOfBodies& bodies, const bool printWarnings );

//! Function to create a copy of a system of bodies that can be used in a separate thread from the original
/*!
 * Function to create a copy of a system of bodies that can be used in a separate thread from the original. Environment models
 * that store data during the propagation (current states, caches, interpolator look-up indices) are copied, and the links
 * between models are recreated for the copy. Environment models without such data are shared between the original and the copy.
 * Model types for which no copy is implemented, and models containing user-defined functions, are shared between original and
 * copy, for which a warning is printed (if printWarnings is true). Such shared models must not be used in multi-threaded
 * propagation, unless they are thread-safe.
 * \param bodies System of bodies that is to be copied
 * \param printWarnings Boolean denoting whether warnings for environment models that are shared between original and copy
 * are to be printed
 * \return Copy of system of bodies
 */
template< typename StateScalarType = double, typename TimeType = double >
SystemOfBodies createSystemOfBodiesCopy( const SystemOfBodies& bodies, const bool printWarnings = true )
{
    SystemOfBodies bodiesCopy = copySystemOfBodiesWithoutFrameProcessing( bodies, printWarnings );
    bodiesCopy.processBodyFrameDefinitions< StateScalarType, TimeType >( );
    return bodiesCopy;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_COPYSYSTEMOFBODIES_H
