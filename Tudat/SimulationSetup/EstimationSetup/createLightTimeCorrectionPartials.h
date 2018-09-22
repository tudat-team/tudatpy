/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_CREATELIGHTTIMECORRECTIONPARTIALS_H
#define TUDAT_CREATELIGHTTIMECORRECTIONPARTIALS_H

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to create a partial objects from list of light time corrections.
/*!
 * Function to create a partial objects from list of light time corrections.
 * \param lightTimeCorrectionList List of light time corrections.
 * \return List of light-time correction partial objects
 */
std::vector< std::shared_ptr< LightTimeCorrectionPartial > > createLightTimeCorrectionPartials(
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrectionList );

}

}

#endif // TUDAT_CREATELIGHTTIMECORRECTIONPARTIALS_H
