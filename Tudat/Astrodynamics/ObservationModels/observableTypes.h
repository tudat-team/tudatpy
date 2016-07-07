/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVABLETYPES_H
#define TUDAT_OBSERVABLETYPES_H

#include <string>

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{


//! Enum for types of observations
enum ObservableType
{
    oneWayRange = 0,
    angular_position = 1,
    position_observable = 2

};

//! Function to get the name (string) associated with a given observable type.
/*!
 * Function to get the name (string) associated with a given observable type.
 * \param observableType Type of observable.
 * \return Name of observable
 */
std::string getObservableName( const ObservableType observableType );

//! Function to get the observable type.ssociated with the name (string) of observable.
/*!
 * Function to get the observable type.ssociated with the name (string) of observable.
 * \param observableName of observable
 * \return observableType Type of observable.
 */
ObservableType getObservableType( const std::string& observableName );

}

}

#endif // TUDAT_OBSERVABLETYPES_H
