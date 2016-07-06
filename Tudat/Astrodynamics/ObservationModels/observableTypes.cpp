/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{

std::string getObservableName( const ObservableType observableType )
{
    std::string observableName;
    switch( observableType )
    {
    case oneWayRange:
        observableName = "OneWayRange";
        break;
    case angular_position:
        observableName = "AngularPosition";
        break;
    case position_observable:
        observableName = "CartesianPosition";
        break;
    default:
        std::cerr<<"Error, could not find observable type "<<observableType<<" when getting name of type"<<std::endl;
    }
    return observableName;
}

ObservableType getObservableType( const std::string& observableName )
{
    ObservableType observableType;

    if( observableName == "OneWayRange" )
    {
        observableType = oneWayRange;
    }
    else if( observableName == "AngularPosition" )
    {
        observableType = angular_position;
    }
    else if( observableName == "CartesianPosition" )
    {
        observableType = position_observable;
    }
    else
    {
        std::cerr<<"Error, could not find observable name "<<observableName<<" when getting type of name"<<std::endl;
    }

    return observableType;
}

}

}
