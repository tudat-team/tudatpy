/*! \file vehicle.cpp
 *  This file contains the definition of the VehicleModel class,
 *  which is a container for the different subsystems and
 *  characteristics of the vehicle. Currently only the
 *  VehicleExternalModel can be set as a property
 *
 *  Path              : Astrodynamics/Bodies/Vehicles/
 *  Version           : 3
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : J. Melman
 *  Affiliation       : Delft University of Technology
 *  E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *  Date created      : 10 September, 2010
 *  Last modified     : 28 September, 2010
 *
 *  References
 *
 *  Notes
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    100910   D. Dirkx                    First version of file
 *    100915   D. Dirkx                    Modified to correct comments, 80
 *                                         lines rule, etc.
 *
 */

#include "vehicle.h"

//! Default constructor.
Vehicle::Vehicle( )
{
    // External model is not set when constructor is called.
    isExternalModelSet_ = 0;
}

//! Default destructor.
Vehicle::~Vehicle( )
{
}

//! Function to set the external model.
void Vehicle::setExternalModel(VehicleExternalModel& externalModel)
{
    // Sets the external mode and the isExternalModelSet_ boolean to true.
    pointerToExternalModel_ = &externalModel;
    isExternalModelSet_ = 1;
}

//! Function to retrieve the externalModel.
VehicleExternalModel* Vehicle::getPointerToExternalModel( )
{
    // Only return the external model if one is set.
    if (isExternalModelSet_ == 1)
    {
        return pointerToExternalModel_;
    }
    else
    {
        return NULL;
    }
}

std::ostream& operator<<( std::ostream& stream, Vehicle& vehicle )
{
    stream << "This is a vehicle; the following properties have been set:" << std::endl;
    if( vehicle.isExternalModelSet_ == true)
    {
        stream << "External model" << std::endl;
    }
    return stream;
}

// End of file.
