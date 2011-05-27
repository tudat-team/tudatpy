/*! \file vehicle.cpp
 *    This file contains the definition of the Vehicle class, which is a
 *    container for the different subsystems and characteristics of the
 *    vehicle. Currently only the VehicleExternalModel can be set as a
 *    property.
 *
 *    Path              : /Astrodynamics/Bodies/Vehicles/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 10 September, 2010
 *    Last modified     : 11 January, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to comments, 80 lines rule, etc.
 *      110112    K. Kumar          Minor comment changes.
 */

// Include statements.
#include "vehicle.h"

// Using declarations.
using std::endl;

//! Default constructor.
Vehicle::Vehicle( )
{
    // External model is not set when constructor is called.
    isExternalModelSet_ = false;
}

//! Default destructor.
Vehicle::~Vehicle( )
{
}

//! Function to set the external model.
void Vehicle::setExternalModel( VehicleExternalModel& externalModel )
{
    // Sets the external mode and the isExternalModelSet_ boolean to true.
    pointerToExternalModel_ = &externalModel;
    isExternalModelSet_ = true;
}

//! Function to retrieve the externalModel.
VehicleExternalModel* Vehicle::getPointerToExternalModel( )
{
    // Only return the external model if one is set.
    if ( isExternalModelSet_ == true )
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
    stream << "This is a vehicle; the following properties have been set: "
           << endl;

    // Check if external model is set and if so, output to stream.
    if ( vehicle.isExternalModelSet_ == true )
    {
        stream << "External model" << endl;
    }

    // Return stream.
    return stream;
}

// End of file.
