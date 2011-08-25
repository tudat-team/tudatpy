/*! \file vehicleExternalModel.cpp
 *    This file contains the definition of the VehicleExternalModel class.
 *
 *    Path              : /Astrodynamics/Bodies/Vehicles/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
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
 *      100915    D. Dirkx          Modified comments, 80 lines rule, etc.
 *      100928    D. Dirkx          Modifications following first checking
 *                                  iteration.
 *      110112    K. Kumar          Minor comment changes.
 */

// Include statements.
#include "vehicleExternalModel.h"

// Using declarations.
using std::endl;

//! Default constructor.
VehicleExternalModel::VehicleExternalModel( )
{
    // No geometry is set when constructor is called.
    isGeometrySet_ = false;
}

//! Default destructor.
VehicleExternalModel::~VehicleExternalModel( )
{
}

//! Function to set external geometry.
void VehicleExternalModel::setVehicleGeometry( GeometricShape& vehicleGeometry )
{
    // Sets the geometry and the isGeometrySet_ boolean to true.
    pointerToVehicleGeometry_ = &vehicleGeometry;
    isGeometrySet_ = true;
}

//! Function to retrieve external geometry.
GeometricShape* VehicleExternalModel::getVehicleExternalGeometry( )
{
    // Only return the geometry when one is set.
    if ( isGeometrySet_ == true )
    {
        return pointerToVehicleGeometry_;
    }

    else
    {
        return NULL;
    }
}

//! Overloaded ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          VehicleExternalModel& vehicleExternalModel )
{
    stream << "This is an external model, the following properties"
           << " have been set:" << endl;

    // Check if external geometry is set and if so, output to stream.
    if ( vehicleExternalModel.isGeometrySet_ == true )
    {
        stream << "Surface geometry" << endl;
    }

    // Return stream.
    return stream;
}

// End of file.
