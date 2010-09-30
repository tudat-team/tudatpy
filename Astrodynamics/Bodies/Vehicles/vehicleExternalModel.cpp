/*! \file vehicleExternalModel.cpp
 * This file contains the definition of the VehicleExternalModel class.
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
 *
 *
 *  Notes
 *
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
 *    100928   D. Dirkx                    Modifications following first
 *                                         checking iteration.
 *
 */

#include "vehicleExternalModel.h"

//! Default constructor.
VehicleExternalModel::VehicleExternalModel( )
{
    // No geometry is set when constructor is called.
    isGeometrySet_ = 0;
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
    isGeometrySet_ = 1;
}

//! Function to retrieve external geometry.
GeometricShape* VehicleExternalModel::getPointerToVehicleGeometry( )
{
    // Only return the geometry when one is set.
    if( isGeometrySet_ == 1 )
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
    stream << "This is an external model, the following properties"<<
              " have been set:" << std::endl;
    if( vehicleExternalModel.isGeometrySet_ == true )
    {
        stream << "Surface geometry" << std::endl;
    }
    return stream;
}

// End of file.
