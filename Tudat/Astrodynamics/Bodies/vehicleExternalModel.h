/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      100928    D. Dirkx          Modifications following first checking iteration.
 *      110112    K. Kumar          Minor comment changes.
 *
 *    References
 *
 */

#ifndef TUDAT_VEHICLE_EXTERNAL_MODELS_H
#define TUDAT_VEHICLE_EXTERNAL_MODELS_H

#include <iostream>
#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
{

//! Vehicle external model class.
/*!
 * Class that contains the properties of the external portion of the
 * vehicle. Current properties that can be set are:
 * Array of pointers to VehiclePartExternalModels, to be able to set
 * different (types of) properties for different vehicle parts.
 */
class VehicleExternalModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    VehicleExternalModel( ) : pointerToVehicleGeometry_( NULL ), isGeometrySet_( false ) { }

    //! Function to set external surface geometry.
    /*!
     * Function to set the external vehicle surface geometry. Geometry object is to be
     * created externally.
     * \param vehicleGeometry pointer to shape that is to be set.
     */
    void setVehicleGeometry( SurfaceGeometry& vehicleGeometry )
    { pointerToVehicleGeometry_ = &vehicleGeometry; isGeometrySet_ = true; }

    //! Get external surface geometry.
    /*!
     * Returns the external surface geometry.
     * \return Pointer to surface geometric shape.
     */
    SurfaceGeometry* getVehicleExternalGeometry( ) { return pointerToVehicleGeometry_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information; prints whether a
     * surface geometry has been set.
     * \param stream Stream object.
     * \param vehicleExternalModel Vehicle external model.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     VehicleExternalModel& vehicleExternalModel )
    {
        stream << "This is an external model, the following properties are set:" << std::endl;
        if ( vehicleExternalModel.isGeometrySet_ == true )
        { stream << "Surface geometry" << std::endl; }
        return stream;
    }

private:

    //! Pointer to vehicle surface geometry.
    /*!
     * Pointer to geometric surface shape that represents the vehicle's external shape.
     */
    SurfaceGeometry* pointerToVehicleGeometry_;

    //! Flag that indicates if geometric shape has been set.
    /*!
     * Boolean that is true if a geometric shape has been set.
     */
    bool isGeometrySet_;
};

} // namespace tudat

#endif // TUDAT_VEHICLE_EXTERNAL_MODELS_H
