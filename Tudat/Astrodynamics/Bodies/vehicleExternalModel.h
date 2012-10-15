/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified comments, 80 lines rule, etc.
 *      100928    D. Dirkx          Modifications following first checking iteration.
 *      110112    K. Kumar          Minor comment changes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#ifndef TUDAT_VEHICLE_EXTERNAL_MODELS_H
#define TUDAT_VEHICLE_EXTERNAL_MODELS_H

#include <iostream>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
{
namespace bodies
{

//! Vehicle external model class.
/*!
 * Class that contains the properties of the external portion of the vehicle. Current properties
 * that can be set are: array of pointers to VehiclePartExternalModels, to be able to set different
 * (types of) properties for different vehicle parts.
 */
class VehicleExternalModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    VehicleExternalModel( )
        : vehicleGeometry_( boost::shared_ptr< SurfaceGeometry >( ) ),
          isGeometrySet_( false )
    { }

    //! Set external surface geometry.
    /*!
     * Sets the external vehicle surface geometry. Geometry object is to be created externally.
     * \param vehicleGeometry pointer to shape that is to be set.
     */
    void setVehicleGeometry( boost::shared_ptr< SurfaceGeometry > vehicleGeometry )
    {
        vehicleGeometry_ = vehicleGeometry;
        isGeometrySet_ = true;
    }

    //! Get external surface geometry.
    /*!
     * Returns the external surface geometry.
     * \return Pointer to surface geometric shape.
     */
    boost::shared_ptr< SurfaceGeometry > getVehicleExternalGeometry( ) { return vehicleGeometry_; }

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
        if ( vehicleExternalModel.isGeometrySet_ )
        {
            stream << "Surface geometry" << std::endl;
        }
        return stream;
    }

protected:

private:

    //! Pointer to vehicle surface geometry.
    /*!
     * Pointer to geometric surface shape that represents the vehicle's external shape.
     */
    boost::shared_ptr< SurfaceGeometry > vehicleGeometry_;

    //! Flag that indicates if geometric shape has been set.
    /*!
     * Boolean that is true if a geometric shape has been set.
     */
    bool isGeometrySet_;
};

} // namespace bodies
} // namespace tudat

#endif // TUDAT_VEHICLE_EXTERNAL_MODELS_H
