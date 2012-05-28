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

#ifndef TUDAT_VEHICLE_MODELS_H
#define TUDAT_VEHICLE_MODELS_H

#include <iostream>

#include "Tudat/Astrodynamics/Bodies/body.h"
#include "Tudat/Astrodynamics/Bodies/vehicleExternalModel.h"

namespace tudat
{
namespace bodies
{

//! Vehicle class.
/*!
 * Class that represents the physical model of the vehicle. Subsystem
 * objects should be created externally and then set by the corresponding
 * function, making the VehicleExternalModel class a container for the different
 * properties.
 */
class Vehicle : public Body
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    Vehicle( )
        : externalModel_( boost::shared_ptr< VehicleExternalModel >( ) ),
          isExternalModelSet_( false )
    { }

    //! Set the external model of the vehicle.
    /*!
     * Sets the external model of the vehicle.
     * \param externalModel Vehicle external model.
     */
    void setExternalModel( boost::shared_ptr< VehicleExternalModel > externalModel )
    {
        externalModel_ = externalModel;
        isExternalModelSet_ = true;
    }

    //! Get external model of the vehicle.
    /*!
     * Returns the external model of the vehicle.
     * \return Pointer to vehicle external model.
     */
    boost::shared_ptr< VehicleExternalModel > getExternalModel( ) { return externalModel_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information; prints what
     * property models have been set for the vehicle model.
     * \param stream Stream object.
     * \param vehicle Vehicle.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, Vehicle& vehicle )
    {
        stream << "This is a vehicle; the following properties have been set: " << std::endl;
        if ( vehicle.isExternalModelSet_ )
        {
            stream << "External model" << std::endl;
        }
        return stream;
    }

protected:

private:

    //! Pointer to external model.
    /*!
     * Pointer to object that represents vehicle exterior.
     */
    boost::shared_ptr< VehicleExternalModel > externalModel_;

    //! Flag that indicates if external model has been set.
    /*!
     * Boolean that is true if an external model has been set.
     */
    bool isExternalModelSet_;
};

} // namespace bodies
} // namespace tudat

#endif // TUDAT_VEHICLE_MODELS_H
