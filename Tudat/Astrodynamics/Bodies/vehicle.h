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
