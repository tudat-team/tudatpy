/*! \file vehicleModels.cpp
 * This file contains the definition of the VehicleModel class, which is a container
 * for the different subsystems and characteristics of the vehicle. Currently only
 * the VehicleExternalModel can be set as a property.
 *
 * Path              : Astrodynamics/Bodies/Vehicles/
 * Version           : 3
 * Check status      : Checked
 *
 * Author            : Dominic Dirkx
 * Affiliation       : TU Delft
 * E-mail address    : D.Dirkx@student.tudelft.nl
 *
 * Checker           : J. Melman
 * Affiliation       : Delft University of Technology
 * E-mail address    : J.C.P.Melman@tudelft.nl
 *
 * Date created      : 10 September, 2010
 * Last modified     : 28 September, 2010
 *
 * References
 *
 * Notes
 *   This class should be derived from the class Body.
 *
 * Copyright (c) 2010 Delft University of Technology.
 *
 * This software is protected by national and international copyright.
 * Any unauthorized use, reproduction or modification is unlawful and
 * will be prosecuted. Commercial and non-private application of the
 * software in any form is strictly prohibited unless otherwise granted
 * by the authors.
 *
 * The code is provided without any warranty; without even the implied
 * warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    100910   D. Dirkx                    First version of file
 *    100915   D. Dirkx                    Modified to correct comments, 80
 *                                         lines rule, etc.
 *    100928   D. Dirkx                    Modifications following first
 *                                         checking iteration.
 *
 */

#ifndef VEHICLEMODELS_H
#define VEHICLEMODELS_H

#include <iostream>
#include "vehicleExternalModels.h"

/*!
 * Class that represents the physical model of the vehicle. Subsystem
 * objects should be created externally and then set by the corresponding
 * function, making the VehicleModel class a container for the different
 * properties.
 */
class VehicleModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    VehicleModel( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~VehicleModel( );

    //! Function to set the externalModel of the vehicle.
    /*!
     * Function to set the external model of the vehicle.
     */
    void setExternalModel( VehicleExternalModel& externalModel );

    //! Function to retrieve the externalModel of the vehicle.
    /*!
     * Function to retrieve the external model of the vehicle.
     */
    VehicleExternalModel* getPointerToExternalModel( );

    //! Overloaded ostream to print class information.
    /*!
     *  Overloaded ostream to print class information; prints what
     *  property models have been set for the vehicle model.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     VehicleModel& vehicleModel );

private:

    /*!
     * Pointer to object that represents vehicle exterior.
     */
    VehicleExternalModel* pointerToExternalModel_;

    /*!
     * Boolean that is true if an external model has been set.
     */
    bool isExternalModelSet_;
};

#endif // VEHICLEMODELS_H

// End of file.
