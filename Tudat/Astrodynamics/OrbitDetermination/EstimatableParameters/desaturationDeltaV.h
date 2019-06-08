/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DESATURATIONDELTAV_H
#define TUDAT_DESATURATIONDELTAV_H

#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's momentum wheel desaturation deltaVs.
/*!
 *  Interface class for estimation of a body's momentum wheel desaturation deltaVs.
 *  Interfaces the estimation with the deltaV values of a momentum wheel desaturation thrust
 *  acceleration object.
 */
class DesaturationDeltaV: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param accelerationModel MomentumWheelDesaturationThrustAcceleration object of which
     *  the deltaV values are a property.
     *  \param associatedBody Name of body undergoing the momentum wheel desaturation acceleration.
     */
    DesaturationDeltaV(
            const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > accelerationModel,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( desaturation_delta_v_values, associatedBody ),
        accelerationModel_( accelerationModel )
    {
        numberOfDeltaVBlocks_  = accelerationModel_->getDeltaVValues( ).size( );
    }

    //! Destructor.
    ~DesaturationDeltaV( ) { }

    //! Get deltaV values for each momentum wheel desaturation maneuver.
    /*!
     *  Get deltaV values for each momentum wheel desaturation maneuver.
     *  \return Vector containing the deltaV values.
     */
    Eigen::VectorXd getParameterValue( )
    {
        std::vector< Eigen::Vector3d > deltaVList = accelerationModel_->getDeltaVValues( );
        Eigen::VectorXd deltaVVector = Eigen::VectorXd( 3 * numberOfDeltaVBlocks_ );

        for( int i = 0; i < numberOfDeltaVBlocks_; i++ )
        {
           deltaVVector.segment( i * 3, 3 ) = deltaVList.at( i );
        }
        return deltaVVector;
    }

    //! Reset deltaV values for all the momentum wheel desaturation maneuvers.
    /*!
     *  Reset deltaV values for all the momentum wheel desaturation maneuvers.
     *  \param parameterValue New deltaV values for the desaturation maneuvers.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {

        std::vector< Eigen::Vector3d > deltaVList;
        for( int i = 0; i < numberOfDeltaVBlocks_; i++ )
        {
            deltaVList.push_back( parameterValue.segment( i * 3, 3 ) );
        }

        accelerationModel_->setDeltaVValues( deltaVList );

    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value (number of desaturation maneuvers times three components
     *  for this parameter)
     */
    int getParameterSize( )
    {
        return 3 * numberOfDeltaVBlocks_;
    }

protected:

private:

    //! Momentum wheel desaturation thrust acceleration.
    const std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > accelerationModel_;

    //! Number of momemtum wheel desaturation maneuvers.
    int numberOfDeltaVBlocks_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_DESATURATIONDELTAV_H
