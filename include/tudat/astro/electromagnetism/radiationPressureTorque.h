/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSURETORQUE_H
#define TUDAT_RADIATIONPRESSURETORQUE_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling radiation pressure torque. Radiation pressure torque a target due to electromagnetic
 * radiation from a source.
 */
class IsotropicPointSourceRadiationPressureTorque: public basic_astrodynamics::TorqueModel
{
public:

    IsotropicPointSourceRadiationPressureTorque(
        const std::shared_ptr< IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration,
        const std::function< Eigen::Vector3d( ) > centerOfMassFunction ):
        radiationPressureAcceleration_( radiationPressureAcceleration ),
        centerOfMassFunction_( centerOfMassFunction )
    {
        radiationPressureAcceleration_->getTargetModel( )->enableTorqueComputation( centerOfMassFunction );
    }

    ~IsotropicPointSourceRadiationPressureTorque( ){ }

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime) override
    {
        radiationPressureAcceleration_->updateMembers( currentTime );
        currentTorque_ = radiationPressureAcceleration_->getTargetModel( )->getCurrentRadiationPressureTorque( );
        if( !radiationPressureAcceleration_->getTargetModel( )->forceFunctionRequiresLocalFrameInputs( ) )
        {
            currentTorque_ = radiationPressureAcceleration_->getTargetRotationFromLocalToGlobalFrameFunction( )( ).inverse( ) * currentTorque_;
        }
    }

    void resetAccelerationModel( const std::shared_ptr< IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration )
    {
        radiationPressureAcceleration_ = radiationPressureAcceleration;
    }

    Eigen::Vector3d getTorque( ) override
    {
        return currentTorque_;
    }


    virtual void resetCurrentTime( ) override
    {
        currentTime_ = TUDAT_NAN;
        radiationPressureAcceleration_->resetCurrentTime( );
    }

protected:

    std::shared_ptr< IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration_;

    const std::function< Eigen::Vector3d( ) > centerOfMassFunction_;

    Eigen::Vector3d currentTorque_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSURETORQUE_H
