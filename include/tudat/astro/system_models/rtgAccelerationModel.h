/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Precise computation of acceleration due to uniform ring or disk, Toshio Fukushima (2010), Celestial Mechanics
 *          and Dynamical Astronomy, 108:339–356.
 */

#ifndef TUDAT_RTGACCELERATIONMODEL_H
#define TUDAT_RTGACCELERATIONMODEL_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <functional>
#include <iostream>

#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{
namespace system_models
{

class RTGAccelerationModel : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
protected:
    //! Typedef for a position-returning function.
    //! (not necessary?) typedef std::function< void( Eigen::Vector3d& ) > StateFunction;

public:
    //! Constructor taking constant parameters for force decay process and rotation/mass functions for bodies.
    /*!
     *  \param bodyFixedForceVectorAtReferenceEpoch Pointer to vector 1x3 with force values along body-fixed frame axes at reference epoch.
     *  \param decayScaleFactor Scale factor value for modelling of force decay process.
     *  \param referenceEpoch Reference epoch for modelling of force decay process.
     *  \param rotationFromBodyFixedToIntegrationFrameFunction Function returning (undergoing) body rotational ephemeris.
     *  \param rotationFromBodyFixedToIntegrationFrameFunction Function returning (undergoing) body current mass.
     *  \return RTG acceleration model.
     */
    RTGAccelerationModel(
            const Eigen::Vector3d& bodyFixedForceVectorAtReferenceEpoch,
            const double decayScaleFactor,
            const double referenceEpoch,
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction,
            const std::function< double( ) > bodyMassFunction):

        bodyFixedForceVectorAtReferenceEpoch_( bodyFixedForceVectorAtReferenceEpoch ),
        forceVectorMagnitudeAtReferenceEpoch_(bodyFixedForceVectorAtReferenceEpoch.norm( )),
        bodyFixedForceUnitVectorAtReferenceEpoch_(bodyFixedForceVectorAtReferenceEpoch/bodyFixedForceVectorAtReferenceEpoch.norm( )),
        decayScaleFactor_( decayScaleFactor ),
        referenceEpoch_( referenceEpoch ),
        rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
        bodyMassFunction_( bodyMassFunction )

    { }

    ~RTGAccelerationModel( ) { }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class members of this class.
     * The potential and laplacian of potential are only updated if the associated flags indicate so.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) && ( currentTime == currentTime) )
        {
            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
            currentTimeDelta_ = currentTime - referenceEpoch_;
            currentDecayTerm_ = std::exp(-decayScaleFactor_ * currentTimeDelta_);
            currentBodyFixedForceVector_ = bodyFixedForceVectorAtReferenceEpoch_ * currentDecayTerm_;
            currentAcceleration_ = rotationToIntegrationFrame_ * currentBodyFixedForceVector_ / bodyMassFunction_();
        }
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

    double getCurrentTimeDelta( ) const
    {
        return currentTimeDelta_;
    }

    Eigen::Vector3d getCurrentBodyFixedForceVector( ) const
    {
        return currentBodyFixedForceVector_;
    }

    void resetForceVectorAtReferenceEpoch(const Eigen::Vector3d& newForceVectorAtReferenceEpoch)
    {
        bodyFixedForceVectorAtReferenceEpoch_ = newForceVectorAtReferenceEpoch;
    }

    void resetForceMagnitudeAtReferenceEpoch(const double newForceVectorMagnitudeAtReferenceEpoch)
    {
        forceVectorMagnitudeAtReferenceEpoch_ = newForceVectorMagnitudeAtReferenceEpoch;
        updateBodyFixedForceVectorAtReferenceEpoch( );
    }

    void updateBodyFixedForceVectorAtReferenceEpoch( )
    {
        bodyFixedForceVectorAtReferenceEpoch_ = forceVectorMagnitudeAtReferenceEpoch_ * bodyFixedForceUnitVectorAtReferenceEpoch_;
    }

    void updateBodyFixedForceUnitVectorAtReferenceEpoch( )
    {
        bodyFixedForceUnitVectorAtReferenceEpoch_ = bodyFixedForceUnitVectorAtReferenceEpoch_ / bodyFixedForceUnitVectorAtReferenceEpoch_.norm();
    }

    Eigen::Vector3d getbodyFixedForceVectorAtReferenceEpoch( ) const
    {
        return bodyFixedForceVectorAtReferenceEpoch_;
    }

    double getForceVectorMagnitudeAtReferenceEpoch( ) const
    {
        return forceVectorMagnitudeAtReferenceEpoch_;
    }

    Eigen::Vector3d getBodyFixedForceUnitVectorAtReferenceEpoch( ) const
    {
        return bodyFixedForceUnitVectorAtReferenceEpoch_;
    }

    double getCurrentDecayTerm( ) const
    {
        return currentDecayTerm_;
    }

    double evaluateBodyMassFunction( ) const
    {
        return bodyMassFunction_( );
    }


private:

    //! Force vector in body-fixed frame at reference epoch
    Eigen::Vector3d bodyFixedForceVectorAtReferenceEpoch_;

    //! Force vector magnitude at reference epoch
    double forceVectorMagnitudeAtReferenceEpoch_;

    //! Force vector in body-fixed frame at reference epoch
    Eigen::Vector3d bodyFixedForceUnitVectorAtReferenceEpoch_;

    //! Scale Factor for force decay process
    double decayScaleFactor_;

    //! Reference epoch for modelled force decay process
    double referenceEpoch_;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Body mass function
    const std::function< double() > bodyMassFunction_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current acceleration in frame fixed to body undergoing acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    //! Delta between current time and reference epoch
    double currentTimeDelta_;

    //! Exponential term (handy for partials computation)
    double currentDecayTerm_;

    //! Body fixed force vector at current epoch incl effect of decay law
    Eigen::Vector3d currentBodyFixedForceVector_;

};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_RTGACCELERATIONMODEL_H
