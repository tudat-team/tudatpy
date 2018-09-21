/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H
#define TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H


#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"

namespace tudat
{

namespace gravitation
{

//! Class to compute the gravitational torque on a body, due to a point mass perturber
/*!
 *  Class to compute the gravitational torque on a body, due to a point mass perturber, and an arbitrary dergee expansion
 *  of the spherical harmonic gravity field of the body on which the torque is acting. This class uses the
 *  SphericalHarmonicsGravitationalAccelerationModel class to compute the acceleration that causes the torque
 */
class SphericalHarmonicGravitationalTorqueModel: public basic_astrodynamics::TorqueModel
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param sphericalHarmonicAcceleration Spherical harmonic acceleration that the body that undergoes the torque exerts
     * on the body that exerts the torque
     * \param rotationToBodyUndergoingTorque Function returning the rotation from an inertial frame to the body-fixed frame of
     * the body undergoing the torque
     * \param perturberMassFunction Function returning the mass of the body exerting the torque.
     */
    SphericalHarmonicGravitationalTorqueModel(
            const std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration,
            const std::function< Eigen::Quaterniond( ) > rotationToBodyUndergoingTorque,
            const std::function< double( ) > perturberMassFunction ):
        sphericalHarmonicAcceleration_( sphericalHarmonicAcceleration ),
        rotationToBodyUndergoingTorque_( rotationToBodyUndergoingTorque ),
        perturberMassFunction_( perturberMassFunction ){ }

    //! Get gravitational torque.
    /*!
     * Returns the gravitational torque. All data required for the computation is taken
     * from member variables, which are set to their latest values by the last call of the
     * updateMembers function.
     * \return Current gravitational torque.
     * \sa updateMembers().
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the gravitational torque model.
    /*!
     * Updates member variables used by the gravitational accfeleration model.
     * Pointers and function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        sphericalHarmonicAcceleration_->updateMembers( currentTime );

        currentTorque_ = -perturberMassFunction_( ) *
                ( ( sphericalHarmonicAcceleration_->getCurrentRelativePosition( ) ).cross(
                      sphericalHarmonicAcceleration_->getAccelerationInBodyFixedFrame( ) ) );
    }

    //! Function to retrieve spherical harmonic acceleration
    /*!
     *  Function to retrieve spherical harmonic acceleration that the body that undergoes the torque exerts on the body that
     *  exerts the torque
     *  \return Spherical harmonic acceleration that the body that undergoes the torque exerts on the body that
     *  exerts the torque
     */
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > getSphericalHarmonicAcceleration( )
    {
        return sphericalHarmonicAcceleration_;
    }


protected:

private:

    //! Spherical harmonic acceleration that the body that undergoes the torque exerts on the body that exerts the torque
    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration_;

    //! Function returning the rotation from an inertial frame to the body-fixed frame of the body undergoing the torque
    std::function< Eigen::Quaterniond( ) > rotationToBodyUndergoingTorque_;

    //! Function returning the mass of the body exerting the torque.
    std::function< double( ) > perturberMassFunction_;

    //! Torque, as computed by last call to updateMembers function
    Eigen::Vector3d currentTorque_;
};

}

}

#endif // TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H
