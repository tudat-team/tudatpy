/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_GRAVITY_DEFORMATION_MODEL_H
#define TUDAT_GRAVITY_DEFORMATION_MODEL_H

#include <vector>
#include <map>
#include <unordered_map>

#include <memory>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{

namespace gravitation
{

class SphericalHarmonicsGravityField;

}  // namespace gravitation

namespace basic_astrodynamics
{

//! Base class for gravity deformation models
/*!
 * Base class for gravity deformation models. Derived classes should contain
 * implementations to perform calculations of gravity field deformations.
 */
class GravityDeformationModel
{
public:
    //! Constructor.
    GravityDeformationModel( ): currentTime_( TUDAT_NAN ), currentDeformation_( Eigen::VectorXd::Zero( 5 ) ) {}

    //! Virtual destructor.
    /*!
     * Virtual destructor, necessary to ensure that derived class destructors get called correctly.
     */
    virtual ~GravityDeformationModel( ) {}

    //! Update member variables used by the gravity deformation model.
    /*!
     * Updates member variables used by the gravity deformation model. In the case of acceleration models
     * containing varying parameters, function-pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function-pointers and updates member variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used by the getAcceleration() function.
     *
     * N.B.: This pure virtual function must be overridden by derived classes!
     * \param currentTime Time at which acceleration model is to be updated.
     */
    virtual void updateMembers( const double currentTime = TUDAT_NAN ) = 0;

    Eigen::VectorXd& getDeformationReference( )
    {
        return currentDeformation_;
    }

    Eigen::VectorXd getDeformation( )
    {
        return currentDeformation_;
    }

    //! Get the current derivative of the gravity-deformation state.
    const Eigen::VectorXd& getCurrentStateDerivative( ) const
    {
        return currentDeformation_;
    }

    void getAccelerationByReference( Eigen::VectorXd& deformation ) const
    {
        deformation = currentDeformation_;
    }

    void addCurrentDeformation( Eigen::VectorXd& deformation ) const
    {
        if( deformation.size( ) != currentDeformation_.size( ) )
        {
            throw std::runtime_error( "Error when adding current gravity deformation, inconsistent sizes." );
        }
        deformation += currentDeformation_;
    }

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the deformation gravity model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

protected:
    //! Previous time to which the gravity deformation model was updated.
    double currentTime_;

    Eigen::VectorXd currentDeformation_;

protected:
private:
};

//! Typedef for the gravity deformation model map.
typedef std::map< std::string, std::vector< std::shared_ptr< GravityDeformationModel > > > GravityDeformationModelMap;

//! Class for Maxwell rheology gravity deformation model.
/*!
 * This class implements a gravity field deformation model assuming a Maxwell rheology.
 */
class MaxwellGravityDeformationModel : public GravityDeformationModel
{
public:
    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector6d& ) > StateFunction;

    //! Create a Maxwell model linked to the deforming body's gravity field.
    /*!
     * The gravity field remains the owner of the variation-free coefficient baseline. At each
     * update this model subtracts that baseline from the current total coefficients and uses only
     * the resulting variation in the Maxwell evolution equation.
     * \param stateOfDeformingBodyFunction Function returning the deforming body's Cartesian state.
     * \param perturbingBody Names of bodies producing the deformation.
     * \param maxwellRelaxationTime Maxwell relaxation time.
     * \param globalRelaxationTime Global relaxation time.
     * \param gravityFieldModel Spherical-harmonic gravity field of the deforming body.
     * \param gravitationalParameterPerturbingBody Functions returning perturber gravitational parameters.
     * \param angularVelocityDeformingBody Function returning deforming-body angular velocity.
     * \param angularVelocityDerivativeDeformingBody Function returning its angular-velocity derivative.
     * \param k2 Degree-two Love number.
     * \param stateOfPerturbingBodyFunction Functions returning perturber Cartesian states.
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixed to integration frame.
     * \param rotationToLocalFrameDerivativeFunction Derivative of the integration-to-body-fixed rotation.
     * \param includeOrder1 Whether degree-two, order-one deformation terms are included.
     * \param includeCentrifugalPotential Whether the centrifugal potential contributes to equilibrium.
     */
    MaxwellGravityDeformationModel(
            const StateFunction stateOfDeformingBodyFunction,
            const std::vector< std::string > perturbingBody,
            const double maxwellRelaxationTime,
            const double globalRelaxationTime,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldModel,
            const std::vector< std::function< double( ) > > gravitationalParameterPerturbingBody,
            const std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody,
            const std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody,
            const double k2,
            const std::vector< StateFunction > stateOfPerturbingBodyFunction = {},
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    []( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const std::function< Eigen::Matrix3d( ) > rotationToLocalFrameDerivativeFunction = []( ) { return Eigen::Matrix3d::Zero( ); },
            const bool includeOrder1 = true,
            const bool includeCentrifugalPotential = false );

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class
     * members of this class.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN );

    void updateEquilibriumDeformation( const double currentTime = TUDAT_NAN );

    std::vector< std::string > getPerturbingBody( ) const
    {
        return perturbingBody_;
    }

    //! Get the current equilibrium gravity coefficients in unnormalized form.
    const Eigen::VectorXd& getCurrentEquilibriumCoefficients( ) const
    {
        return equilibriumCoefficients_;
    }

    //! Get the current time derivative of the equilibrium gravity coefficients in unnormalized form.
    const Eigen::VectorXd& getCurrentEquilibriumCoefficientDerivative( ) const
    {
        return derivativeEquilibriumCoefficients_;
    }

protected:
private:
    Eigen::Vector6d stateOfDeformingBody_;

    const std::vector< std::string > perturbingBody_;

    std::vector< Eigen::Vector6d > stateOfPerturbingBody_;

    const double maxwellRelaxationTime_;

    const double globalRelaxationTime_;

    //! Gravity field providing both current and variation-free coefficient values.
    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldModel_;

    const std::vector< std::function< double( ) > > gravitationalParameterPerturbingBody_;

    //! Love number k2
    const double k2_;

    Eigen::VectorXd equilibriumCoefficients_;

    Eigen::VectorXd derivativeEquilibriumCoefficients_;

    //! Current unnormalised gravity-coefficient variation relative to the gravity field baseline.
    Eigen::VectorXd currentCoefficientVariation_;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    std::function< Eigen::Matrix3d( ) > rotationToBodyFixedDerivativeFunction_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    std::vector< Eigen::Vector3d > currentRelativePosition_;

    std::vector< Eigen::Vector3d > currentRelativeVelocity_;

    std::vector< Eigen::Vector6d > currentInertialRelativeState_;

    //! Function returning the state of the body undergoing deformation
    StateFunction stateOfDeformingBodyFunction_;

    //! Function returning the state of the body causing the deformation
    std::vector< StateFunction > stateOfPerturbingBodyFunction_;

    //! Current body-fixed longitude of the perturbing body
    std::vector< double > currentLongitude_;

    //! Current body-fixed latitude of the perturbing body
    std::vector< double > currentLatitude_;

    std::vector< double > currentLongitudeDerivative_;
    std::vector< double > currentLatitudeDerivative_;

    const bool includeOrder1_;
    const bool includeCentrifugalPotential_;

    std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody_;
    std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif  // TUDAT_GRAVITY_DEFORMATION_MODEL_H
