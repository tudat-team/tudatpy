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
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{
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
private:
    //! Typedef for coefficient-matrix-returning function.
    typedef std::function< Eigen::MatrixXd( ) > CoefficientMatrixReturningFunction;

public:
    //! Typedef for a position-returning function.
    typedef std::function< void( Eigen::Vector6d& ) > StateFunction;

    //! Constructor taking position-functions for bodies, and constant parameters of spherical
    //! harmonics expansion.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, gravitational-parameter functions and equatorial radius of the
     * body exerting the acceleration, constant coefficient matrices for the spherical harmonics
     * expansion, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * provided functions to retrieve the current gravitational parameters, equatorial radius and
     * coefficient matrices. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter Function returning the current gravitational parameter [m^3 s^-2].
     * \param anEquatorialRadius A (constant) equatorial radius [m].
     * \param aCosineHarmonicCoefficientMatrix A (constant) cosine harmonic coefficient matrix.
     * \param aSineHarmonicCoefficientMatrix A (constant) sine harmonic coefficient matrix.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    MaxwellGravityDeformationModel(
            const StateFunction stateOfDeformingBodyFunction,
            const std::vector< std::string > perturbingBody,
            const double maxwellRelaxationTime,
            const double globalRelaxationTime,
            const std::function< double( ) > gravitationalParameterDeformingBody,
            const std::vector< std::function< double( ) > > gravitationalParameterPerturbingBody,
            const double referenceRadius,
            const std::function< Eigen::Vector3d( ) > angularVelocityDeformingBody,
            const std::function< Eigen::Vector3d( ) > angularVelocityDerivativeDeformingBody,
            const double k2,
            CoefficientMatrixReturningFunction cosineCoefficients,
            CoefficientMatrixReturningFunction sineCoefficients,
            const std::vector< StateFunction > stateOfPerturbingBodyFunction = {},
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                    []( ) { return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const std::function< Eigen::Matrix3d( ) > rotationToLocalFrameDerivativeFunction = []( ) { return Eigen::Matrix3d::Zero( ); },
            const Eigen::VectorXd staticCoefficients = Eigen::VectorXd::Zero( 3 ),
            const bool includeOrder1 = true,
            const bool includeCentrifugalPotential = false ):
        GravityDeformationModel( ), perturbingBody_( perturbingBody ), maxwellRelaxationTime_( maxwellRelaxationTime ),
        globalRelaxationTime_( globalRelaxationTime ), gravitationalParameterDeformingBody_( gravitationalParameterDeformingBody ),
        gravitationalParameterPerturbingBody_( gravitationalParameterPerturbingBody ), referenceRadius_( referenceRadius ), k2_( k2 ),
        staticCoefficients_( staticCoefficients ), getCosineHarmonicsCoefficients( cosineCoefficients ),
        getSineHarmonicsCoefficients( sineCoefficients ),
        rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
        rotationToBodyFixedDerivativeFunction_( rotationToLocalFrameDerivativeFunction ),
        stateOfDeformingBodyFunction_( stateOfDeformingBodyFunction ), stateOfPerturbingBodyFunction_( stateOfPerturbingBodyFunction ),
        includeOrder1_( includeOrder1 ), includeCentrifugalPotential_( includeCentrifugalPotential ),
        angularVelocityDeformingBody_( angularVelocityDeformingBody ),
        angularVelocityDerivativeDeformingBody_( angularVelocityDerivativeDeformingBody )
    {
        unsigned int numberPerturbingBodies = perturbingBody.size( );
        stateOfPerturbingBody_.resize( numberPerturbingBodies );
        currentRelativePosition_.resize( numberPerturbingBodies );
        currentRelativeVelocity_.resize( numberPerturbingBodies );
        currentInertialRelativeState_.resize( numberPerturbingBodies );
        currentLongitude_.resize( numberPerturbingBodies );
        currentLatitude_.resize( numberPerturbingBodies );
        currentLongitudeDerivative_.resize( numberPerturbingBodies );
        currentLatitudeDerivative_.resize( numberPerturbingBodies );

        // Initialise nominal coefficient values
        nominalCoefficients_ = Eigen::VectorXd::Zero( 5 );
        nominalCoefficients_[ 0 ] = getCosineHarmonicsCoefficients( )( 2, 0 );
        nominalCoefficients_[ 1 ] = getCosineHarmonicsCoefficients( )( 2, 1 );
        nominalCoefficients_[ 2 ] = getCosineHarmonicsCoefficients( )( 2, 2 );
        nominalCoefficients_[ 3 ] = getSineHarmonicsCoefficients( )( 2, 1 );
        nominalCoefficients_[ 4 ] = getSineHarmonicsCoefficients( )( 2, 2 );
        std::cout << "original nominal coefficients: " << nominalCoefficients_.transpose( ) << std::endl;

        equilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );
        derivativeEquilibriumCoefficients_ = Eigen::VectorXd::Zero( 5 );

        // Tranform to **unnormalised** coefficients
        staticCoefficients_[ 0 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        staticCoefficients_[ 1 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
        staticCoefficients_[ 2 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
        staticCoefficients_[ 3 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
        staticCoefficients_[ 4 ] *= basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        std::cout << "done creating MaxwellDeformation model" << std::endl;
    }

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

    const std::function< double( ) > gravitationalParameterDeformingBody_;

    const std::vector< std::function< double( ) > > gravitationalParameterPerturbingBody_;

    const double referenceRadius_;

    //! Love number k2
    const double k2_;

    Eigen::VectorXd equilibriumCoefficients_;

    Eigen::VectorXd derivativeEquilibriumCoefficients_;

    Eigen::VectorXd nominalCoefficients_;

    Eigen::VectorXd staticCoefficients_;

    //! Matrix of cosine coefficients.
    /*!
     * Matrix containing coefficients of cosine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd cosineHarmonicCoefficients;

    //! Matrix of sine coefficients.
    /*!
     * Matrix containing coefficients of sine terms for spherical harmonics expansion.
     */
    Eigen::MatrixXd sineHarmonicCoefficients;

    //! Pointer to function returning cosine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the cosine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getCosineHarmonicsCoefficients;

    //! Pointer to function returning sine harmonics coefficients matrix.
    /*!
     * Pointer to function that returns the current coefficients of the sine terms of the
     * spherical harmonics expansion.
     */
    const CoefficientMatrixReturningFunction getSineHarmonicsCoefficients;

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
