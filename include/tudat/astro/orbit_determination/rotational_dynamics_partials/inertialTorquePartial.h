/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INERTIALTORQUEPARTIALS_H
#define TUDAT_INERTIALTORQUEPARTIALS_H

#include "tudat/astro/gravitation/secondDegreeGravitationalTorque.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/torquePartial.h"
#include "tudat/astro/orbit_determination/rotational_dynamics_partials/inertiaTensorPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/simulation/environment_setup/rigidBodyProperties.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to compute the partial derivatives of the inertial torque w.r.t. states and parameters
class InertialTorquePartial : public TorquePartial
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param angularVelocityFunction Function that returns body's body-fixed angular velocity vector
     * \param inertiaTensorFunction Function that returns body's inertia tensor
     * \param rigidBodyProperties Rigid body properties of body undergoing torque (inertia tensor depends on gravitational parameter
     * and degree two gravity field coefficients only if these are derived from the gravity field)
     * \param acceleratedBody Name of body undergoing torque
     * \param inertiaTensorDerivativeFunction Function that returns time derivative of body's inertia tensor (empty if
     * this derivative is not included in the inertial torque)
     * \param tidalGravityFieldVariations Tidal gravity field variations of the body, contributing to the inertia tensor derivative
     * \param isTidalGravityFieldVariationInterpolated Per entry of tidalGravityFieldVariations, whether it is evaluated through an
     * interpolator (if so, it does not depend on the propagated rotational state, and is only used for partials w.r.t. parameters)
     * \param tidalLoveNumberPartialInterfaces Objects to compute partials of tidal gravity field variations w.r.t. Love numbers
     * \param tidalLoveNumberRatePartialInterfaces Objects to compute partials of time derivatives of tidal gravity field variations
     * w.r.t. Love numbers
     * \param rotationToBaseFrameFunction Function returning rotation from body-fixed to inertial frame (required if tidal gravity
     * field variations are provided)
     */
    InertialTorquePartial(
            const std::function< Eigen::Vector3d( ) > angularVelocityFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const std::shared_ptr< simulation_setup::RigidBodyProperties > rigidBodyProperties,
            const std::string acceleratedBody,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorDerivativeFunction = std::function< Eigen::Matrix3d( ) >( ),
            const std::vector< std::shared_ptr< gravitation::SolidBodyTideGravityFieldVariations > >& tidalGravityFieldVariations =
                    std::vector< std::shared_ptr< gravitation::SolidBodyTideGravityFieldVariations > >( ),
            const std::vector< bool >& isTidalGravityFieldVariationInterpolated = std::vector< bool >( ),
            const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >& tidalLoveNumberPartialInterfaces =
                    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >( ),
            const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >&
                    tidalLoveNumberRatePartialInterfaces =
                            std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >( ),
            const std::function< Eigen::Quaterniond( ) > rotationToBaseFrameFunction = std::function< Eigen::Quaterniond( ) >( ) ):
        TorquePartial( acceleratedBody, acceleratedBody, basic_astrodynamics::torque_free ),
        angularVelocityFunction_( angularVelocityFunction ), inertiaTensorFunction_( inertiaTensorFunction ),
        getInertiaTensorNormalizationFactor_(
                simulation_setup::getInertiaTensorNormalizationFunctionOfGravityDerivedProperties( rigidBodyProperties ) ),
        bodyGravitationalParameterFunction_(
                simulation_setup::getGravitationalParameterFunctionOfGravityDerivedProperties( rigidBodyProperties ) ),
        inertiaTensorDerivativeFunction_( inertiaTensorDerivativeFunction ), tidalGravityFieldVariations_( tidalGravityFieldVariations ),
        isTidalGravityFieldVariationInterpolated_( isTidalGravityFieldVariationInterpolated ),
        tidalLoveNumberPartialInterfaces_( tidalLoveNumberPartialInterfaces ),
        tidalLoveNumberRatePartialInterfaces_( tidalLoveNumberRatePartialInterfaces ),
        rotationToBaseFrameFunction_( rotationToBaseFrameFunction )
    {
        if( isTidalGravityFieldVariationInterpolated_.size( ) != tidalGravityFieldVariations_.size( ) )
        {
            throw std::runtime_error( "Error in inertial torque partial, inconsistent size of tidal variation interpolation flags" );
        }
        currentPartialDerivativeWrtQuaternion_.setZero( );
    }

    ~InertialTorquePartial( ) {}

    //! Function type computing partials of the tidal part of the inertia tensor and of its time derivative w.r.t. a parameter
    //! (one matrix per parameter entry)
    typedef std::function< void( std::vector< Eigen::Matrix3d >&, std::vector< Eigen::Matrix3d >& ) > TidalInertiaTensorPartialFunction;

    //! Function to retrieve a function computing partials of the tidal part of the inertia tensor (and its derivative) w.r.t. a
    //! double parameter
    /*!
     *  Function to retrieve a function computing partials of the tidal part of the inertia tensor (and its time derivative) w.r.t. a
     *  double parameter (gravitational parameter of deformed or deforming body).
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of partial function and number of parameter entries (0 for no dependency).
     */
    std::pair< TidalInertiaTensorPartialFunction, int > getTidalInertiaTensorPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function to retrieve a function computing partials of the tidal part of the inertia tensor (and its derivative) w.r.t. a
    //! vector parameter
    /*!
     *  Function to retrieve a function computing partials of the tidal part of the inertia tensor (and its time derivative) w.r.t. a
     *  vector parameter (tidal Love numbers).
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of partial function and number of parameter entries (0 for no dependency).
     */
    std::pair< TidalInertiaTensorPartialFunction, int > getTidalInertiaTensorPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function to retrieve the current partials of the tidal part of the inertia tensor w.r.t. the quaternion entries (set by update)
    const std::vector< Eigen::Matrix3d >& getCurrentTidalInertiaTensorPartialsWrtQuaternion( )
    {
        return currentTidalInertiaTensorPartialsWrtQuaternion_;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current torque.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    //! Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the orientation of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. Orientation of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                          const bool addContribution = 1,
                                          const int startRow = 0,
                                          const int startColumn = 0 )
    {
        // Torque depends on orientation only through tidal gravity field variations
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialDerivativeWrtQuaternion_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialDerivativeWrtQuaternion_;
        }
    }

    //! Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body.
    /*!
     *  Function for calculating the partial of the torque w.r.t. the angular velocity of the accelerated body and
     *  adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of torque w.r.t. angular velocity of body
     *  undergoing torque where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtRotationalVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                                 const bool addContribution = 1,
                                                 const int startRow = 0,
                                                 const int startColumn = 3 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialDerivativeWrtAngularVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialDerivativeWrtAngularVelocity_;
        }
    }

    //! Update partial model to current time
    /*!
     * Update partial model to current time, retrieving information from consituent functions
     * \param currentTime Time to which partial is to be updated
     */
    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentAngularVelocityVector_ = angularVelocityFunction_( );
            currentAngularVelocityCrossProductMatrix_ = linear_algebra::getCrossProductMatrix( currentAngularVelocityVector_ );

            // These functions are not defined for all bodies, and are only used by partials that check for their existence
            if( getInertiaTensorNormalizationFactor_ != nullptr )
            {
                currentInertiaTensorNormalizationFactor_ = getInertiaTensorNormalizationFactor_( );
            }
            if( bodyGravitationalParameterFunction_ != nullptr )
            {
                currentGravitationalParameter_ = bodyGravitationalParameterFunction_( );
            }

            currentInertiaTensor_ = inertiaTensorFunction_( );
            currentInverseInertiaTensor_ = currentInertiaTensor_.inverse( );

            currentPartialDerivativeWrtAngularVelocity_ =
                    -linear_algebra::getCrossProductMatrix( currentAngularVelocityVector_ ) * currentInertiaTensor_ +
                    linear_algebra::getCrossProductMatrix( currentInertiaTensor_ * currentAngularVelocityVector_ );

            // Add contribution of inertia tensor derivative (-dI/dt * omega) to torque, including the dependency of dI/dt on omega
            // through tidal gravity field variations
            currentPartialDerivativeWrtQuaternion_.setZero( );
            currentTidalInertiaTensorPartialsWrtQuaternion_.clear( );
            if( inertiaTensorDerivativeFunction_ != nullptr )
            {
                currentPartialDerivativeWrtAngularVelocity_ -= inertiaTensorDerivativeFunction_( );
                for( unsigned int i = 0; i < tidalGravityFieldVariations_.size( ); i++ )
                {
                    if( isTidalGravityFieldVariationInterpolated_.at( i ) )
                    {
                        continue;
                    }
                    std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > coefficientRatePartials =
                            tidalGravityFieldVariations_.at( i )
                                    ->calculateSphericalHarmonicsCorrectionsTimeDerivativePartialWrtAngularVelocity( currentTime, 2 );
                    for( int j = 0; j < 3; j++ )
                    {
                        currentPartialDerivativeWrtAngularVelocity_.block( 0, j, 3, 1 ) -=
                                getInertiaTensorPartialFromDegreeTwoCoefficients( coefficientRatePartials.at( j ).first,
                                                                                  coefficientRatePartials.at( j ).second ) *
                                currentAngularVelocityVector_;
                    }
                }

                // Add contribution of tidal variations to partial w.r.t. orientation (through inertia tensor and its derivative)
                if( tidalGravityFieldVariations_.size( ) > 0 )
                {
                    std::vector< Eigen::Matrix3d > rotationMatrixDerivativesWrtQuaternion( 4 );
                    linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
                            linear_algebra::convertQuaternionToVectorFormat( rotationToBaseFrameFunction_( ) ),
                            rotationMatrixDerivativesWrtQuaternion );
                    std::vector< Eigen::Matrix3d > rotationToBodyFixedFramePartials;
                    for( unsigned int j = 0; j < 4; j++ )
                    {
                        rotationToBodyFixedFramePartials.push_back( rotationMatrixDerivativesWrtQuaternion.at( j ).transpose( ) );
                    }

                    currentTidalInertiaTensorPartialsWrtQuaternion_.assign( 4, Eigen::Matrix3d::Zero( ) );
                    for( unsigned int i = 0; i < tidalGravityFieldVariations_.size( ); i++ )
                    {
                        if( isTidalGravityFieldVariationInterpolated_.at( i ) )
                        {
                            continue;
                        }
                        std::vector< std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > coefficientPartials, coefficientRatePartials;
                        tidalGravityFieldVariations_.at( i )->calculateSphericalHarmonicsCorrectionsAndTimeDerivativesPartialWrtOrientation(
                                currentTime,
                                rotationToBodyFixedFramePartials,
                                std::bind( &InertialTorquePartial::getTidalForcingHessian,
                                           this,
                                           tidalGravityFieldVariations_.at( i )->getDeformedBodyReferenceRadius( ),
                                           std::placeholders::_1,
                                           std::placeholders::_2,
                                           std::placeholders::_3,
                                           std::placeholders::_4 ),
                                2,
                                coefficientPartials,
                                coefficientRatePartials );
                        for( unsigned int j = 0; j < 4; j++ )
                        {
                            const Eigen::Matrix3d inertiaTensorPartial = getInertiaTensorPartialFromDegreeTwoCoefficients(
                                    coefficientPartials.at( j ).first, coefficientPartials.at( j ).second );
                            currentTidalInertiaTensorPartialsWrtQuaternion_.at( j ) += inertiaTensorPartial;
                            currentPartialDerivativeWrtQuaternion_.block( 0, j, 3, 1 ) +=
                                    -currentAngularVelocityCrossProductMatrix_ * inertiaTensorPartial * currentAngularVelocityVector_ -
                                    getInertiaTensorPartialFromDegreeTwoCoefficients( coefficientRatePartials.at( j ).first,
                                                                                      coefficientRatePartials.at( j ).second ) *
                                            currentAngularVelocityVector_;
                        }
                    }
                }

                for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
                {
                    tidalLoveNumberPartialInterfaces_.at( i )->update( currentTime );
                }
                for( unsigned int i = 0; i < tidalLoveNumberRatePartialInterfaces_.size( ); i++ )
                {
                    tidalLoveNumberRatePartialInterfaces_.at( i )->update( currentTime );
                }
            }
            currentTime_ = currentTime;
        }
    }

protected:
    //! Function to reset the Love number partial interfaces to the current time
    void resetCurrentTimeOfMemberObjects( ) override
    {
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->resetCurrentTime( );
        }
        for( unsigned int i = 0; i < tidalLoveNumberRatePartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberRatePartialInterfaces_.at( i )->resetCurrentTime( );
        }
    }

    //! Function to update the parameter partials of the Love number partial interfaces to the current time
    void updateParameterPartialsOfMemberObjects( ) override
    {
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->updateParameterPartials( );
        }
        for( unsigned int i = 0; i < tidalLoveNumberRatePartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberRatePartialInterfaces_.at( i )->updateParameterPartials( );
        }
    }

    //! Function to compute the Hessian of the potential of a single spherical harmonic tidal forcing term
    /*!
     * Function to compute the Hessian (in body-fixed Cartesian coordinates) of the potential of a single (geodesy-normalized, unit
     * coefficient) spherical harmonic term, with gravitational parameter equal to the reference radius
     * \param referenceRadius Reference radius of the spherical harmonic expansion
     * \param degree Degree of term
     * \param order Order of term
     * \param isCosineTerm Boolean denoting whether the cosine (if true) or sine (if false) term is used
     * \param bodyFixedPosition Body-fixed position at which Hessian is to be computed
     * \return Hessian of single spherical harmonic term
     */
    Eigen::Matrix3d getTidalForcingHessian( const double referenceRadius,
                                            const int degree,
                                            const int order,
                                            const bool isCosineTerm,
                                            const Eigen::Vector3d& bodyFixedPosition );

    //! Function to compute partial of the inertia tensor w.r.t. a (geodesy-normalized) degree two cosine or sine coefficient
    /*!
     * Function to compute partial of the inertia tensor w.r.t. a (geodesy-normalized) degree two cosine or sine coefficient
     * \param order Order of coefficient
     * \param isCosineCoefficient Boolean denoting whether partial is w.r.t. cosine (if true) or sine (if false) coefficient
     * \return Partial of inertia tensor w.r.t. coefficient
     */
    Eigen::Matrix3d getInertiaTensorPartialWrtDegreeTwoCoefficient( const int order, const bool isCosineCoefficient );

    //! Function to compute the (partial of the) inertia tensor due to (partials of) degree two coefficients
    /*!
     * Function to compute the (partial of the) inertia tensor due to (partials of) geodesy-normalized coefficients, or of its time
     * derivative due to (partials of) coefficient rates (first row is degree two, first column is order zero)
     * \param cosineCoefficientPartials (Partials of) cosine coefficients or their rates
     * \param sineCoefficientPartials (Partials of) sine coefficients or their rates
     * \return (Partial of) inertia tensor or its time derivative
     */
    Eigen::Matrix3d getInertiaTensorPartialFromDegreeTwoCoefficients( const Eigen::MatrixXd& cosineCoefficientPartials,
                                                                      const Eigen::MatrixXd& sineCoefficientPartials );

    //! Function to compute partial of torque w.r.t. a parameter, through the tidal part of the inertia tensor and its derivative
    /*!
     * Function to compute partial of torque w.r.t. a parameter, through the tidal part of the inertia tensor and its time derivative
     * \param partialMatrix Computed partial of torque w.r.t. parameter (returned by reference)
     * \param tidalPartialFunction Function computing partials of tidal part of inertia tensor and its time derivative
     */
    void wrtParameterThroughTidalInertiaTensor( Eigen::MatrixXd& partialMatrix,
                                                const TidalInertiaTensorPartialFunction& tidalPartialFunction );

    //! Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. non-mode coupled tidal Love numbers
    /*!
     * Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. non-mode coupled tidal Love numbers
     * \param inertiaTensorPartials Partials of inertia tensor, per parameter entry (returned by reference)
     * \param inertiaTensorDerivativePartials Partials of inertia tensor derivative, per parameter entry (returned by reference)
     * \param coefficientPartialFunctions Functions returning partials of coefficients w.r.t. Love numbers, per order
     * \param coefficientRatePartialFunctions Functions returning partials of coefficient rates w.r.t. Love numbers, per order
     * \param orders Orders of Love numbers
     * \param sumOrders Boolean denoting whether a single Love number is used for all orders
     * \param parameterSize Size of parameter
     */
    void tidalInertiaTensorPartialsWrtLoveNumbers(
            std::vector< Eigen::Matrix3d >& inertiaTensorPartials,
            std::vector< Eigen::Matrix3d >& inertiaTensorDerivativePartials,
            const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientPartialFunctions,
            const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >&
                    coefficientRatePartialFunctions,
            const std::vector< int >& orders,
            const bool sumOrders,
            const int parameterSize );

    //! Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. mode-coupled tidal Love numbers
    /*!
     * Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. mode-coupled tidal Love numbers
     * \param inertiaTensorPartials Partials of inertia tensor, per parameter entry (returned by reference)
     * \param inertiaTensorDerivativePartials Partials of inertia tensor derivative, per parameter entry (returned by reference)
     * \param coefficientPartialFunctions Functions returning partials of forcing w.r.t. Love numbers, per Love number
     * \param coefficientRatePartialFunctions Functions returning partials of forcing rates w.r.t. Love numbers, per Love number
     * \param responseIndices Index of response degree and order (in responseDegreeOrders) per Love number
     * \param responseDegreeOrders Degrees and orders of responses
     * \param parameterSize Size of parameter
     */
    void tidalInertiaTensorPartialsWrtModeCoupledLoveNumbers(
            std::vector< Eigen::Matrix3d >& inertiaTensorPartials,
            std::vector< Eigen::Matrix3d >& inertiaTensorDerivativePartials,
            const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientPartialFunctions,
            const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >&
                    coefficientRatePartialFunctions,
            const std::vector< int >& responseIndices,
            const std::vector< std::pair< int, int > >& responseDegreeOrders,
            const int parameterSize );

    //! Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. gravitational parameter of the body
    /*!
     * Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. gravitational parameter of the body
     * undergoing the torque. The tidal contributions are independent of this parameter, while the rest of the inertia tensor (and
     * its derivative) is proportional to it, so that these partials are minus the tidal contributions divided by the parameter.
     * \param inertiaTensorPartials Partial of inertia tensor (returned by reference)
     * \param inertiaTensorDerivativePartials Partial of inertia tensor derivative (returned by reference)
     */
    void tidalInertiaTensorPartialsWrtGravitationalParameter( std::vector< Eigen::Matrix3d >& inertiaTensorPartials,
                                                              std::vector< Eigen::Matrix3d >& inertiaTensorDerivativePartials );

    //! Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. gravitational parameter of a body
    //! causing tidal deformation
    /*!
     * Function to compute partials of tidal part of inertia tensor (and its derivative) w.r.t. gravitational parameter of a body
     * causing tidal deformation
     * \param inertiaTensorPartials Partial of inertia tensor (returned by reference)
     * \param inertiaTensorDerivativePartials Partial of inertia tensor derivative (returned by reference)
     * \param variationAndBodyIndices List of indices of tidal variations (in tidalGravityFieldVariations_) and of the body causing
     * deformation in these variations
     */
    void tidalInertiaTensorPartialsWrtDeformingBodyGravitationalParameter(
            std::vector< Eigen::Matrix3d >& inertiaTensorPartials,
            std::vector< Eigen::Matrix3d >& inertiaTensorDerivativePartials,
            const std::vector< std::pair< int, int > >& variationAndBodyIndices );

    //! Function to compute partial of torque w.r.t. mean moment of inertia
    /*!
     * Function to compute partial of torque w.r.t. mean moment of inertia
     * \param momentOfInertiaPartial Computed partial of torque w.r.t. mean moment of inertia (returned by reference)
     */
    void wrtMeanMomentOfInertia( Eigen::MatrixXd& momentOfInertiaPartial );

    //! Function to compute partial of torque w.r.t. gravitational parameter
    /*!
     * Function to compute partial of torque w.r.t. gravitational parameter
     * \param gravitationalParameterPartial Computed partial of torque w.r.t. gravitational parameter (returned by reference)
     * \param tidalPartialFunction Function computing partials of tidal part of inertia tensor and its time derivative (empty if
     * none)
     */
    void wrtGravitationalParameter( Eigen::MatrixXd& gravitationalParameterPartial,
                                    const TidalInertiaTensorPartialFunction& tidalPartialFunction );

    //! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic cosine coefficients
     * (returned by reference)
     * \param c20Index for degree=2,order=0 coefficient
     * \param c21Index for degree=2,order=1 coefficient
     * \param c22Index for degree=2,order=2 coefficient
     */
    void wrtCosineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                              const int c20Index,
                                                              const int c21Index,
                                                              const int c22Index );

    //! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
    /*!
     * Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
     * \param sphericalHarmonicCoefficientPartial Computed partial of torque w.r.t. spherical harmonic sine coefficients
     * (returned by refe
     * \param s21Index for degree=2,order=1 coefficient
     * \param s22Index for degree=2,order=2 coefficient)
     */
    void wrtSineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                            const int s21Index,
                                                            const int s22Index );

    //! Function returning body angular velocity vector in body fixed frame.
    std::function< Eigen::Vector3d( ) > angularVelocityFunction_;

    //! Function returning body inertia tensor
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    //! Function the inertia tensor normalization factor
    std::function< double( ) > getInertiaTensorNormalizationFactor_;

    //!  Function returning body gravitational parameter
    std::function< double( ) > bodyGravitationalParameterFunction_;

    //! Function returning time derivative of body inertia tensor (empty if not included in torque)
    std::function< Eigen::Matrix3d( ) > inertiaTensorDerivativeFunction_;

    //! Tidal gravity field variations of the body, contributing to the inertia tensor derivative
    std::vector< std::shared_ptr< gravitation::SolidBodyTideGravityFieldVariations > > tidalGravityFieldVariations_;

    //! Per entry of tidalGravityFieldVariations_, whether it is evaluated through an interpolator (no rotational state dependency)
    std::vector< bool > isTidalGravityFieldVariationInterpolated_;

    //! Objects to compute partials of tidal gravity field variations w.r.t. Love numbers
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > tidalLoveNumberPartialInterfaces_;

    //! Objects to compute partials of time derivatives of tidal gravity field variations w.r.t. Love numbers
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > tidalLoveNumberRatePartialInterfaces_;

    //! Function returning rotation from body-fixed to inertial frame
    std::function< Eigen::Quaterniond( ) > rotationToBaseFrameFunction_;

    //! Current partials of tidal part of inertia tensor w.r.t. quaternion entries (empty if no tidal variations)
    std::vector< Eigen::Matrix3d > currentTidalInertiaTensorPartialsWrtQuaternion_;

    //! Current partial derivative of torque w.r.t. quaternion (due to tidal gravity field variations)
    Eigen::Matrix< double, 3, 4 > currentPartialDerivativeWrtQuaternion_;

    //! Cache for computing Hessians of spherical harmonic tidal forcing terms
    basic_mathematics::SphericalHarmonicsCache tidalForcingHessianCache_;

    //! Current angular velocity vector
    Eigen::Vector3d currentAngularVelocityVector_;

    //! Current cross-product matrix from angular velocity vector
    Eigen::Matrix3d currentAngularVelocityCrossProductMatrix_;

    //! Current inertia tensor
    Eigen::Matrix3d currentInertiaTensor_;

    //! Current inverse inertia tensor
    Eigen::Matrix3d currentInverseInertiaTensor_;

    //! Current partial derivative w.r.t. angular velocity vector
    Eigen::Matrix3d currentPartialDerivativeWrtAngularVelocity_;

    //! Current inertia tensor normalization
    double currentInertiaTensorNormalizationFactor_;

    //! Current gravitational parameter
    double currentGravitationalParameter_;
};

}  // namespace acceleration_partials

}  // namespace tudat

#endif  // TUDAT_INERTIALTORQUEPARTIALS_H
