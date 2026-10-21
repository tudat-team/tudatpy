/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/inertialTorquePartial.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicPartialFunctions.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case gravitational_parameter: {
                if( bodyGravitationalParameterFunction_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error when getting partial of inertial torque w.r.t. gravitational parameter, gravitational parameter "
                            "function not found." );
                }
                partialFunction =
                        std::make_pair( std::bind( &InertialTorquePartial::wrtGravitationalParameter, this, std::placeholders::_1 ), 1 );

                break;
            }
            case mean_moment_of_inertia: {
                if( getInertiaTensorNormalizationFactor_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error when getting partial of inertial torque w.r.t. mean moment of inertia, inertia tensor normalization "
                            "function not found." );
                }
                partialFunction =
                        std::make_pair( std::bind( &InertialTorquePartial::wrtMeanMomentOfInertia, this, std::placeholders::_1 ), 1 );

                break;
            }
            default:
                break;
        }
    }
    else if( parameter->getParameterName( ).first == gravitational_parameter && inertiaTensorDerivativeFunction_ != nullptr )
    {
        // Check if parameter is gravitational parameter of body causing tidal deformation
        std::vector< std::pair< int, int > > variationAndBodyIndices;
        for( unsigned int i = 0; i < tidalGravityFieldVariations_.size( ); i++ )
        {
            std::vector< std::string > deformingBodies = tidalGravityFieldVariations_.at( i )->getDeformingBodies( );
            for( unsigned int j = 0; j < deformingBodies.size( ); j++ )
            {
                if( deformingBodies.at( j ) == parameter->getParameterName( ).second.first )
                {
                    variationAndBodyIndices.push_back( std::make_pair( i, j ) );
                }
            }
        }

        if( variationAndBodyIndices.size( ) > 0 )
        {
            partialFunction = std::make_pair( std::bind( &InertialTorquePartial::wrtGravitationalParameterOfDeformingBody,
                                                         this,
                                                         std::placeholders::_1,
                                                         variationAndBodyIndices ),
                                              1 );
        }
    }
    return partialFunction;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case spherical_harmonics_cosine_coefficient_block: {
                // Cast parameter object to required type.
                std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );

                int c20Index, c21Index, c22Index;
                coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );

                if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
                {
                    if( getInertiaTensorNormalizationFactor_ == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when getting partial of inertial torque w.r.t. cosine sh parameters, inertia tensor normalization "
                                "function not found." );
                    }
                    partialFunction =
                            std::make_pair( std::bind( &InertialTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody,
                                                       this,
                                                       std::placeholders::_1,
                                                       c20Index,
                                                       c21Index,
                                                       c22Index ),
                                            coefficientsParameter->getParameterSize( ) );
                }

                break;
            }
            case spherical_harmonics_sine_coefficient_block: {
                // Cast parameter object to required type.

                std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );

                int s21Index, s22Index;
                coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );

                if( s21Index >= 0 || s22Index >= 0 )
                {
                    if( ( getInertiaTensorNormalizationFactor_ == nullptr ) )
                    {
                        throw std::runtime_error(
                                "Error when getting partial of inertial torque w.r.t. sine sh parameters, inertia tensor normalization "
                                "function not found." );
                    }
                    partialFunction = std::make_pair( std::bind( &InertialTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 s21Index,
                                                                 s22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }

                break;
            }
            case full_degree_tidal_love_number:
            case single_degree_variable_tidal_love_number:
            case mode_coupled_tidal_love_numbers:
                partialFunction = getTidalLoveNumberPartialFunction( parameter );
                break;
            default:
                break;
        }
    }
    return partialFunction;
}

//! Function to compute partial of torque w.r.t. mean moment of inertia
void InertialTorquePartial::wrtMeanMomentOfInertia( Eigen::MatrixXd& momentOfInertiaPartial )
{
    momentOfInertiaPartial.block( 0, 0, 3, 1 ) = -currentAngularVelocityCrossProductMatrix_ * currentInertiaTensorNormalizationFactor_ *
            UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_MEAN_MOMENT * currentAngularVelocityVector_;
}

//! Function to compute partial of torque w.r.t. gravitational parameter
void InertialTorquePartial::wrtGravitationalParameter( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = -currentAngularVelocityCrossProductMatrix_ * currentInertiaTensor_ *
            currentAngularVelocityVector_ / currentGravitationalParameter_;

    // The contribution of the tidal forcing to the inertia tensor is independent of the gravitational parameter (mass ratio is
    // inversely proportional to it), so it is removed from the above
    for( unsigned int i = 0; i < tidalGravityFieldVariations_.size( ); i++ )
    {
        for( unsigned int j = 0; j < tidalGravityFieldVariations_.at( i )->getDeformingBodies( ).size( ); j++ )
        {
            std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficientPartials =
                    tidalGravityFieldVariations_.at( i )
                            ->calculateSphericalHarmonicsCorrectionsPartialWrtDeformingBodyGravitationalParameter( currentTidalPartialTime_,
                                                                                                                   j );
            gravitationalParameterPartial.block( 0, 0, 3, 1 ) -=
                    tidalGravityFieldVariations_.at( i )->getDeformingBodyMasses( ).at( j )( ) *
                    getTorquePartialFromDegreeTwoCoefficients( coefficientPartials.first, coefficientPartials.second ) /
                    currentGravitationalParameter_;
        }
    }

    // Add contribution of inertia tensor derivative. The tidal contribution is independent of the gravitational parameter, while
    // the other contributions are proportional to it.
    if( inertiaTensorDerivativeFunction_ != nullptr )
    {
        Eigen::Vector3d nonTidalInertiaTensorDerivativeTorque = -inertiaTensorDerivativeFunction_( ) * currentAngularVelocityVector_;
        for( unsigned int i = 0; i < tidalGravityFieldVariations_.size( ); i++ )
        {
            std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficientRates =
                    tidalGravityFieldVariations_.at( i )->calculateSphericalHarmonicsCorrectionsTimeDerivative( currentTidalPartialTime_ );
            nonTidalInertiaTensorDerivativeTorque -=
                    getTorquePartialFromDegreeTwoCoefficientRates( coefficientRates.first, coefficientRates.second );
        }
        gravitationalParameterPartial.block( 0, 0, 3, 1 ) += nonTidalInertiaTensorDerivativeTorque / currentGravitationalParameter_;
    }
}

//! Function to compute the Hessian of the potential of a single spherical harmonic tidal forcing term
Eigen::Matrix3d InertialTorquePartial::getTidalForcingHessian( const double referenceRadius,
                                                               const int degree,
                                                               const int order,
                                                               const bool isCosineTerm,
                                                               const Eigen::Vector3d& bodyFixedPosition )
{
    // Create new cache (with second derivatives) if current cache does not support degree
    if( tidalForcingHessianCache_.getMaximumDegree( ) < degree + 2 )
    {
        tidalForcingHessianCache_ = basic_mathematics::SphericalHarmonicsCache( degree + 2, degree + 2 );
        tidalForcingHessianCache_.getLegendreCache( ).setComputeSecondDerivatives( 1 );
    }

    Eigen::MatrixXd unitCoefficients = Eigen::MatrixXd::Zero( degree + 1, degree + 1 );
    unitCoefficients( degree, order ) = 1.0;
    const Eigen::MatrixXd zeroCoefficients = Eigen::MatrixXd::Zero( degree + 1, degree + 1 );
    return computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration( bodyFixedPosition,
                                                                             referenceRadius,
                                                                             referenceRadius,
                                                                             isCosineTerm ? unitCoefficients : zeroCoefficients,
                                                                             isCosineTerm ? zeroCoefficients : unitCoefficients,
                                                                             tidalForcingHessianCache_ );
}

//! Function to compute partial of the inertia tensor w.r.t. a (geodesy-normalized) degree two cosine or sine coefficient
Eigen::Matrix3d InertialTorquePartial::getInertiaTensorPartialWrtDegreeTwoCoefficient( const int order, const bool isCosineCoefficient )
{
    Eigen::Matrix3d unscaledInertiaTensorPartial;
    if( isCosineCoefficient )
    {
        unscaledInertiaTensorPartial = ( order == 0 )
                ? UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20
                : ( ( order == 1 ) ? UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 : UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 );
    }
    else
    {
        unscaledInertiaTensorPartial = ( order == 0 )
                ? Eigen::Matrix3d::Zero( )
                : ( ( order == 1 ) ? UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 : UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 );
    }
    return basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, order ) * currentInertiaTensorNormalizationFactor_ *
            unscaledInertiaTensorPartial;
}

//! Function to compute partial of the inertia tensor derivative contribution to the torque w.r.t. a degree two coefficient rate
Eigen::Vector3d InertialTorquePartial::getTorquePartialWrtDegreeTwoCoefficientRate( const int order, const bool isCosineCoefficient )
{
    return -getInertiaTensorPartialWrtDegreeTwoCoefficient( order, isCosineCoefficient ) * currentAngularVelocityVector_;
}

//! Function to compute the partial of the torque due to (partials of) degree two coefficients, through the inertia tensor
Eigen::Vector3d InertialTorquePartial::getTorquePartialFromDegreeTwoCoefficients( const Eigen::MatrixXd& cosineCoefficientPartials,
                                                                                  const Eigen::MatrixXd& sineCoefficientPartials )
{
    Eigen::Matrix3d inertiaTensorPartial = Eigen::Matrix3d::Zero( );
    for( int order = 0; order < std::min( 3, static_cast< int >( cosineCoefficientPartials.cols( ) ) ); order++ )
    {
        inertiaTensorPartial += getInertiaTensorPartialWrtDegreeTwoCoefficient( order, true ) * cosineCoefficientPartials( 0, order ) +
                getInertiaTensorPartialWrtDegreeTwoCoefficient( order, false ) * sineCoefficientPartials( 0, order );
    }
    return -currentAngularVelocityCrossProductMatrix_ * inertiaTensorPartial * currentAngularVelocityVector_;
}

//! Function to compute the inertia tensor derivative contribution to the torque (or its partial) from degree two coefficient rates
Eigen::Vector3d InertialTorquePartial::getTorquePartialFromDegreeTwoCoefficientRates( const Eigen::MatrixXd& cosineCoefficientRates,
                                                                                      const Eigen::MatrixXd& sineCoefficientRates )
{
    Eigen::Vector3d torquePartial = Eigen::Vector3d::Zero( );
    for( int order = 0; order < std::min( 3, static_cast< int >( cosineCoefficientRates.cols( ) ) ); order++ )
    {
        torquePartial += getTorquePartialWrtDegreeTwoCoefficientRate( order, true ) * cosineCoefficientRates( 0, order ) +
                getTorquePartialWrtDegreeTwoCoefficientRate( order, false ) * sineCoefficientRates( 0, order );
    }
    return torquePartial;
}

//! Function to create a function computing partial of torque w.r.t. tidal Love numbers (through the inertia tensor derivative)
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > InertialTorquePartial::getTidalLoveNumberPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    // Retrieve partials of degree two coefficients and their rates w.r.t. parameter from all relevant tidal variations
    std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > > coefficientPartialFunctions;
    std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > > coefficientRatePartialFunctions;
    for( int computeRates = 0; computeRates < 2; computeRates++ )
    {
        const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >& loveNumberInterfaces =
                computeRates ? tidalLoveNumberRatePartialInterfaces_ : tidalLoveNumberPartialInterfaces_;
        for( unsigned int i = 0; i < loveNumberInterfaces.size( ); i++ )
        {
            std::pair< int, std::pair< int, int > > currentTidalPartialOutput =
                    loveNumberInterfaces.at( i )->setParameterPartialFunction( parameter, 2, 2 );
            if( currentTidalPartialOutput.first > 0 )
            {
                ( computeRates ? coefficientRatePartialFunctions : coefficientPartialFunctions )
                        .push_back( std::bind( &orbit_determination::TidalLoveNumberPartialInterface::getCurrentVectorParameterPartial,
                                               loveNumberInterfaces.at( i ),
                                               parameter,
                                               currentTidalPartialOutput.second ) );
            }
        }
    }

    if( coefficientPartialFunctions.size( ) > 0 || coefficientRatePartialFunctions.size( ) > 0 )
    {
        if( parameter->getParameterName( ).first == estimatable_parameters::mode_coupled_tidal_love_numbers )
        {
            std::shared_ptr< estimatable_parameters::ModeCoupledTidalLoveNumber > tidalLoveNumber =
                    std::dynamic_pointer_cast< estimatable_parameters::ModeCoupledTidalLoveNumber >( parameter );
            partialFunction.first = std::bind( &InertialTorquePartial::wrtModeCoupledTidalLoveNumbers,
                                               this,
                                               std::placeholders::_1,
                                               coefficientPartialFunctions,
                                               coefficientRatePartialFunctions,
                                               tidalLoveNumber->getResponseIndices( ),
                                               tidalLoveNumber->getResponseDegreeOrders( ),
                                               parameter->getParameterSize( ) );
        }
        else
        {
            std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > tidalLoveNumber =
                    std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > >( parameter );
            partialFunction.first = std::bind( &InertialTorquePartial::wrtTidalLoveNumbers,
                                               this,
                                               std::placeholders::_1,
                                               coefficientPartialFunctions,
                                               coefficientRatePartialFunctions,
                                               tidalLoveNumber->getOrders( ),
                                               tidalLoveNumber->getSumOrders( ),
                                               parameter->getParameterSize( ) );
        }
        partialFunction.second = parameter->getParameterSize( );
    }
    return partialFunction;
}

//! Function to compute partial of torque w.r.t. (non-mode coupled) tidal Love numbers
void InertialTorquePartial::wrtTidalLoveNumbers(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientPartialFunctions,
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientRatePartialFunctions,
        const std::vector< int >& orders,
        const bool sumOrders,
        const int parameterSize )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, parameterSize );

    // Add contributions through inertia tensor (computeRates = 0) and its time derivative (computeRates = 1)
    for( int computeRates = 0; computeRates < 2; computeRates++ )
    {
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& partialFunctions =
                computeRates ? coefficientRatePartialFunctions : coefficientPartialFunctions;
        for( unsigned int i = 0; i < partialFunctions.size( ); i++ )
        {
            std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > partialsPerOrder = partialFunctions.at( i )( );
            for( unsigned int j = 0; j < orders.size( ); j++ )
            {
                const int singleOrderPartialSize = partialsPerOrder.at( j ).cols( );
                Eigen::MatrixXd cosineCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
                Eigen::MatrixXd sineCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
                for( int k = 0; k < singleOrderPartialSize; k++ )
                {
                    cosineCoefficientPartials( 0, orders.at( j ) ) = partialsPerOrder.at( j )( 0, k );
                    sineCoefficientPartials( 0, orders.at( j ) ) = partialsPerOrder.at( j )( 1, k );
                    partialMatrix.block( 0, ( sumOrders ? 0 : j * singleOrderPartialSize ) + k, 3, 1 ) += computeRates
                            ? getTorquePartialFromDegreeTwoCoefficientRates( cosineCoefficientPartials, sineCoefficientPartials )
                            : getTorquePartialFromDegreeTwoCoefficients( cosineCoefficientPartials, sineCoefficientPartials );
                }
            }
        }
    }
}

//! Function to compute partial of torque w.r.t. mode-coupled tidal Love numbers
void InertialTorquePartial::wrtModeCoupledTidalLoveNumbers(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientPartialFunctions,
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& coefficientRatePartialFunctions,
        const std::vector< int >& responseIndices,
        const std::vector< std::pair< int, int > >& responseDegreeOrders,
        const int parameterSize )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, parameterSize );

    // Add contributions through inertia tensor (computeRates = 0) and its time derivative (computeRates = 1)
    for( int computeRates = 0; computeRates < 2; computeRates++ )
    {
        const std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > >& partialFunctions =
                computeRates ? coefficientRatePartialFunctions : coefficientPartialFunctions;
        for( unsigned int i = 0; i < partialFunctions.size( ); i++ )
        {
            std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > forcingPartials = partialFunctions.at( i )( );
            for( unsigned int j = 0; j < responseIndices.size( ); j++ )
            {
                // Only degree two responses contribute to the inertia tensor
                const std::pair< int, int > responseDegreeOrder = responseDegreeOrders.at( responseIndices.at( j ) );
                if( responseDegreeOrder.first == 2 )
                {
                    Eigen::MatrixXd cosineCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
                    Eigen::MatrixXd sineCoefficientPartials = Eigen::MatrixXd::Zero( 1, 3 );
                    cosineCoefficientPartials( 0, responseDegreeOrder.second ) = forcingPartials.at( j )( 0 );
                    sineCoefficientPartials( 0, responseDegreeOrder.second ) = forcingPartials.at( j )( 1 );
                    partialMatrix.block( 0, j, 3, 1 ) += computeRates
                            ? getTorquePartialFromDegreeTwoCoefficientRates( cosineCoefficientPartials, sineCoefficientPartials )
                            : getTorquePartialFromDegreeTwoCoefficients( cosineCoefficientPartials, sineCoefficientPartials );
                }
            }
        }
    }
}

//! Function to compute partial of torque w.r.t. gravitational parameter of a body causing tidal deformation
void InertialTorquePartial::wrtGravitationalParameterOfDeformingBody( Eigen::MatrixXd& gravitationalParameterPartial,
                                                                      const std::vector< std::pair< int, int > >& variationAndBodyIndices )
{
    // Tidal forcing due to a single body is proportional to its gravitational parameter
    gravitationalParameterPartial.block( 0, 0, 3, 1 ).setZero( );
    for( unsigned int i = 0; i < variationAndBodyIndices.size( ); i++ )
    {
        std::shared_ptr< gravitation::SolidBodyTideGravityFieldVariations > tidalVariation =
                tidalGravityFieldVariations_.at( variationAndBodyIndices.at( i ).first );
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficientRates =
                tidalVariation->calculateSphericalHarmonicsCorrectionsTimeDerivativeDueToBody(
                        currentTidalPartialTime_, variationAndBodyIndices.at( i ).second, 2 );
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficientPartials =
                tidalVariation->calculateSphericalHarmonicsCorrectionsPartialWrtDeformingBodyGravitationalParameter(
                        currentTidalPartialTime_, variationAndBodyIndices.at( i ).second );
        gravitationalParameterPartial.block( 0, 0, 3, 1 ) +=
                getTorquePartialFromDegreeTwoCoefficientRates( coefficientRates.first, coefficientRates.second ) /
                        tidalVariation->getDeformingBodyMasses( ).at( variationAndBodyIndices.at( i ).second )( ) +
                getTorquePartialFromDegreeTwoCoefficients( coefficientPartials.first, coefficientPartials.second );
    }
}

//! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
void InertialTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                                 const int c20Index,
                                                                                 const int c21Index,
                                                                                 const int c22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( c20Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c20Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 * currentAngularVelocityVector_;
    }

    if( c21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 * currentAngularVelocityVector_;
    }

    if( c22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 * currentAngularVelocityVector_;
    }
}

//! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
void InertialTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                               const int s21Index,
                                                                               const int s22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( s21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 * currentAngularVelocityVector_;
    }

    if( s22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentAngularVelocityCrossProductMatrix_ *
                currentInertiaTensorNormalizationFactor_ * UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 * currentAngularVelocityVector_;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
