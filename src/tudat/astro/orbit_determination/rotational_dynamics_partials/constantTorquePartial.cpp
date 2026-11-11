/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/constantTorquePartial.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to add the dependency of the tidal part of the inertia tensor on a parameter to a partial function
template< typename ParameterType >
void ConstantTorquePartial::addTidalInertiaTensorPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameter,
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& partialFunction )
{
    if( inertialTorquePartial_ != nullptr )
    {
        std::pair< InertialTorquePartial::TidalInertiaTensorPartialFunction, int > tidalPartialFunction =
                inertialTorquePartial_->getTidalInertiaTensorPartialFunction( parameter );
        if( tidalPartialFunction.second > 0 )
        {
            partialFunction = std::make_pair( std::bind( &ConstantTorquePartial::wrtParameterThroughTidalInertiaTensor,
                                                         this,
                                                         std::placeholders::_1,
                                                         partialFunction.first,
                                                         tidalPartialFunction.first ),
                                              tidalPartialFunction.second );
        }
    }
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > ConstantTorquePartial::getParameterPartialFunction(
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
                // Function is only defined if the inertia tensor is derived from the gravity field
                if( bodyGravitationalParameterFunction_ == nullptr )
                {
                    break;
                }

                partialFunction =
                        std::make_pair( std::bind( &ConstantTorquePartial::wrtGravitationalParameter, this, std::placeholders::_1 ), 1 );
                break;
            }
            case mean_moment_of_inertia: {
                if( ( getInertiaTensorNormalizationFactor_ == nullptr ) )
                {
                    throw std::runtime_error(
                            "Error when getting partial of inertial torque w.r.t. mean moment of inertia, inertia tensor normalization "
                            "function not found." );
                }
                partialFunction =
                        std::make_pair( std::bind( &ConstantTorquePartial::wrtMeanMomentOfInertia, this, std::placeholders::_1 ), 1 );
                break;
            }
            default:
                break;
        }
    }
    addTidalInertiaTensorPartialFunction( parameter, partialFunction );
    return partialFunction;
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > ConstantTorquePartial::getParameterPartialFunction(
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
                    // Function is only defined if the inertia tensor is derived from the gravity field
                    if( getInertiaTensorNormalizationFactor_ == nullptr )
                    {
                        break;
                    }
                    partialFunction =
                            std::make_pair( std::bind( &ConstantTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody,
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
                    // Function is only defined if the inertia tensor is derived from the gravity field
                    if( getInertiaTensorNormalizationFactor_ == nullptr )
                    {
                        break;
                    }
                    partialFunction = std::make_pair( std::bind( &ConstantTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 s21Index,
                                                                 s22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }

                break;
            }
            default:
                break;
        }
    }
    addTidalInertiaTensorPartialFunction( parameter, partialFunction );
    return partialFunction;
}

//! Function to compute partial of torque w.r.t. a parameter, including the dependency of the tidal part of the inertia tensor
void ConstantTorquePartial::wrtParameterThroughTidalInertiaTensor(
        Eigen::MatrixXd& partialMatrix,
        const std::function< void( Eigen::MatrixXd& ) > nonTidalPartialFunction,
        const InertialTorquePartial::TidalInertiaTensorPartialFunction tidalPartialFunction )
{
    std::vector< Eigen::Matrix3d > inertiaTensorPartials, inertiaTensorDerivativePartials;
    tidalPartialFunction( inertiaTensorPartials, inertiaTensorDerivativePartials );
    if( nonTidalPartialFunction != nullptr )
    {
        nonTidalPartialFunction( partialMatrix );
    }
    else
    {
        partialMatrix.setZero( 3, inertiaTensorPartials.size( ) );
    }

    for( unsigned int i = 0; i < inertiaTensorPartials.size( ); i++ )
    {
        partialMatrix.block( 0, i, 3, 1 ) -= inertiaTensorPartials.at( i ) * currentInverseInertiaTensor_ * currentTotalTorque_;
    }
}

//! Function to compute partial of torque w.r.t. mean moment of inertia
void ConstantTorquePartial::wrtMeanMomentOfInertia( Eigen::MatrixXd& momentOfInertiaPartial )
{
    momentOfInertiaPartial.block( 0, 0, 3, 1 ) = -currentInertiaTensorNormalizationFactor_ *
            UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_MEAN_MOMENT * currentInverseInertiaTensor_ * currentTotalTorque_;
}

//! Function to compute partial of torque w.r.t. gravitational parameter
void ConstantTorquePartial::wrtGravitationalParameter( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) =
            -currentInertiaTensor_ / bodyGravitationalParameterFunction_( ) * currentInverseInertiaTensor_ * currentTotalTorque_;
}

//! Function to compute partial of torque w.r.t. spherical harmonic cosine coefficients
void ConstantTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                                 const int c20Index,
                                                                                 const int c21Index,
                                                                                 const int c22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( c20Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c20Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) * currentInertiaTensorNormalizationFactor_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20 * currentInverseInertiaTensor_ * currentTotalTorque_;
    }

    if( c21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentInertiaTensorNormalizationFactor_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21 * currentInverseInertiaTensor_ * currentTotalTorque_;
    }

    if( c22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, c22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentInertiaTensorNormalizationFactor_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22 * currentInverseInertiaTensor_ * currentTotalTorque_;
    }
}

//! Function to compute partial of torque w.r.t. spherical harmonic sine coefficients
void ConstantTorquePartial::wrtSineSphericalHarmonicCoefficientsOfCentralBody( Eigen::MatrixXd& sphericalHarmonicCoefficientPartial,
                                                                               const int s21Index,
                                                                               const int s22Index )
{
    sphericalHarmonicCoefficientPartial.setZero( );

    if( s21Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s21Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 ) * currentInertiaTensorNormalizationFactor_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21 * currentInverseInertiaTensor_ * currentTotalTorque_;
    }

    if( s22Index >= 0 )
    {
        sphericalHarmonicCoefficientPartial.block( 0, s22Index, 3, 1 ) =
                -basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ) * currentInertiaTensorNormalizationFactor_ *
                UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22 * currentInverseInertiaTensor_ * currentTotalTorque_;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
