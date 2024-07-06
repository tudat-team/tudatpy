/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATABLEPARAMETERS_H
#define TUDAT_ESTIMATABLEPARAMETERS_H

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <map>


#include <type_traits>
#include <memory>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace estimatable_parameters
{

//! List of parameters that can be estimated by the orbit determination code.
enum EstimatebleParametersEnum
{
    arc_wise_initial_body_state,
    initial_body_state,
    initial_rotational_body_state,
    initial_mass_state,
    gravitational_parameter,
    constant_drag_coefficient,
    radiation_pressure_coefficient,
    arc_wise_radiation_pressure_coefficient,
    spherical_harmonics_cosine_coefficient_block,
    spherical_harmonics_sine_coefficient_block,
    constant_rotation_rate,
    rotation_pole_position,
    constant_additive_observation_bias,
    arcwise_constant_additive_observation_bias,
    constant_relative_observation_bias,
    arcwise_constant_relative_observation_bias,
    ppn_parameter_gamma,
    ppn_parameter_beta,
    ground_station_position,
    equivalence_principle_lpi_violation_parameter,
    empirical_acceleration_coefficients,
    arc_wise_empirical_acceleration_coefficients,
    full_degree_tidal_love_number,
    single_degree_variable_tidal_love_number,
    direct_dissipation_tidal_time_lag,
    mean_moment_of_inertia,
    arc_wise_constant_drag_coefficient,
    periodic_spin_variation,
    polar_motion_amplitude,
    core_factor,
    free_core_nutation_rate,
    desaturation_delta_v_values,
    scaled_longitude_libration_amplitude,
    constant_thrust_magnitude_parameter,
    constant_specific_impulse,
    constant_time_drift_observation_bias,
    arc_wise_time_drift_observation_bias,
    constant_time_observation_bias,
    arc_wise_time_observation_bias,
    inverse_tidal_quality_factor,
    yarkovsky_parameter,
    custom_estimated_parameter,
    reference_point_position,
    polynomial_gravity_field_variation_amplitudes,
    periodic_gravity_field_variation_amplitudes
};

std::string getParameterTypeString( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
/*!
 * Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
 * \param parameterType Parameter identifier.
 * \return True if parameter is an initial dynamical state.
 */
bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given (non-dynamical) parameter is a double or vector parameter.
/*!
 * Function to determine whether the given (non-dynamical) parameter is a double or vector parameter.
 * \param parameterType Parameter identifier.
 * \return True if parameter is a double parameter.
 */
bool isDoubleParameter( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given (non-dynamical) parameter influences a body's orientation.
/*!
 * Function to determine whether the given (non-dynamical) parameter influences a body's orientation.
 * \param parameterType Parameter identifier.
 * \return True if parameter is a property of rotation model
 */
bool isParameterRotationMatrixProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter influences an observation link directly
/*!
 * Function to determine whether the given parameter influences an observation link directly, such as observation biases or
 * clock parameters
 * \param parameterType Parameter identifier.
 * \return True if parameter is a property of an observation link
 */
bool isParameterObservationLinkProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter influences an observation time directly
/*!
 * Function to determine whether the given parameter influences an observation time directly, such as clock parameters
 * \param parameterType Parameter identifier.
 * \return True if parameter is a time property of an observation link
 */
bool isParameterObservationLinkTimeProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter influences a body's tidal gravity field variations.
/*!
 * Function to determine whether the given parameter influences a body's tidal gravity field variations.
 * \param parameterType Parameter identifier.
 * \return True if parameter influences a body's tidal gravity field variations.
 */
bool isParameterTidalProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter influences a body's non-tidal gravity field variations.
bool isParameterNonTidalGravityFieldVariationProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter represents an arc-wise initial dynamical state.
/*!
 * Function to determine whether the given parameter represents an arc-wise initial dynamical state.
 * \param parameterType Parameter identifier.
 * \return True if parameter is an arc-wise initial dynamical state.
 */
bool isParameterArcWiseInitialStateProperty( const EstimatebleParametersEnum parameterType );

//! Typedef for full parameter identifier.
typedef std::pair< EstimatebleParametersEnum, std::pair< std::string, std::string > > EstimatebleParameterIdentifier;


class CustomAccelerationPartialSettings
{
public:

    CustomAccelerationPartialSettings(
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const basic_astrodynamics::AvailableAcceleration accelerationType ):
        bodyUndergoingAcceleration_( bodyUndergoingAcceleration ),
        bodyExertingAcceleration_( bodyExertingAcceleration ),
        accelerationType_( accelerationType )
        { }

    virtual ~CustomAccelerationPartialSettings( ){ }

    std::string bodyUndergoingAcceleration_;

    std::string bodyExertingAcceleration_;

    basic_astrodynamics::AvailableAcceleration accelerationType_;

    bool accelerationMatches(
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const basic_astrodynamics::AvailableAcceleration accelerationType )
    {
        return ( accelerationType_ == accelerationType ) &&
        ( bodyUndergoingAcceleration_ == bodyUndergoingAcceleration ) &&
        ( bodyExertingAcceleration_ == bodyExertingAcceleration );
    }
};

class NumericalAccelerationPartialSettings: public CustomAccelerationPartialSettings
{
public:

    NumericalAccelerationPartialSettings(
        const Eigen::VectorXd& parameterPerturbation,
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& environmentUpdateSettings =
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >( ) ):
        CustomAccelerationPartialSettings( bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationType ),
        parameterPerturbation_( parameterPerturbation ), environmentUpdateSettings_( environmentUpdateSettings ){ }

    Eigen::VectorXd parameterPerturbation_;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentUpdateSettings_;
};

class AnalyticalAccelerationPartialSettings: public CustomAccelerationPartialSettings
{
public:

    AnalyticalAccelerationPartialSettings(
        const std::function< Eigen::MatrixXd( const double, const Eigen::Vector3d& ) > accelerationPartialFunction,
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const basic_astrodynamics::AvailableAcceleration accelerationType ):
        CustomAccelerationPartialSettings( bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationType ),
    accelerationPartialFunction_( accelerationPartialFunction ){ }

    std::function< Eigen::MatrixXd( const double, const Eigen::Vector3d& ) > accelerationPartialFunction_;
};


inline std::shared_ptr< CustomAccelerationPartialSettings > analyticalAccelerationPartialSettings(
    const std::function< Eigen::MatrixXd( const double, const Eigen::Vector3d& ) > accelerationPartialFunction,
    const std::string bodyUndergoingAcceleration,
    const std::string bodyExertingAcceleration,
    const basic_astrodynamics::AvailableAcceleration accelerationType )
{
    return std::make_shared< AnalyticalAccelerationPartialSettings >(
        accelerationPartialFunction, bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationType );
}

inline std::shared_ptr< CustomAccelerationPartialSettings > numericalAccelerationPartialSettings(
    const Eigen::VectorXd& parameterPerturbation,
    const std::string bodyUndergoingAcceleration,
    const std::string bodyExertingAcceleration,
    const basic_astrodynamics::AvailableAcceleration accelerationType,
    const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& environmentUpdateSettings  =
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >( ) )
{
    return std::make_shared< NumericalAccelerationPartialSettings >(
        parameterPerturbation, bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationType, environmentUpdateSettings );
}

//! Base class for a parameter that is to be estimated.
/*!
 *  Base class for a parameter that is to be estimated. A separate derived class is to be made for each type of parameter
 *  (i.e. gravitational parameter, initial translational state, etc. ).
 */
template< typename ParameterType >
class EstimatableParameter
{

public:
    //! Constructor.
    /*!
     *  Constructor taking parameter name and associated body. All parameters are identified by a these two variables.
     *  Any additional information that may be required for uniquely defining a parameter is to be defined in the derived class.
     *  \param parameterName Enum value defining the type of the parameter.
     *  \param associatedBody Name of body associated with patameters
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    EstimatableParameter( const EstimatebleParametersEnum parameterName,
                          const std::string& associatedBody,
                          const std::string& pointOnBodyId = ""  ):
        parameterName_( std::make_pair( parameterName, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor.
    virtual ~EstimatableParameter( ) { }

    //! Pure virtual function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    virtual ParameterType getParameterValue( ) = 0;

    //! Pure virtual function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param parameterValue to which the parameter is to be set.
     */
    virtual void setParameterValue( const ParameterType parameterValue ) = 0;

    //! Function to retrieve the type and associated body of the parameter.
    /*!
     *  Function to retrieve the type and associated body of the parameter.
     *  \return Identifier of parameter as a pair of parameter type and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier getParameterName( ) { return parameterName_; }

    virtual std::string getParameterDescription( )
    {
        std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "of (" + parameterName_.second.first;
        if( parameterName_.second.second == "" )
        {
            parameterDescription += ").";
        }
        else
        {
            parameterDescription += ", " + parameterName_.second.second + ").";
        }
        return parameterDescription;
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Pure virtual function to retrieve the size of the parameter (i.e. 1 for double parameters)
     *  \return Size of parameter value.
     */
    virtual int getParameterSize( ) = 0;

    //! Function to return additional identifier for parameter
    /*!
     *  Function to return additional identifier for parameter, beyond information stored in parameterName_, default
     *  none.
     *  \return Additional identifier for parameter (default empty string).
     */
    virtual std::string getSecondaryIdentifier( )
    {
        return "";
    }

    //! Function to retrieve size of constraint to be applied on parameter
    /*!
     * Function to retrieve size of constraint to be applied on parameter, zero by default. Can be overridden in derived class
     * \return Size of constraint to be applied on parameter
     */
    virtual int getConstraintSize( )
    {
        return 0;
    }

    //! Function to retrieve multiplier for parameter linear constraint
    /*!
     * Function to retrieve multiplier for parameter linear constraint, empty by default. Can be overridden in derived class
     * \return Multiplier for parameter linear constraint
     */
    virtual Eigen::MatrixXd getConstraintStateMultipler( )
    {
        return Eigen::MatrixXd::Zero( 0, 0 );
    }

    //! Function to retrieve right-hand side for parameter linear constraint
    /*!
     * Function to retrieve right-hand side for parameter linear constraint, empty by default. Can be overridden in derived class
     * \return Right-hand side for parameter linear constraint
     */
    virtual Eigen::VectorXd getConstraintRightHandSide( )
    {
        return Eigen::VectorXd::Zero( 0 );
    }

    virtual void throwExceptionIfNotFullyDefined( ){ }

    std::vector< std::shared_ptr< CustomAccelerationPartialSettings > > getCustomPartialSettings( )
    {
        return customPartialSettings_;
    }

    bool hasCustomAccelerationPartialSettings( const std::string bodyUndergoingAcceleration,
                                               const std::string bodyExertingAcceleration,
                                               const basic_astrodynamics::AvailableAcceleration accelerationType )
    {
        bool hasCustomSettings = false;
        for( unsigned int i = 0; i < customPartialSettings_.size( ); i++ )
        {
            if( customPartialSettings_.at( i )->accelerationMatches( bodyUndergoingAcceleration, bodyExertingAcceleration, accelerationType ) == true )
            {
                hasCustomSettings = true;
            }
        }
        return hasCustomSettings;
    }

    void addCustomPartialSettings( const std::shared_ptr< CustomAccelerationPartialSettings > customPartialSettings )
    {
        customPartialSettings_.push_back( customPartialSettings );
        checkCustomPartialSettings( );
    }

    void setCustomPartialSettings( const std::vector< std::shared_ptr< CustomAccelerationPartialSettings > > customPartialSettings )
    {
        customPartialSettings_ = customPartialSettings;
        checkCustomPartialSettings( );
    }

protected:

    //! Identifier of parameter.
    EstimatebleParameterIdentifier parameterName_;

    std::vector< std::shared_ptr< CustomAccelerationPartialSettings > > customPartialSettings_;

    void checkCustomPartialSettings( )
    {
        for( unsigned int i = 0; i < customPartialSettings_.size( ); i++ )
        {
            for( unsigned int j = i + 1; j < customPartialSettings_.size( ); j++ )
            {
                if( customPartialSettings_.at( i )->accelerationMatches(
                    customPartialSettings_.at( j )->bodyUndergoingAcceleration_,
                    customPartialSettings_.at( j )->bodyExertingAcceleration_,
                    customPartialSettings_.at( j )->accelerationType_  ) )
                {
                    throw std::runtime_error( "Error, found multiple custom partials for acceleration on (" +
                                                  customPartialSettings_.at( j )->bodyUndergoingAcceleration_ +") due to (" +
                                                  customPartialSettings_.at( j )->bodyExertingAcceleration_ + ") of type: " +
                                                  basic_astrodynamics::getAccelerationModelName( customPartialSettings_.at( j )->accelerationType_ ) );
                }
            }

        }
    }


};

class CustomEstimatableParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:
    //! Constructor.
    /*!
     *  Constructor taking parameter name and associated body. All parameters are identified by a these two variables.
     *  Any additional information that may be required for uniquely defining a parameter is to be defined in the derived class.
     *  \param parameterName Enum value defining the type of the parameter.
     *  \param associatedBody Name of body associated with patameters
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    CustomEstimatableParameter(
        const std::string& customId,
        const int parameterSize,
        const std::function< Eigen::VectorXd( ) > getParameterFunction,
        const std::function< void( const Eigen::VectorXd& ) > setParameterFunction ):
        EstimatableParameter< Eigen::VectorXd >( custom_estimated_parameter, "", customId ),
        parameterSize_( parameterSize ),
        getParameterFunction_( getParameterFunction ),
        setParameterFunction_( setParameterFunction ){ }

    //! Virtual destructor.
    ~CustomEstimatableParameter( ) { }

    //! Pure virtual function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd currentParameter = getParameterFunction_( );
        if( currentParameter.rows( ) != parameterSize_ )
        {
            throw std::runtime_error( "Error, when getting cutom parameter value, size is incompatible " +
                                          std::to_string( currentParameter.rows( ) ) + ", " +
                                          std::to_string( parameterSize_ ) );
        }
        return currentParameter;
    }

    //! Pure virtual function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param parameterValue to which the parameter is to be set.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        if( parameterValue.rows( ) != parameterSize_ )
        {
            throw std::runtime_error( "Error, when getting cutom parameter value, size is incompatible " +
                                      std::to_string( parameterValue.rows( ) ) + ", " +
                                      std::to_string( parameterSize_ ) );
        }
        setParameterFunction_( parameterValue );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Pure virtual function to retrieve the size of the parameter (i.e. 1 for double parameters)
     *  \return Size of parameter value.
     */
    int getParameterSize( )
    {
        return parameterSize_;
    }



protected:

    int parameterSize_;

    std::function< Eigen::VectorXd( ) > getParameterFunction_;

    std::function< void( const Eigen::VectorXd& ) > setParameterFunction_;

};

class CustomAccelerationPartialCalculator
{
public:

    CustomAccelerationPartialCalculator( ){ }

    virtual ~CustomAccelerationPartialCalculator( ){ }

    virtual Eigen::MatrixXd computePartial(
        const double currentTime, const Eigen::Vector3d& currentAcceleration,
        const std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel ) = 0;

protected:
//    std::string bodyUndergoingAcceleration_;
//
//    std::string bodyExertingAcceleration_;
//
//    basic_astrodynamics::AvailableAcceleration accelerationType_;
};


class NumericalAccelerationPartialWrtStateCalculator: public CustomAccelerationPartialCalculator
{
public:

    NumericalAccelerationPartialWrtStateCalculator(
        const Eigen::VectorXd& bodyStatePerturbations,
        const std::function< Eigen::Vector6d( ) > bodyStateGetFunction,
        const std::function< void( const Eigen::Vector6d& ) > bodyStateSetFunction,
        const std::function< void( const double ) > environmentUpdateFunction ):
        bodyStatePerturbations_( bodyStatePerturbations ),
        bodyStateGetFunction_( bodyStateGetFunction ),
        bodyStateSetFunction_( bodyStateSetFunction ),
        environmentUpdateFunction_( environmentUpdateFunction )
    {
        if( bodyStatePerturbations_.rows( ) != 6 )
        {
            throw std::runtime_error( "Error when making numerical acceleration partial for initial state parameter, sizes are inconsistent: " +
                                      std::to_string( bodyStatePerturbations_.rows( ) ) + ", " + std::to_string( 6 ) );
        }
        currentPartial_ = Eigen::MatrixXd::Zero( 3, 6 );

    }

    Eigen::MatrixXd computePartial(
        const double currentTime, const Eigen::Vector3d& currentAcceleration,
        const std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel )
    {
        currentPartial_.setZero( );

        Eigen::Vector6d nominalState = bodyStateGetFunction_( );
        Eigen::Vector6d perturbedState;

        // Compute state partial by numerical difference
        Eigen::Vector3d upperturbedAcceleration, downperturbedAcceleration;
        for( unsigned int i = 0; i < 6; i++ )
        {
            // Perturb state upwards
            perturbedState = nominalState;
            perturbedState( i ) += bodyStatePerturbations_( i );

            // Update environment/acceleration to perturbed state.
            accelerationModel->resetCurrentTime( );
            bodyStateSetFunction_( perturbedState );
            environmentUpdateFunction_( currentTime );
            accelerationModel->updateMembers( currentTime );

            // Retrieve perturbed acceleration.
            upperturbedAcceleration = accelerationModel->getAcceleration( );

            // Perturb state downwards
            perturbedState = nominalState;
            perturbedState( i ) -= bodyStatePerturbations_( i );

            // Update environment/acceleration to perturbed state.
            accelerationModel->resetCurrentTime( );
            bodyStateSetFunction_( perturbedState );
            environmentUpdateFunction_( currentTime );
            accelerationModel->updateMembers( currentTime );

            // Retrieve perturbed acceleration.
            downperturbedAcceleration = accelerationModel->getAcceleration( );

            // Compute partial
            currentPartial_.block( 0, i, 3, 1 ) =
                ( upperturbedAcceleration - downperturbedAcceleration ) / ( 2.0 * bodyStatePerturbations_( i ) );
        }
        return currentPartial_;
    }

protected:

    Eigen::VectorXd bodyStatePerturbations_;

    const std::function< Eigen::Vector6d( ) > bodyStateGetFunction_;

    const std::function< void( const Eigen::Vector6d& ) > bodyStateSetFunction_;

    std::function< void( const double ) > environmentUpdateFunction_;

    Eigen::Matrix< double, 3, 6 > currentPartial_;
};

class AnalyticalAccelerationPartialCalculator: public CustomAccelerationPartialCalculator
{
public:

    AnalyticalAccelerationPartialCalculator(
        const std::function<Eigen::MatrixXd( const double, const Eigen::Vector3d & )> accelerationPartialFunction,
        const int parameterSize,
        const std::string parameterName ) :
        accelerationPartialFunction_( accelerationPartialFunction ), parameterSize_( parameterSize ),
        parameterName_( parameterName )
    {
    }

    Eigen::MatrixXd computePartial(
        const double currentTime, const Eigen::Vector3d &currentAcceleration,
        const std::shared_ptr<basic_astrodynamics::AccelerationModel3d> accelerationModel )
    {
        Eigen::MatrixXd currentAccelerationPartial = accelerationPartialFunction_( currentTime, currentAcceleration );
        if ( currentAccelerationPartial.cols( ) != parameterSize_)
        {
            throw std::runtime_error( "Error when making numerical acceleration partial for parameter " +
                                      parameterName_ + ", sizes are inconsistent: " +
                                      std::to_string( currentAccelerationPartial.cols( )) + ", " +
                                      std::to_string( parameterSize_));
        }
        return currentAccelerationPartial;
    }

protected:


    std::function<Eigen::MatrixXd( const double, const Eigen::Vector3d & )> accelerationPartialFunction_;

    int parameterSize_;

    std::string parameterName_;
};

class CustomSingleAccelerationPartialCalculatorSet
{
public:
    CustomSingleAccelerationPartialCalculatorSet( ){ }

    std::map< estimatable_parameters::EstimatebleParameterIdentifier,
        std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > > customInitialStatePartials_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier,
        std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > > customDoubleParameterPartials_;

    std::map< estimatable_parameters::EstimatebleParameterIdentifier,
        std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > > customVectorParameterPartials_;

    int getNumberOfCustomPartials( )
    {
        return customInitialStatePartials_.size( ) + customDoubleParameterPartials_.size( ) + customVectorParameterPartials_.size( );
    }
};


//! Function to determine if an initial state parameter is a single- or multi-arc parameter
/*!
 *  Function to determine if an initial state parameter is a single- or multi-arc parameter. Function throws an error, if
 *  input is not an initial state parameter
 *  \param parameterToCheck Parameter object for which the check is to be performed.
 *  \return True of parameter is single-arc, false if multi-arc
 */
template< typename ParameterType >
bool isDynamicalParameterSingleArc(
        const std::shared_ptr< EstimatableParameter< Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > parameterToCheck )
{
    bool flag = -1;
    switch( parameterToCheck->getParameterName( ).first )
    {
    case arc_wise_initial_body_state:
    {
        flag = false;
        break;
    }
    case initial_body_state:
    {
        flag = true;
        break;
    }
    case initial_rotational_body_state:
    {
        flag = true;
        break;
    }
    case initial_mass_state:
    {
        flag = true;
        break;
    }
    default:
        throw std::runtime_error( "Error when checking single/multi-arc dynamical parameter, parameter not identified" );
    }
    return flag;

}



} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_ESTIMATABLEPARAMETERS_H
