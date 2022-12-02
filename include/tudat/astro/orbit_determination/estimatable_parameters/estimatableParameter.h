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

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <Eigen/Geometry>

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
    arc_wise_time_drift_observation_bias
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

//! Function to determine whether the given parameter influences a body's tidal gravity field variations.
/*!
 * Function to determine whether the given parameter influences a body's tidal gravity field variations.
 * \param parameterType Parameter identifier.
 * \return True if parameter influences a body's tidal gravity field variations.
 */
bool isParameterTidalProperty( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given parameter represents an arc-wise initial dynamical state.
/*!
 * Function to determine whether the given parameter represents an arc-wise initial dynamical state.
 * \param parameterType Parameter identifier.
 * \return True if parameter is an arc-wise initial dynamical state.
 */
bool isParameterArcWiseInitialStateProperty( const EstimatebleParametersEnum parameterType );

//! Typedef for full parameter identifier.
typedef std::pair< EstimatebleParametersEnum, std::pair< std::string, std::string > > EstimatebleParameterIdentifier;



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



protected:

    //! Identifier of parameter.
    EstimatebleParameterIdentifier parameterName_;
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
