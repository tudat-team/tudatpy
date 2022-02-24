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
    scaled_longitude_libration_amplitude
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

//extern template class EstimatableParameter< double >;
//extern template class EstimatableParameter< Eigen::VectorXd >;

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template class EstimatableParameter< Eigen::Matrix< long double, Eigen::Dynamic, 1 > >;
//#endif

//! Container class for all parameters that are to be estimated.
/*!
 *  Container class for all parameters that are to be estimated. Class is templated with the scalar type used for the
 *  estimation of any initial dynamical states that may be included
 */
template< typename InitialStateParameterType = double >
class EstimatableParameterSet
{
public:

    //! Constructor of parameter set.
    /*!
     *  Constructor of parameter set.
     *  \param estimatedDoubleParameters List of double parameters that are estimated.
     *  \param estimatedVectorParameters List of vector parameters that are estimated.
     *  \param estimateInitialStateParameters List of initial dynamical states that are to be estimated.
     */
    EstimatableParameterSet(
            const std::vector< std::shared_ptr< EstimatableParameter< double > > >& estimatedDoubleParameters,
            const std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& estimatedVectorParameters,
            const std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >& estimateInitialStateParameters =
            ( std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix
              < InitialStateParameterType, Eigen::Dynamic, 1 > > > >( ) ) ):
        estimatedDoubleParameters_( estimatedDoubleParameters ), estimatedVectorParameters_( estimatedVectorParameters ),
        totalConstraintSize_( 0 )
    {
        // Initialize total number of parameters to 0.
        estimatedParameterSetSize_ = 0;
        initialDynamicalStateParameterSize_ = 0;
        initialDynamicalSingleArcStateParameterSize_ = 0;
        initialDynamicalMultiArcStateParameterSize_ = 0;

        // Iterate over all double parameters and add to parameter size.
        for( unsigned int i = 0; i < estimateInitialStateParameters.size( ); i++ )
        {
            if( isDynamicalParameterSingleArc( estimateInitialStateParameters[ i ] ) )
            {
                estimateSingleArcInitialStateParameters_.push_back( estimateInitialStateParameters[ i ] );
            }
            else
            {
                estimateMultiArcInitialStateParameters_.push_back( estimateInitialStateParameters[ i ] );
            }
        }

        estimateInitialStateParameters_ = estimateSingleArcInitialStateParameters_;
        estimateInitialStateParameters_.insert(
                    estimateInitialStateParameters_.end( ), estimateMultiArcInitialStateParameters_.begin( ),
                    estimateMultiArcInitialStateParameters_.end( ) );

        for( unsigned int i = 0; i < estimateSingleArcInitialStateParameters_.size( ); i++ )
        {
            initialStateParameters_[ estimatedParameterSetSize_ ] = estimateSingleArcInitialStateParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimateInitialStateParameters_[ i ]->getParameterSize( ) ) );
            totalConstraintSize_ += estimateInitialStateParameters_[ i ]->getConstraintSize( );

            initialDynamicalSingleArcStateParameterSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
            initialSingleArcStateParameters_[ estimatedParameterSetSize_ ] = estimateSingleArcInitialStateParameters_[ i ];

            initialDynamicalStateParameterSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
            estimatedParameterSetSize_ += estimateSingleArcInitialStateParameters_[ i ]->getParameterSize( );
        }


        for( unsigned int i = 0; i < estimateMultiArcInitialStateParameters_.size( ); i++ )
        {
            initialStateParameters_[ estimatedParameterSetSize_ ] = estimateMultiArcInitialStateParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( ) ) );

            initialDynamicalMultiArcStateParameterSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
            initialMultiArcStateParameters_[ estimatedParameterSetSize_ ] = estimateMultiArcInitialStateParameters_[ i ];

            initialDynamicalStateParameterSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
            estimatedParameterSetSize_ += estimateMultiArcInitialStateParameters_[ i ]->getParameterSize( );
        }


        // Iterate over all double parameters and add to parameter size and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            doubleParameters_[ estimatedParameterSetSize_ ] = estimatedDoubleParameters_[ i ];
            totalConstraintSize_ += estimatedDoubleParameters_[ i ]->getConstraintSize( );

            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_, 1 ) );
            estimatedParameterSetSize_++;
        }

        // Iterate over all vector parameter, add to total number of parameters and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            vectorParameters_[ estimatedParameterSetSize_ ] = estimatedVectorParameters_[ i ];
            totalConstraintSize_ += estimatedVectorParameters_[ i ]->getConstraintSize( );

            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimatedVectorParameters_[ i ]->getParameterSize( ) ) );
            estimatedParameterSetSize_ += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        totalParameterSetSize_ = estimatedParameterSetSize_;
    }

    //! Function to return the total number of parameter values (including consider parameters)
    /*!
     *  Function to return the total number of parameter values (including consider parameters)
     *  \return Size of parameter vector (including consider parameters)
     */
    int getParameterSetSize( )
    {
        return totalParameterSetSize_;
    }

    //! Function to return the total number of parameter values (excluding consider parameters).
    /*!
     *  Function to return the total number of parameter values (excluding consider parameters)
     *  \return Size of parameter vector (excluding consider parameters)
     */
    int getEstimatedParameterSetSize( )
    {
        return estimatedParameterSetSize_;
    }

    //! Function to return the total number of initial state values that are estimated.
    /*!
     *  Function to return the total number of initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalStateParameterSize( )
    {
        return initialDynamicalStateParameterSize_;
    }

    //! Function to return the total number of single-arc initial state values that are estimated.
    /*!
     *  Function to return the total number of single-arc initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalSingleArcStateParameterSize( )
    {
        return initialDynamicalSingleArcStateParameterSize_;
    }

    //! Function to return the total number of multi-arc initial state values that are estimated.
    /*!
     *  Function to return the total number of multi-arc initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalMultiArcStateParameterSize( )
    {
        return initialDynamicalMultiArcStateParameterSize_;
    }

    //! Function that returns a vector containing all current parameter values
    /*!
     *  Function that returns a vector containing all current parameter values. The total vector starts with the initial
     *  state parameters, followed by the double and vector parameters, respectively.
     *  Initial state, double and vector parameter values are concatenated in the order in which they are set in the
     *  estimateInitialStateParameters_, doubleParameters_ and vectorParameters_ members.
     *  \return Vector containing all parameter values
     */
    template< typename ParameterScalar >
    Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 > getFullParameterValues( )
    {
        Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >  parameterValues =
                Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >::Zero( totalParameterSetSize_ );

        int currentStartIndex = 0;

        // Retrieve initial state parameter values.
        for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ) =
                    estimateInitialStateParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
        }

        // Retrieve double parameter values.
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            parameterValues( currentStartIndex ) = static_cast< ParameterScalar >(
                        estimatedDoubleParameters_[ i ]->getParameterValue( ) );
            currentStartIndex++;
        }

        // Retrieve vector parameter values.
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ) =
                    estimatedVectorParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        return parameterValues;
    }

    //! Function to reset all parameter values.
    /*!
     *  Function to reset all parameter values.
     *  \param newParameterValues New parameter values. Order of values in vector must be same order as return vector of getFullParameterValues
     */
    template< typename ParameterScalar >
    void resetParameterValues( const Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >& newParameterValues )
    {
        // Check input consistency
        if( newParameterValues.rows( ) != totalParameterSetSize_ )
        {
            throw std::runtime_error( "Error when resetting parameters of parameter set, given vector has size " +
                                      std::to_string( newParameterValues.rows( ) ) +
                                      ", while internal size is " + std::to_string( totalParameterSetSize_ ) );
        }
        else
        {
            int currentStartIndex = 0;

            for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
            {
                estimateInitialStateParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ).
                            template cast< InitialStateParameterType >( ) );
                currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
            }

            // Set double parameter values.
            for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
            {
                estimatedDoubleParameters_[ i ]->setParameterValue( static_cast< double >( newParameterValues( currentStartIndex ) ) );
                currentStartIndex++;
            }

            // Set vector parameter values.

            for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
            {
                estimatedVectorParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ).
                            template cast< double >( ) );

                currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
            }
        }
    }

    //! Function to retrieve double parameter objects.
    /*!
     *  Function to retrieve double parameter objects.
     *  \return Map containing all double parameter objects, with map key start index of parameter in total vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< double > > > getDoubleParameters( )
    {
        return doubleParameters_;
    }

    //! Function to retrieve vector parameter objects.
    /*!
     *  Function to retrieve vector parameter objects.
     *  \return Map containing all vector parameter objects, with map key start index of parameter in total vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getVectorParameters( )
    {
        return vectorParameters_;
    }

    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialStateParameters( )
    {
        return initialStateParameters_;
    }

    //! Function to retrieve all single-arc initial state parameter objects.
    /*!
     *  Function to retrieve all single-arc initial state parameter objects.
     *  \return Map containing all single-arc initial state parameter objects, with map key start index of parameter in total
     *  vector.
     */
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialSingleArcStateParameters( )
    {
        return initialSingleArcStateParameters_;
    }

    //! Function to retrieve all multi-arc initial state parameter objects.
    /*!
     *  Function to retrieve all multi-arc initial state parameter objects.
     *  \return Map containing all multi-arc initial state parameter objects, with map key start index of parameter in total
     * vector.
     */
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getInitialMultiArcStateParameters( )
    {
        return initialMultiArcStateParameters_;
    }

    std::vector< std::shared_ptr< EstimatableParameter< double > > > getEstimatedDoubleParameters( )
    {
        return estimatedDoubleParameters_;
    }

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getEstimatedVectorParameters( )
    {
        return estimatedVectorParameters_;
    }

    //! Function to retrieve the start index and size of (a) parameters(s) with a given identifier
    /*!
     * Function to retrieve the start index and size of (a) parameters(s) with a given identifier
     * \param requiredParameterId Parameter identifier that is to be searched in full list of patameters
     * \return List of start indices and sizes of parameters corresponding to requiredParameterId
     */
    std::vector< std::pair< int, int > > getIndicesForParameterType(
            const EstimatebleParameterIdentifier requiredParameterId )
    {
        std::vector< std::pair< int, int > > typeIndices;

        for( auto parameterIterator : initialSingleArcStateParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : initialMultiArcStateParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : doubleParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        for( auto parameterIterator : vectorParameters_ )
        {
            if( parameterIterator.second->getParameterName( ) == requiredParameterId )
            {
                typeIndices.push_back( std::make_pair( parameterIterator.first, parameterIterator.second->getParameterSize( ) ) );
            }
        }

        return typeIndices;
    }

    //! Function to get list of initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedInitialStateParameters( )
    {
        return estimateInitialStateParameters_;
    }

    //! Function to get list of single-arc initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of single-arc initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedSingleArcInitialStateParameters( )
    {
        return estimateSingleArcInitialStateParameters_;
    }

    //! Function to get list of multi-arc initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of multi-arc initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedMultiArcInitialStateParameters( )
    {
        return estimateMultiArcInitialStateParameters_;
    }

    //! Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
    /*!
     *  Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
     *  \return List of start indices and sizes (map keys) of estimated parameters.
     */
    std::vector< std::pair< int, int > > getParametersIndices( )
    {
        return parameterIndices_;
    }

    //! Function to retrieve total multiplier and right-hand side for parameter estimation linear constraint
    /*!
     * Function to retrieve total multiplier and right-hand side for parameter estimation linear constraint
     * \param constraintStateMultiplier Multiplier for parameter linear constraint
     * \param constraintRightHandSide Right-hand side for parameter linear constraint
     */
    void getConstraints( Eigen::MatrixXd& constraintStateMultiplier, Eigen::VectorXd& constraintRightHandSide )
    {
        // Resize constraint elements
        constraintStateMultiplier.setZero( totalConstraintSize_, estimatedParameterSetSize_ );
        constraintRightHandSide.setZero( totalConstraintSize_, 1 );

        // Iterate over all state parameters
        int currentConstraintRow = 0;
        int currentConstraintSize = 0;
        for( auto parameterIterator = initialStateParameters_.begin( ); parameterIterator != initialStateParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }

        }

        // Iterate over all double parameters
        for( auto parameterIterator = doubleParameters_.begin( ); parameterIterator != doubleParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }
        }

        // Iterate over all vector parameters
        for( auto parameterIterator = vectorParameters_.begin( ); parameterIterator != vectorParameters_.end( );
             parameterIterator++ )
        {
            // Add constraint if of non-zero size
            currentConstraintSize = parameterIterator->second->getConstraintSize( );
            if( currentConstraintSize > 0 )
            {
                constraintStateMultiplier.block(
                            currentConstraintRow, parameterIterator->first, currentConstraintSize,
                            parameterIterator->second->getParameterSize( )  ) =
                        parameterIterator->second->getConstraintStateMultipler( );
                constraintRightHandSide.segment( currentConstraintRow, currentConstraintSize ) =
                        parameterIterator->second->getConstraintRightHandSide( );

                currentConstraintRow += currentConstraintSize;
            }
        }
    }

    //! Total size of linear constraint that is to be applied during estimation
    /*!
     * Total size of linear constraint that is to be applied during estimation
     * \return Size of linear constraint that is to be applied during estimation
     */
    int getConstraintSize( )
    {
        return totalConstraintSize_;
    }

protected:

    //! Total size of all initial dynamical states that are to be estimated.
    int initialDynamicalStateParameterSize_;

    //! Total size of all initial single-arc dynamical states that are to be estimated.
    int initialDynamicalSingleArcStateParameterSize_;

    //! Total size of all initial multi-arc dynamical states that are to be estimated.
    int initialDynamicalMultiArcStateParameterSize_;

    //! Total number of parameter values (including currently non yet implemented consider parameters).
    int totalParameterSetSize_;

    //! Total number of estimated parameter values (excluding currently non yet implemented consider parameters).
    int estimatedParameterSetSize_;

    //! List of start indices and sizes (map keys) of estimated parameters.
    /*!
     * List of start indices and sizes (map keys) of estimated parameters, in order of vector
     * estimateInitialStateParameters_, followed by estimatedDoubleParameters_, followed by estimatedVectorParameters_.
     */
    std::vector< std::pair< int, int > > parameterIndices_;

    //! List of double parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< double > > > estimatedDoubleParameters_;

    //! List of vector parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedVectorParameters_;

    //! List of initial dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateInitialStateParameters_;

    //! List of initial single-arc dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateSingleArcInitialStateParameters_;

    //! List of initial multi-arc dynamical states that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateMultiArcInitialStateParameters_;

    //! Map of double parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr< EstimatableParameter< double > > > doubleParameters_;

    //! Map of vector parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters_;

    //! Map of initial dynamical states that are to be estimated, with start index in total parameter vector as key.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialStateParameters_;

    //! Size of linear constraint that is to be applied during estimation
    int totalConstraintSize_;

    //! Map containing all single-arc initial state parameter objects, with map key start index of parameter in total vector.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialSingleArcStateParameters_;

    //! Map containing all multi-arc initial state parameter objects, with map key start index of parameter in total vector.
    std::map< int, std::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialMultiArcStateParameters_;

};

//! Function to create a subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter
/*!
 *  Function to create a subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter.
 *  All non-dynamical state parameters are copied from the input to the output
 *  \param parametersToEstimate Total set of parameters
 *  \param getSingleArcParameters Boolean denoting whether single- or multi-arc initial state parameters are to be kept
 *  \return Subset of all estimated parameters, with either only single-arc or multi-arc initial state parameter
 */
template< typename StateScalarType >
std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > createEstimatableParameterSetArcSubSet(
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const bool getSingleArcParameters)
{
    if( getSingleArcParameters )
    {
        return std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                    parametersToEstimate->getEstimatedDoubleParameters( ),
                    parametersToEstimate->getEstimatedVectorParameters( ),
                    parametersToEstimate->getEstimatedSingleArcInitialStateParameters( ) );
    }
    else
    {
        return std::make_shared< estimatable_parameters::EstimatableParameterSet< StateScalarType > >(
                    parametersToEstimate->getEstimatedDoubleParameters( ),
                    parametersToEstimate->getEstimatedVectorParameters( ),
                    parametersToEstimate->getEstimatedMultiArcInitialStateParameters( )  );
    }
}

template< typename InitialStateParameterType >
void printEstimatableParameterEntries(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< int, std::shared_ptr<
            EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialStateParameters =
            estimatableParameters->getInitialStateParameters( );
    std::map< int, std::shared_ptr<
            EstimatableParameter< double > > > doubleParameters = estimatableParameters->getDoubleParameters( );
    std::map< int, std::shared_ptr<
            EstimatableParameter< Eigen::VectorXd > > > vectorParameters = estimatableParameters->getVectorParameters( );

    std::cout << "Parameter start index, Parameter definition" << std::endl;
    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< Eigen::Matrix<
         InitialStateParameterType, Eigen::Dynamic, 1 > > > >::const_iterator parameterIterator = initialStateParameters.begin( );
         parameterIterator != initialStateParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }

    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< double > > >::const_iterator
         parameterIterator = doubleParameters.begin( );
         parameterIterator != doubleParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }

    for( typename  std::map< int, std::shared_ptr<  EstimatableParameter< Eigen::VectorXd > > >::const_iterator
         parameterIterator = vectorParameters.begin( );
         parameterIterator != vectorParameters.end( ); parameterIterator++ )
    {
        std::cout << parameterIterator->first << ", " << parameterIterator->second->getParameterDescription( ) << std::endl;
    }
    std::cout << std::endl;
}

//! Function to get the list of names of bodies for which initial translational dynamical state is estimated.
/*!
 *  Function to get the list of names of bodies for which initial translational dynamical state is estimated.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return List of names of bodies for which initial state is estimated.
 */
template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithTranslationalStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

//! Function to get the list of names of bodies for which initial rotational dynamical state is estimated.
/*!
 *  Function to get the list of names of bodies for which initial rotational dynamical state is estimated.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return List of names of bodies for which initial rotational state is estimated.
 */
template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithRotationalStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithMassStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of bodies for which the translational state is estimated in a multi-arc fashion
/*!
 * Function to retrieve the list of bodies for which the translational state is estimated in a multi-arc fashion
 * \param estimatableParameters Full set of estimated parameters
 * \return List of parameters (with body names as keys) used for the multi-arc estimation of initial translational state
 */
template< typename InitialStateParameterType >
std::map< std::string, std::shared_ptr< EstimatableParameter<
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  > > >
getListOfBodiesWithTranslationalMultiArcStateToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< std::string, std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  > > > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state )
        {
            bodiesToEstimate[ initialDynamicalParameters.at( i )->getParameterName( ).second.first ] =
                    initialDynamicalParameters.at( i );
        }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of bodies for which an initial dynamical state is to be estimated
/*!
 * Function to retrieve the list of bodies for which an initial dynamical state is to be estimated
 * \param estimatableParameters Full set of estimated parameters
 * \return List of bodies for which an initial dynamical state is estimated.
 */
template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::string > > getListOfBodiesToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::map< propagators::IntegratedStateType, std::vector< std::string > > bodiesToEstimate;

    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )  ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            bodiesToEstimate[ propagators::translational_state ].push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state ) )
         {
             bodiesToEstimate[ propagators::rotational_state ].push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
         }
    }

    return bodiesToEstimate;
}

//! Function to retrieve the list of translational state parameters from full parameter list
/*!
 * Function to retrieve the list of translational state parameters from full parameter list
 * \param estimatableParameters Full set of estimated parameters
 * \return List of translational state parameters
 */
template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfTranslationalStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > translationalStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )  ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            translationalStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return translationalStateParameters;
}

//! Function to retrieve the list of rotational state parameters from full parameter list
/*!
 * Function to retrieve the list of rotational state parameters from full parameter list
 * \param estimatableParameters Full set of estimated parameters
 * \return List of rotational state parameters
 */
template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfRotationalStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > rotationalStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state )   )
        {
            rotationalStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return rotationalStateParameters;
}

template< typename InitialStateParameterType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > getListOfMassStateParametersToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > massStateParameters;

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state )   )
        {
            massStateParameters.push_back(  initialDynamicalParameters.at( i ) );
        }
    }

    return massStateParameters;
}

//! Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
/*!
 *  Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return Map containing dynamics type (key) and vector of pairs: list of bodies (first in pair) with reference point
 *  identifier (second in pair; empty if not relevant) for which given dynamics type is estimated.
 */
template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
getListOfInitialDynamicalStateParametersEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStateParametersEstimate;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state ) ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == arc_wise_initial_body_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::translational_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_rotational_body_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::rotational_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
        else if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_mass_state ) )
        {
            initialDynamicalStateParametersEstimate[ propagators::body_mass_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }
    }

    return initialDynamicalStateParametersEstimate;
}

//! Function to get initial state vector of estimated dynamical states.
/*!
 *  Function to get initial state vector of estimated dynamical states (i.e. presently estimated state at propagation
 *  start time.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return State vector of estimated dynamics at propagation start time.
 */
template< typename InitialStateParameterType = double >
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getInitialStateVectorOfBodiesToEstimate(
        const std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Initialize state vector.
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateVector =
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero(
                estimatableParameters->getInitialDynamicalStateParameterSize( ), 1 );

    int vectorSize = 0;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( isParameterDynamicalPropertyInitialState( initialDynamicalParameters.at( i )->getParameterName( ).first ) )
        {
            int currentParameterSize = initialDynamicalParameters.at( i )->getParameterSize( );
            initialStateVector.block( vectorSize, 0, currentParameterSize, 1 ) = initialDynamicalParameters.at( i )->getParameterValue( );

            vectorSize += currentParameterSize;
        }
    }

    return initialStateVector.block( 0, 0, vectorSize, 1 );
}

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_ESTIMATABLEPARAMETERS_H
