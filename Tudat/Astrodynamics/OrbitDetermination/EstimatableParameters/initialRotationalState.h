/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INITIALROTATIONALSTATE_H
#define TUDAT_INITIALROTATIONALSTATE_H

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of an initial rotational state.
template< typename InitialStateParameterType = double >
class InitialRotationalStateParameter: public EstimatableParameter<
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:

    //! Constructor, sets initial value of rotational state.
    /*!
     * Constructor, sets initial value of rotational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialRotationalState Current value of initial state (w.r.t. centralBody).
     * \param inertiaTensorFunction Function returning the value of the inertia tensor.
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialRotationalStateParameter(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialRotationalState,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >(
            initial_rotational_body_state, associatedBody ),
        initialRotationalState_( initialRotationalState ), centralBody_( centralBody ),
        frameOrientation_( frameOrientation ), inertiaTensorFunction_( inertiaTensorFunction )
    { }

    //! Function to get the current value of initial state w.r.t. centralBody.
    /*!
     * Function to get the current value of initial state w.r.t. centralBody.
     * \return The current value of initial state w.r.t. centralBody.
     */
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialRotationalState_;
    }

    //! Function to reset the current value of initial state w.r.t. centralBody.
    /*!
     * Function to reset the current value of initial state w.r.t. centralBody.
     * \param parameterValue The new value of initial state w.r.t. centralBody.
     */
    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialRotationalState_ = parameterValue;
        initialRotationalState_.segment( 0, 4 ).normalize( );
    }

    //! Function to retrieve the size of the parameter (always 6).
    /*!
     *  Function to retrieve the size of the parameter (always 6).
     *  \return Size of parameter value (always 6).
     */
    int getParameterSize( )
    {
        return 7;
    }

    //! Function to get the name of the body w.r.t. which the initial state is to be estimated.
    /*!
     * Function to get the name of the body w.r.t. which the initial state is to be estimated.
     * \return Name of the body w.r.t. which the initial state is to be estimated.
     */
    std::string getCentralBody( )
    {
        return centralBody_;
    }

    //! Function to retrieve size of constraint to be applied on parameter
    /*!
     * Function to retrieve size of constraint to be applied on parameter, equal to 1 for this parameter (quaternion norm).
     * \return Size of constraint to be applied on parameter
     */
    int getConstraintSize( )
    {
        return 1;
    }

    //! Function to retrieve multiplier for parameter linear constraint
    /*!
     * Function to retrieve multiplier for parameter linear constraint.
     * \return Multiplier for parameter linear constraint
     */
    Eigen::MatrixXd getConstraintStateMultipler( )
    {
        Eigen::MatrixXd constraintsMatrix = Eigen::MatrixXd::Zero( 1, 7 );
        constraintsMatrix.block( 0, 0, 1, 4 ) = initialRotationalState_.segment( 0, 4 ).template cast< double >( ).transpose( );
        return constraintsMatrix;
    }

    //! Function to retrieve right-hand side for parameter linear constraint
    /*!
     * Function to retrieve right-hand side for parameter linear constraint.
     * \return Right-hand side for parameter linear constraint
     */
    Eigen::VectorXd getConstraintRightHandSide( )
    {
        return Eigen::VectorXd::Zero( 1 );
    }

    //! Function to retrieve function that returns the current inertia tensor
    /*!
     * Function to retrieve function that returns the current inertia tensor
     * \return Function that returns the current inertia tensor
     */
    std::function< Eigen::Matrix3d( ) > getBodyInertiaTensorFunction( )
    {
        return inertiaTensorFunction_;
    }

private:

    //! Current value of initial state (w.r.t. centralBody)
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialRotationalState_;

    //!  Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

    //! Function that returns the current inertia tensor
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_INITIALROTATIONALSTATE_H
