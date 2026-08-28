/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTRAINED_TRANSLATIONAL_STATE_H
#define TUDAT_CONSTRAINED_TRANSLATIONAL_STATE_H

#include <Eigen/Core>
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"

namespace tudat
{
namespace estimatable_parameters
{

//! Constrained initial translational state parameter for Line of Variations (LOV) differential correction.
template< typename InitialStateParameterType = double >
class ConstrainedTranslationalStateParameter : public InitialTranslationalStateParameter< InitialStateParameterType >
{
public:
    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state. Explicitly forwards to the base class constructor with
     * parameterType set to constrained_initial_body_state, so that this parameter correctly identifies itself as distinct
     * from a plain InitialTranslationalStateParameter (e.g. via getParameterName()).
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTranslationalState Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    ConstrainedTranslationalStateParameter( const std::string& associatedBody,
                                            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialTranslationalState,
                                            const std::string& centralBody = "SSB",
                                            const std::string& frameOrientation = "ECLIPJ2000" ):
        InitialTranslationalStateParameter< InitialStateParameterType >( associatedBody,
                                                                         initialTranslationalState,
                                                                         centralBody,
                                                                         frameOrientation,
                                                                         constrained_initial_body_state )
    {}

    //! Function to get the size of the constraint
    /*!
     * \return Size of the constraint (1 for the orthogonality condition).
     */
    int getConstraintSize( ) override
    {
        return 1;
    }

    //! Function to get the multiplier of the parameter constraint
    /*!
     * \return Multiplier of the parameter constraint (transposed weak direction).
     * NOTE: The spelling 'Multipler' intentionally matches the base EstimatableParameter class.
     */
    Eigen::MatrixXd getConstraintStateMultipler( ) override
    {
        // Initialize 1xN matrix (N=6 for translational state) for the constraint multiplier
        Eigen::MatrixXd constraintMultiplier = Eigen::MatrixXd::Zero( 1, this->getParameterSize( ) );

        // Populate the row with the transposed weak direction
        constraintMultiplier.row( 0 ) = weakDirection_.transpose( );

        return constraintMultiplier;
    }

    //! Function to get the right-hand side of the parameter constraint
    /*!
     * \return Right-hand side of the parameter constraint (0 for exact orthogonality).
     */
    Eigen::VectorXd getConstraintRightHandSide( ) override
    {
        // For exact orthogonality, dot product = 0
        return Eigen::VectorXd::Zero( 1 );
    }

    //! Function to dynamically set the weak direction used for the constraint
    /*!
     * \param newWeakDirection New weak direction vector evaluated at the current state.
     */
    void setWeakDirection( const Eigen::Matrix< double, 6, 1 >& newWeakDirection )
    {
        weakDirection_ = newWeakDirection;
    }

    //! Function to retrieve the currently set weak direction
    /*!
     * \return Current weak direction vector.
     */
    Eigen::Matrix< double, 6, 1 > getWeakDirection( ) const
    {
        return weakDirection_;
    }

private:
    //! Weak direction vector restricting the differential correction.
    Eigen::Matrix< double, 6, 1 > weakDirection_ = Eigen::Matrix< double, 6, 1 >::Zero( );
};

}  // namespace estimatable_parameters
}  // namespace tudat

#endif  // TUDAT_CONSTRAINED_TRANSLATIONAL_STATE_H