/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120911    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TEST_STATE_DERIVATIVE_MODELS_H
#define TUDAT_TEST_STATE_DERIVATIVE_MODELS_H

#include <cmath>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/StateDerivativeModels/stateDerivativeModel.h"

namespace tudat
{
namespace unit_tests
{

//! Typedef for a 2x3 matrix state.
typedef Eigen::Matrix< double, 2, 3 > Matrix2x3d;

//! Typedef for a Vector5d state.
typedef Eigen::Matrix< double, 5, 1 > Vector5d;

//! Time (independent variable).
static double time;

//! Vector2d state.
static Eigen::Vector2d vector2dState;

//! Get Vector2d state.
Eigen::Vector2d getVector2dState( ) { return vector2dState; }

//! Vector3d state.
static Eigen::Vector3d vector3dState;

//! Get Vector3d state.
Eigen::Vector3d getVector3dState( ) { return vector3dState; }

//! Matrix2d state.
static Eigen::Matrix2d matrix2dState;

//! Get Matrix2d state.
Eigen::Matrix2d getMatrix2dState( ) { return matrix2dState; }

//! Compute Vector2d state derivative.
/*!
 * Computes the state derivative of the Vector2d state using a test function,
 * given the independent variable and Vector2d state as input.
 * \param independentVariable Independent variable (time).
 * \param state Vector2d state.
 * \return State derivative of the Vector2d state.
 */
Eigen::Vector2d computeVector2dStateDerivative( const double independentVariable,
                                                const Eigen::Vector2d& state )
{
    // Declare state derivative.
    Eigen::Vector2d stateDerivative;

    // Compute state derivative.
    stateDerivative( 0 ) = state( 0 ) * state( 0 ) - independentVariable;
    stateDerivative( 1 ) = std::sqrt( std::fabs( state( 1 ) ) )
            - independentVariable * independentVariable;

    // Return state derivative.
    return stateDerivative;
}

//! Compute Vector3d state derivative.
/*!
 * Computes the state derivative of the Vector3d state using a test function,
 * given the independent variable and Vector3d state as input.
 * \param independentVariable Independent variable (time).
 * \param state Vector3d state.
 * \return State derivative of the Vector3d state.
 */
Eigen::Vector3d computeVector3dStateDerivative( const double independentVariable,
                                                const Eigen::Vector3d& state )
{
    // Declare state derivative.
    Eigen::Vector3d stateDerivative;

    // Compute state derivative.
    stateDerivative( 0 ) = state( 0 ) * state( 1 );
    stateDerivative( 1 ) = state( 2 ) /  independentVariable;
    stateDerivative( 2 ) = state( 2 ) * state( 2 );

    // Return state derivative.
    return stateDerivative;
}

//! Update matrix data.
/*!
 * Updates data stored in data repository for matrix-based composite state. The data repository in
 * this case consists of the time, vector2dState, and matrix2dstate (static) variables available in
 * this file. Essentially, this function disassembles the composite state and updates these
 * variables to the latest values. The indices used in this function must correspond to the indices
 * provided to the state derivative model, e.g., the CompositeStateDerivativeModel object.
 * \param independentVariable Independent variable (time)
 * \param compositeState Composite state as 2x3 matrix.
 * \sa CompositeStateDerivativeModel.
 */
void updateMatrixData( const double independentVariable, const Matrix2x3d& compositeState )
{
    // Update time.
    time = independentVariable;

    // Update Vector2d state.
    vector2dState = compositeState.block( 0, 0, 2, 1 );

    // Update Matrix2d state.
    matrix2dState = compositeState.block( 0, 1, 2, 2 );
}

//! Update vector data.
/*!
 * Updates data stored in data repository for vector-based composite state. The data repository in
 * this case consists of the time, vector3dState, and vector3dState (static) variables available in
 * this file. Essentially, this function disassembles the composite state and updates these
 * variables to the latest values. The indices used in this function must correspond to the indices
 * provided to the state derivative model, e.g., the CompositeStateDerivativeModel object.
 * \param independentVariable Independent variable (time)
 * \param compositeState Composite state as Vector5d.
 * \sa CompositeStateDerivativeModel.
 */
void updateVectorData( const double independentVariable, const Vector5d& compositeState )
{
    // Update time.
    time = independentVariable;

    // Update Vector3d state.
    vector3dState = compositeState.segment( 0, 3 );

    // Update Vector2d state.
    vector2dState = compositeState.segment( 3, 2 );
}

//! Coupled Matrix2d state derivative model.
/*!
 * This class is used to compute the state derivative of the Matrix2d state. This class serves as
 * an example of how the CompositeStateDerivative model class can be used to construct composite
 * state derivative models with coupling between the part states. In this case, the class depends
 * on the Vector2d state, which is fetched using the function pointer provided through the
 * constructor. This is a derived class of the StateDerivativeModel class.
 */
class CoupledMatrix2dStateDerivativeModel
        : public state_derivative_models::StateDerivativeModel< double, Eigen::Matrix2d >
{
private:

    //! Typedef for Eigen::Vector2d returning function.
    typedef boost::function< Eigen::Vector2d( ) > Vector2dReturningFunction;

public:

    //! Constructor taking a pointer to a Vector2d-returning function.
    /*!
     * Constructor taking a pointer to a Vector2d-returning function, which is used to fetch the
     * current Vector2d, necessary to compute the state derivative of the Matrix2d state.
     * \param aVector2dReturningFunction Pointer to a function returning the current Vector2d
     *          state.
     */
    CoupledMatrix2dStateDerivativeModel(
            const Vector2dReturningFunction aVector2dReturningFunction );

    //! Compute state derivative.
    /*!
     * Computes state derivative model for the Matrix2d state. This state derivative model
     * illustrates how the CompositeStateDerivative model class can be used to construct composite
     * state derivative models that are coupled.
     * \param independentVariable Independent variable (time).
     * \param state Matrix2d state.
     * \return The state derivative of the Matrix2d state.
     */
    Eigen::Matrix2d computeStateDerivative(
            const double independentVariable, const Eigen::Matrix2d& state );

protected:

private:

    //! Pointer to function returning Vector2d state.
    Vector2dReturningFunction getCurrentVector2dState;

    //! Current Vector2d state stored locally.
    Eigen::Vector2d currentVector2d;
};

//! Constructor taking a pointer to a Vector2d-returning function.
CoupledMatrix2dStateDerivativeModel::CoupledMatrix2dStateDerivativeModel(
        const Vector2dReturningFunction aVector2dReturningFunction )
    : getCurrentVector2dState( aVector2dReturningFunction )
{
    // Update locally stored Vector2d state to the current value.
    currentVector2d = getCurrentVector2dState( );
}

//! Compute state derivative.
Eigen::Matrix2d CoupledMatrix2dStateDerivativeModel::computeStateDerivative(
        const double independentVariable, const Eigen::Matrix2d& state )
{
    // Declare state derivative.
    Eigen::Matrix2d stateDerivative;

    // Update locally stored Vector2d state to the current value.
    currentVector2d = getCurrentVector2dState( );

    // Compute state derivative.
    stateDerivative( 0, 0 ) = state( 0, 0 ) - independentVariable + currentVector2d( 0 );
    stateDerivative( 1, 0 ) = state( 0, 0 ) * state( 1, 0 )
            - std::sqrt( std::fabs( currentVector2d( 1 ) ) );
    stateDerivative( 0, 1 ) = std::sqrt( std::fabs( independentVariable ) );
    stateDerivative( 1, 1 ) = state( 1, 1 ) / state( 0, 0 )  - independentVariable
            * currentVector2d( 1 ) * currentVector2d( 1 );

    // Return state derivative.
    return stateDerivative;
}

//! Typedef for a shared-pointer to a coupled matrix3d state derivative model.
typedef boost::shared_ptr< CoupledMatrix2dStateDerivativeModel >
CoupledMatrix2dStateDerivativeModelPointer;


} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_STATE_DERIVATIVE_MODELS_H
