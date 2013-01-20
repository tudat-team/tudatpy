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
 *      120817    K. Kumar          File created.
 *      120817    K. Kumar          Templated class and completed Doxygen documentation.
 *      120823    K. Kumar          Rewrote as derived class of templated StateDerivativeModel base
 *                                  class.
 *      120824    K. Kumar          Split off functionality for part model, to be used with
 *                                  CompositeStateDerivativeModel
 *                                  (cartesianStateDerivativePartModel.h).
 *      120913    K. Kumar          Updated model to work with and without
 *                                  CompositeStateDerivativeModel; added frame transformation
 *                                  functionality.
 *
 *    References
 *
 *    Notes
 *      Even though this class is fully templatized to work with generic data types for the
 *      composite state (derivative), the use of the Zero(), rows(), and segment() functions are
 *      Eigen-specific. At present, these functions must be available in any other data types
 *      used.
 *
 *      Currently, the transformNothing() and updateNothing() provide the means to use the
 *      CartesianStateDerivativeModel class without frame transformations and in conjunction with
 *      the CompositeStateDerivativeModel class; in general, the compiler should inline these
 *      functions so that there are not function table lookups performed. This can be improved in
 *      future perhaps by making use of the Boost::Lambda library.
 *
 */

#ifndef TUDAT_CARTESIAN_STATE_DERIVATIVE_MODEL_H
#define TUDAT_CARTESIAN_STATE_DERIVATIVE_MODEL_H

#include <utility>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/StateDerivativeModels/stateDerivativeModel.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace state_derivative_models
{

//! Perform no frame transformation.
/*!
 * Returns the input acceleration without performing any frame transformation. This dummy function
 * is used as the default argument of the constructor of the CartesianStateDerivativeModel class.
 * \param Input acceleration.
 * \return Output acceleration (same as input).
 * \sa CartesianStateDerivativeModel.
 */
template< typename AccelerationType >
inline AccelerationType transformNothing( const AccelerationType& acceleration )
{
    return acceleration;
}

//! Perform no update.
/*!
 * Performs no update; empty function body.  This dummy function can be used to pass a dummy update
 * function to the constructor of the CartesianStateDerivativeModel class, when used in conjunction
 * with the CompositeStateDerivativeModel class, to ensure that the update is only performed once.
 * \param independentVariable Current independent variable value (e.g., time).
 * \param compositeState Current composite state.
 * \sa CartesianStateDerivativeModel, CompositeStateDerivativeModel.
 */
template< typename IndependentVariableType, typename CartesianStateType >
inline void updateNothing( const IndependentVariableType independentVariable,
                           const CartesianStateType& cartesianState )
{
}

//! Cartesian state derivative model class.
/*!
 * Templated class that generates a Cartesian state derivative model based on a list of
 * acceleration models provided by the user. The user can in addition provide a list of frame
 * transformations per acceleration model. The state derivative model is constructed by looping
 * through the list of acceleration models to generate a cumulative (summed) state derivative,
 * including any specified frame transformations, for a given current independent variable
 * (typically, but not necessarilly time) and Cartesian state.
 * This class can be used in conjunction with the numerical integrators in Tudat,
 * to integrate equations of motion, constructed based on given acceleration models. The user is
 * also required to pass an update-function through the constructor that updates a user-defined
 * data repository where all dependent variables are set, which can be then accessed by the
 * acceleration models.
 * \tparam IndependentVariableType Data type for independent variable, e.g., time, (default is
 *          double).
 * \tparam CartesianStateType Data type for Cartesian state (default is Eigen::Vector6d).
 * \tparam AccelerationType Data type for Cartesian acceleration (default is Eigen::Vector3d).
 * \tparam AccelerationModelType Type of acceleration models provided by user (default is
 *          acceleration models that return an AccelerationType acceleration).
 */
template< typename IndependentVariableType = double,
          typename CartesianStateType = basic_mathematics::Vector6d,
          typename AccelerationType = Eigen::Vector3d, typename AccelerationModelType
          = basic_astrodynamics::AccelerationModel< AccelerationType > >
class CartesianStateDerivativeModel
        : public StateDerivativeModel< IndependentVariableType, CartesianStateType >
{
private:

    //! Typedef for a shared-pointer to an acceleration model.
    typedef boost::shared_ptr< AccelerationModelType > AccelerationModelPointer;

    //! Typedef for Cartesian state-derivative type.
    typedef CartesianStateType CartesianStateDerivativeType;

    //! Typedef for pointer to a set-function that updates independent variable and state data.
    typedef boost::function< void ( const IndependentVariableType, const CartesianStateType& ) >
    IndependentVariableAndStateUpdateFunction;

    //! Typedef for a reference frame transformation function.
    typedef boost::function< AccelerationType ( const AccelerationType& ) >
    ReferenceFrameTransformationFunction;

public:

    //! Typedef for a vector of shared-pointers to acceleration models.
    typedef std::vector< AccelerationModelPointer > AccelerationModelPointerVector;

    //! Typedef for list of reference frame transformation functions.
    typedef std::vector< ReferenceFrameTransformationFunction >
    ListOfReferenceFrameTransformations;

    //! Typedef for acceleration model/frame transformation pair.
    typedef std::pair< AccelerationModelPointer, ListOfReferenceFrameTransformations >
    AccelerationFrameTransformationPair;

    //! Typedef for list of acceleration model/frame transformation pairs.
    typedef std::vector< AccelerationFrameTransformationPair >
    ListOfAccelerationFrameTransformationPairs;

    //! Constructor taking list of acceleration models, and pointer to a function to update
    //! independent variable and state.
    /*!
     * Constructor taking list of acceleration models, and pointer to a function to update
     * independent variable and state, held externally in user-defined data repository. This
     * constructor ensures that the default dummy frame transformation transformNothing() is used.
     * \param aListOfAccelerations List of acceleration models.
     * \param anIndependentVariableAndStateUpdateFunction Pointer to a function to update
     *          independent variable and state.
     * \sa transformNothing().
     */
    CartesianStateDerivativeModel(
            const AccelerationModelPointerVector& aListOfAccelerations,
            const IndependentVariableAndStateUpdateFunction
            anIndependentVariableAndStateUpdateFunction );

    //! Constructor taking list of acceleration models with associated frame transformations, and
    //! pointer to a function to update independent variable and state.
    /*!
     * Constructor taking list of acceleration models, with a list of frame transformations per
     * acceleration model, and pointer to a function to update independent variable and state, held
     * externally in user-defined data repository.
     * \param aListOfAccelerationFrameTransformationPairs List of acceleration models with
     *          associated frame transformations. Note that the transformations are applied in the
     *          order in which they are provide in the in vectors, i.e. reading the frame
     *          transformation equation involving multiple transformations backwards.
     * \param anIndependentVariableAndStateUpdateFunction Pointer to a function to update
     *          independent variable and state.
     */
    CartesianStateDerivativeModel( const ListOfAccelerationFrameTransformationPairs&
                                   aListOfAccelerationFrameTransformationPairs,
                                   const IndependentVariableAndStateUpdateFunction
                                   anIndependentVariableAndStateUpdateFunction )
        : listOfAccelerationFrameTransformationPairs(
              aListOfAccelerationFrameTransformationPairs ),
          updateIndependentVariableAndState( anIndependentVariableAndStateUpdateFunction )
    { }

    //! Compute Cartesian state derivative.
    /*!
     * Computes the Cartesian state derivative based on the list of acceleration models and
     * associated frame transformations provided through the constructor.
     * \param independentVariable Current independent variable value.
     * \param cartesianState Current Cartesian state.
     * \return Computed Cartesian State derivative vector.
     */
    CartesianStateDerivativeType computeStateDerivative(
            const IndependentVariableType independentVariable,
            const CartesianStateType& cartesianState );

protected:

private:

    //! List of acceleration model/frame transformation pairs.
    /*!
     * List of pairs of shared-pointers to acceleration model and associated lists of reference
     * frame transformation functions. This list is used by computeStateDerivative() to compute the
     * overall Cartesian state derivative.
     * \sa computeStateDerivative().
     */
    ListOfAccelerationFrameTransformationPairs listOfAccelerationFrameTransformationPairs;

    //! Pointer to update function to update independent variable and state.
    /*!
     * Pointer to a function that updates user-defined data repository containing independent
     * variable and state data to the current values.
     */
    const IndependentVariableAndStateUpdateFunction updateIndependentVariableAndState;
};

//! Constructor taking list of acceleration models, and pointer to a function to update independent
//! variable and state.
template< typename IndependentVariableType, typename CartesianStateType, typename AccelerationType,
          typename AccelerationModelType >
CartesianStateDerivativeModel< IndependentVariableType, CartesianStateType,
AccelerationType, AccelerationModelType >::CartesianStateDerivativeModel(
        const AccelerationModelPointerVector& aListOfAccelerations,
        const IndependentVariableAndStateUpdateFunction
        anIndependentVariableAndStateUpdateFunction )
    : updateIndependentVariableAndState( anIndependentVariableAndStateUpdateFunction )
{
    // Loop through list of acceleration models.
    for ( unsigned int i = 0; i < aListOfAccelerations.size( ); i++ )
    {
        // Assign reference frame transformation list with dummy function.
        ListOfReferenceFrameTransformations listOfReferenceFrameTransformations
                = boost::assign::list_of( &transformNothing< AccelerationType > );

        // Add pair of acceleration model and list of associated reference frame
        // transformations to list.
        listOfAccelerationFrameTransformationPairs.push_back(
                    std::make_pair( aListOfAccelerations.at( i ),
                                    listOfReferenceFrameTransformations ) );
    }
}

//! Compute Cartesian state derivative.
template< typename IndependentVariableType, typename CartesianStateType, typename AccelerationType,
          typename AccelerationModelType >
CartesianStateType CartesianStateDerivativeModel< IndependentVariableType, CartesianStateType,
AccelerationType, AccelerationModelType >::computeStateDerivative(
        const IndependentVariableType independentVariable,
        const CartesianStateType& cartesianState )
{
    // Update data.
    updateIndependentVariableAndState( independentVariable, cartesianState );

    // Declare Cartesian state derivative size.
    unsigned int stateDerivativeSize = cartesianState.rows( );

    // Declare Cartesian state derivative of the same size as Cartesian state.
    CartesianStateDerivativeType cartesianStateDerivative
            = CartesianStateDerivativeType::Zero( stateDerivativeSize );

    // Set derivative of position components to current Cartesian velocity.
    cartesianStateDerivative.segment( 0, stateDerivativeSize / 2 )
            = cartesianState.segment( stateDerivativeSize / 2, stateDerivativeSize / 2 );

    // Loop through list of acceleration/frame transformation pairs.
    for ( unsigned int i = 0; i < listOfAccelerationFrameTransformationPairs.size( ); i++ )
    {
        // Update class members for current acceleration model.
        listOfAccelerationFrameTransformationPairs.at( i ).first->updateMembers( );

        // Get acceleration for current acceleration model.
        AccelerationType acceleration = listOfAccelerationFrameTransformationPairs.at( i )
                .first->getAcceleration( );

        // Loop througb list of frame transformations and apply to computed acceleration.
        for ( unsigned j = 0;
              j < listOfAccelerationFrameTransformationPairs.at( i ).second.size( ); j++ )
        {
            acceleration = listOfAccelerationFrameTransformationPairs.at( i ).second.at( j )(
                        acceleration );
        }

        // Add transformed acceleration to state derivative.
        cartesianStateDerivative.segment( stateDerivativeSize / 2, stateDerivativeSize / 2 )
                += acceleration;
    }

    // Return assembled state derivative.
    return cartesianStateDerivative;
}

//! Typedef for a 6D Cartesian state derivative model.
typedef CartesianStateDerivativeModel< > CartesianStateDerivativeModel6d;

//! Typedef for shared-pointer to CartesianStateDerivativeModel6d object.
typedef boost::shared_ptr< CartesianStateDerivativeModel6d >
CartesianStateDerivativeModel6dPointer;

//! Typedef for a 4D Cartesian state derivative model.
typedef CartesianStateDerivativeModel< double, Eigen::Vector4d, Eigen::Vector2d >
CartesianStateDerivativeModel4d;

//! Typedef for shared-pointer to CartesianStateDerivativeModel4d object.
typedef boost::shared_ptr< CartesianStateDerivativeModel4d >
CartesianStateDerivativeModel4dPointer;

} // namespace state_derivative_models
} // namespace tudat

#endif // TUDAT_CARTESIAN_STATE_DERIVATIVE_MODEL_H
