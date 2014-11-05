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
 *      130805    J. Geul           Copied file from cartesianStateDerivativeModel.h (by K. Kumar)
 *      130830    J. Geul           Free function acceleration to stateDerivative.
 *      131114    J. Geul           Added extra and expanded syncing calls.
 *      130106    J. Geul           Prepare to commit, added more comments.
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
 *      OrbitalStateDerivativeModel class without frame transformations and in conjunction with
 *      the CompositeStateDerivativeModel class; in general, the compiler should inline these
 *      functions so that there are not function table lookups performed. This can be improved in
 *      future perhaps by making use of the Boost::Lambda library.
 *
 *      In contrast to the cartesianStateDerivativeModel a 4 argument update function is used. 
 *      This is to also include the totalAccelerations and stateDerivative in the sync / body
 *      class. It however remains possible to bind a 2 argument updating function of a body, as
 *      the extra parameters are simply ignored.      
 *
 */

#ifndef TUDAT_ORBITAL_STATE_DERIVATIVE_MODEL_H
#define TUDAT_ORBITAL_STATE_DERIVATIVE_MODEL_H

#include <utility>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/stateDerivativeModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace tudat
{
namespace state_derivative_models
{

//! Perform no frame transformation.
/*!
 * Returns the input acceleration without performing any frame transformation. This dummy function
 * is used as the default argument of the constructor of the OrbitalStateDerivativeModel class.
 * \param acceleration.
 * \return acceleration (same as input).
 * \sa OrbitalStateDerivativeModel.
 */
template< typename AccelerationType >
inline AccelerationType transformNothing( const AccelerationType& acceleration )
{
    return acceleration;
}

//! Perform no update.
/*!
 * Performs no update; empty function body.  This dummy function can be used to pass a dummy update
 * function to the constructor of the OrbitalStateDerivativeModel class, when used in conjunction
 * with the CompositeStateDerivativeModel class, to ensure that the update is only performed once.
 * \param independentVariable Current independent variable value (e.g., time).
 * \param orbitalState Current composite state.
 * \sa OrbitalStateDerivativeModel, CompositeStateDerivativeModel.
 */
template< typename IndependentVariableType, typename OrbitalStateType >
inline void updateNothing( const IndependentVariableType independentVariable,
                           const OrbitalStateType& orbitalState )
{
}

//! Orbital state derivative model class.
/*!
 * Templated class that generates a Orbital state derivative model based on a list of
 * acceleration models provided by the user. The user can in addition provide a list of frame
 * transformations per acceleration model. The state derivative model is constructed by looping
 * through the list of acceleration models to generate a cumulative (summed) state derivative,
 * including any specified frame transformations, for a given current independent variable
 * (typically, but not necessarilly time) and Orbital state.
 * This class can be used in conjunction with the numerical integrators in Tudat,
 * to integrate equations of motion, constructed based on given acceleration models. The user is
 * also required to pass an update-function through the constructor that updates a user-defined
 * data repository where all dependent variables are set, which can be then accessed by the
 * acceleration models.
 * \tparam IndependentVariableType Data type for independent variable, e.g., time, (default is
 *          double).
 * \tparam OrbitalStateType Data type for Orbital state (default is Eigen::Vector6d).
 * \tparam AccelerationType Data type for Orbital acceleration (default is Eigen::Vector3d).
 * \tparam AccelerationModelType Type of acceleration models provided by user (default is
 *          acceleration models that return an AccelerationType acceleration).
 */
template< typename IndependentVariableType = double,
    typename OrbitalStateType = Eigen::VectorXd,
    typename AccelerationType = Eigen::Vector3d, 
    typename AccelerationModelType = basic_astrodynamics::AccelerationModel< AccelerationType > >
class OrbitalStateDerivativeModel
        : public StateDerivativeModel< IndependentVariableType, OrbitalStateType >
{
private:

    //! Typedef for a shared-pointer to an acceleration model.
    typedef boost::shared_ptr< AccelerationModelType > AccelerationModelPointer;

    //! Typedef for Orbital state-derivative type.
    typedef OrbitalStateType OrbitalStateDerivativeType;

    //! Typedef for pointer to a set-function that updates independent variable and state data.
    typedef boost::function< void ( const IndependentVariableType, const OrbitalStateType&, 
	const OrbitalStateType&, const AccelerationType& ) > IndependentVariableAndStateUpdateFunction;

    //! Typedef for a reference frame transformation function.
    typedef boost::function< AccelerationType ( const AccelerationType& ) >
        ReferenceFrameTransformationFunction;

    //! Typedef for a function (eg planetary equations) that relates accelerations to
    //! the state derivative
    typedef boost::function< OrbitalStateDerivativeType ( const OrbitalStateType&, 
	const AccelerationType& ) > AccelerationsToStateDerivativeFunction;

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
     * \param accelerationsToStateDerivativeFunction Pointer to a function to function that maps
     *          accelerations to state state derivative.
     * \sa transformNothing().
     */
    OrbitalStateDerivativeModel(
            const AccelerationModelPointerVector& aListOfAccelerations,
            const IndependentVariableAndStateUpdateFunction
            anIndependentVariableAndStateUpdateFunction,
            const AccelerationsToStateDerivativeFunction 
	    accelerationsToStateDerivativeFunction );

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
     * \param accelerationsToStateDerivativeFunction boost function that relates the accelerations in
     *          aListOfAccelerations to the derivative of the state, e.g. Gauss planetary equations.
     */
    OrbitalStateDerivativeModel( const ListOfAccelerationFrameTransformationPairs&
	aListOfAccelerationFrameTransformationPairs,
	const IndependentVariableAndStateUpdateFunction	anIndependentVariableAndStateUpdateFunction,
	const AccelerationsToStateDerivativeFunction accelerationsToStateDerivativeFunction )
        : listOfAccelerationFrameTransformationPairs(
              aListOfAccelerationFrameTransformationPairs ),
        updateIndependentVariableAndState( anIndependentVariableAndStateUpdateFunction ),
        accelerationsToStateDerivativeFunction( accelerationsToStateDerivativeFunction )
    { }

    //! Compute Orbital state derivative.
    /*!
     * Computes the Orbital state derivative based on the list of acceleration models and
     * associated frame transformations provided through the constructor.
     * \param independentVariable Current independent variable value.
     * \param orbitalState Current Orbital state.
     * \return Computed Orbital State derivative vector.
     */
    OrbitalStateDerivativeType computeStateDerivative(
            const IndependentVariableType independentVariable,
            const OrbitalStateType& orbitalState );

protected:

private:

    //! List of acceleration model/frame transformation pairs.
    /*!
     * List of pairs of shared-pointers to acceleration model and associated lists of reference
     * frame transformation functions. This list is used by computeStateDerivative() to compute the
     * overall Orbital state derivative.
     * \sa computeStateDerivative().
     */
    ListOfAccelerationFrameTransformationPairs listOfAccelerationFrameTransformationPairs;

    //! Pointer to update function to update independent variable and state.
    /*!
     * Pointer to a function that updates user-defined data repository containing independent
     * variable and state data to the current values.
     */
    const IndependentVariableAndStateUpdateFunction updateIndependentVariableAndState;

    //! Pointer to update function to update independent variable and state.
    /*!
     * Pointer to a function that updates user-defined data repository containing independent
     * variable and state data to the current values.
     */
    const AccelerationsToStateDerivativeFunction accelerationsToStateDerivativeFunction;
    
};

//! Constructor taking list of acceleration models, and pointer to a function to update independent
template< typename IndependentVariableType, typename OrbitalStateType, typename AccelerationType,
          typename AccelerationModelType >
OrbitalStateDerivativeModel< IndependentVariableType, OrbitalStateType,
AccelerationType, AccelerationModelType >::OrbitalStateDerivativeModel(
        const AccelerationModelPointerVector& aListOfAccelerations,
        const IndependentVariableAndStateUpdateFunction anIndependentVariableAndStateUpdateFunction,
	const AccelerationsToStateDerivativeFunction accelerationsToStateDerivativeFunction )
    : updateIndependentVariableAndState( anIndependentVariableAndStateUpdateFunction ),
      accelerationsToStateDerivativeFunction( accelerationsToStateDerivativeFunction )
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

//! Compute Orbital state derivative.
template< typename IndependentVariableType, typename OrbitalStateType, typename AccelerationType,
          typename AccelerationModelType >
OrbitalStateType OrbitalStateDerivativeModel< IndependentVariableType, OrbitalStateType,
AccelerationType, AccelerationModelType >::computeStateDerivative(
        const IndependentVariableType independentVariable,
        const OrbitalStateType& orbitalState)
{

    // Declare total acceleration.
    AccelerationType totalAcceleration = AccelerationType::Zero( );

    // Declare state derivative of the same size as state.
    OrbitalStateType stateDerivative = orbitalState * 0;

    // Update data for first time. Although totalAcceleration and
    // stateDerivative are zero at the moment, the idependentVariable
    // and orbitalState could be used by the different models.
    // Therefore it makes sense to update it before and after.
    updateIndependentVariableAndState( independentVariable, orbitalState, 
	stateDerivative, totalAcceleration );

    // Loop through list of acceleration/frame transformation pairs.
    for ( unsigned int i = 0; i < listOfAccelerationFrameTransformationPairs.size( ); i++ )
    {
        // Update class members for current acceleration model.
        listOfAccelerationFrameTransformationPairs.at( i ).first->updateMembers( );

        // Get acceleration for current acceleration model.
        AccelerationType acceleration = listOfAccelerationFrameTransformationPairs.at( i )
                .first->getAcceleration( );

        // Loop through list of frame transformations and apply to computed acceleration.
        for ( unsigned j = 0;
              j < listOfAccelerationFrameTransformationPairs.at( i ).second.size( ); j++ )
        {
            acceleration = listOfAccelerationFrameTransformationPairs.at( i ).second.at( j )(
                        acceleration );
        }

        // Add transformed acceleration to total acceleration
        totalAcceleration += acceleration;
    }
    stateDerivative = 
	accelerationsToStateDerivativeFunction( orbitalState, totalAcceleration );

    // Update data a second time, so also stateDerivative and
    // totalAcceleration are in sync as well.
    updateIndependentVariableAndState( independentVariable, orbitalState, 
	stateDerivative, totalAcceleration );

    // Return assembled state derivative.
    return stateDerivative;
}

//! Typedef for shared-pointer to OrbitalStateDerivativeModel object.
typedef OrbitalStateDerivativeModel< > OrbitalStateDerivativeModelType;

//! Typedef for shared-pointer to OrbitalStateDerivativeModel object.
typedef boost::shared_ptr< OrbitalStateDerivativeModelType > OrbitalStateDerivativeModelPointer;

} // namespace state_derivative_models
} // namespace tudat

#endif // TUDAT_ORBITAL_STATE_DERIVATIVE_MODEL_H
