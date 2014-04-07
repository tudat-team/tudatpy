/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      120817    K. Kumar          Completed Doxygen documentation.
 *      120913    K. Kumar          Rewrote class to work with generic matrices.
 *      140127    E. Brandon        Corrected doxygen documentation.
 *
 *    References
 *
 *    Notes
 *      Even though this class is fully templatized to work with generic data types for the
 *      composite state (derivative), the use of the Zero(), rows(), cols(), and block() functions
 *      are Eigen-specific. At present, these functions must be available in any other data types
 *      used.
 *
 */

#ifndef TUDAT_COMPOSITE_STATE_DERIVATIVE_MODEL_H
#define TUDAT_COMPOSITE_STATE_DERIVATIVE_MODEL_H

#include <map>
#include <utility>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/StateDerivativeModels/stateDerivativeModel.h"

namespace tudat
{
namespace state_derivative_models
{

//! Composite state derivative model class.
/*!
 * Templated class that generates a composite state derivative model based on a map of pointers to
 * state-derivative-model functions provided by the user. The map indexes the composite state
 * so that the correct state derivative functions are called for each part state. This can be used,
 * for instance, to numerically integrate the state of multiple satellites around a central body,
 * each subject to their own state-derivative model, or the state of a satellite and the associated
 * State Transition Matrix. The class has been set up in a general fashion such that it can be used
 * in many other simulation scenarios too.
 * \tparam IndependentVariableType Data type for independent variable, e.g., time, (default is
 *          double).
 * \tparam CompositeStateType Data type for composite state (default is Eigen::MatrixXd).
 * \tparam PartStateType Data type for individual part states that are part of the composite state
 *          (default is same data type as CompositeStateType).
 */
template< typename IndependentVariableType = double, typename CompositeStateType = Eigen::MatrixXd,
          typename PartStateType = CompositeStateType >
class CompositeStateDerivativeModel
        : public StateDerivativeModel< IndependentVariableType, CompositeStateType >
{
private:

    //! Typedef for composite state derivative type.
    typedef CompositeStateType CompositeStateDerivativeType;

    //! Typedef for part state derivative type.
    typedef PartStateType PartStateDerivativeType;

    //! Typedef for block indices used to define part states in composite state matrix
    //! (1st index=start row, 2nd index = start column, 3rd index = number of rows,
    //! 4th index = number of columns).
    typedef boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int >
    StateBlockIndices;

    //! Typedef for segment indices used to define part states in composite state vector
    //! (1st index=start, 2nd index=length).
    typedef std::pair< unsigned int, unsigned int > StateSegmentIndices;

    //! Typedef for a pointer to a function that updates independent variable and composite state.
    typedef boost::function< void( const IndependentVariableType, const CompositeStateType& ) >
    IndependentVariableAndStateUpdateFunction;

    //! Typedef for a pointer to a function that evaluates the state derivative corresponding to a
    //! part of the composite state derivative.
    typedef boost::function< PartStateDerivativeType(
            const IndependentVariableType, const PartStateType& ) > PartStateDerivativeFunction;

public:

    //! Typedef for state-derivative model map for matrix composite state.
    /*!
     * Typedef for state-derivative model map, that maps part states (matrices) in the composite
     * state matrix to the associated state derivative functions.
     */
    typedef std::map< StateBlockIndices, PartStateDerivativeFunction > StateDerivativeModelMap;

    //! Typedef for state-derivative model map for vector composite state.
    /*!
     * Typedef for state-derivative model map, that maps part states (vectors) in the composite
     * state vector to the associated state derivative functions.
     */
    typedef std::map< StateSegmentIndices, PartStateDerivativeFunction >
    VectorStateDerivativeModelMap;

    //! Constructor taking a state-derivative model map (matrix) and an update function.
    /*!
     * Constructor taking a state-derivative model map, that maps part states (matrices) in the
     * composite state matrix with associated state derivative functions, and an update function
     * that updates the values of the independent variable and the composite state, as well
     * as any dependent variables not included in the state, in the data repository created
     * externally by the user.
     * \param aStateDerivativeModelMap A state derivative model map for matrix-based composite
                states.
     * \param anUpdateIndependentVariableAndStateFunction A function to update the independent
     *          variable and composite state in the user's data repository.
     */
    CompositeStateDerivativeModel( const StateDerivativeModelMap& aStateDerivativeModelMap,
                                   const IndependentVariableAndStateUpdateFunction
                                   anUpdateIndependentVariableAndStateFunction )
        : stateDerivativeModelMap( aStateDerivativeModelMap ),
          updateIndependentVariableAndState( anUpdateIndependentVariableAndStateFunction )
    { }

    //! Constructor taking a state-derivative model map (vector) and an update function.
    /*!
     * Constructor taking a state-derivative model map, that maps part states (vectors) in the
     * composite state vector with associated state derivative functions, and an update function
     * that updates the values of the independent variable and the composite state, as well
     * as any dependent variables not included in the state, in the data repository created
     * externally by the user.
     * \param aVectorStateDerivativeModelMap A state derivative model map for vector-based composite
                states.
     * \param anUpdateIndependentVariableAndStateFunction A function to update the independent
     *          variable and composite state in the user's data repository.
     */
    CompositeStateDerivativeModel( const VectorStateDerivativeModelMap&
                                   aVectorStateDerivativeModelMap,
                                   const IndependentVariableAndStateUpdateFunction
                                   anUpdateIndependentVariableAndStateFunction );

    //! Compute state derivative.
    /*!
     * Computes the state derivative based on the state-derivative model map provided through the
     * constructor. This function generates the composite state derivative, given a current value
     * of the independent variable and the composite state.
     * \param independentVariable Current independent variable value.
     * \param compositeState Current composite state.
     * \return Composite state derivative.
     */
    CompositeStateDerivativeType computeStateDerivative(
            const IndependentVariableType independentVariable,
            const CompositeStateType& compositeState );

protected:

private:

    //! State-derivative model map.
    /*!
     * State-derivative model map that maps part states in the composite state with associated
     * state derivative models.
     */
    StateDerivativeModelMap stateDerivativeModelMap;

    //! Pointer to update function.
    /*!
     * Pointer to a function that updates the user's data repository with the current independent
     * variable and state data to the current values.
     */
    const IndependentVariableAndStateUpdateFunction updateIndependentVariableAndState;
};

//! Constructor taking a state-derivative model map (vector) and an update function.
template< typename IndependentVariableType, typename CompositeStateType, typename PartStateType >
CompositeStateDerivativeModel< IndependentVariableType, CompositeStateType, PartStateType >::
CompositeStateDerivativeModel( const VectorStateDerivativeModelMap& aVectorStateDerivativeModelMap,
                               const IndependentVariableAndStateUpdateFunction
                               anUpdateIndependentVariableAndStateFunction )
    : updateIndependentVariableAndState( anUpdateIndependentVariableAndStateFunction )
{
    // Loop through vector state derivative model map and copy to matrix state derivative model
    // map.
    for ( typename VectorStateDerivativeModelMap::const_iterator
          iteratorVectorStateDerivativeModelMap = aVectorStateDerivativeModelMap.begin( );
          iteratorVectorStateDerivativeModelMap != aVectorStateDerivativeModelMap.end( );
          iteratorVectorStateDerivativeModelMap++ )
    {
        // Make tuple for input to .block function from provided input to .segment function.
        stateDerivativeModelMap[ boost::make_tuple(
                    iteratorVectorStateDerivativeModelMap->first.first, 0,
                    iteratorVectorStateDerivativeModelMap->first.second, 1 ) ]
                = iteratorVectorStateDerivativeModelMap->second;
    }
}

//! Compute state derivative.
template< typename IndependentVariableType, typename CompositeStateType, typename PartStateType >
CompositeStateType CompositeStateDerivativeModel< IndependentVariableType, CompositeStateType,
PartStateType >::computeStateDerivative( const IndependentVariableType independentVariable,
                                         const CompositeStateType& compositeState )
{
    // Update to current data.
    updateIndependentVariableAndState( independentVariable, compositeState );

    // Declare composite state derivative to be the same size as the composite state.
    CompositeStateType compositeStateDerivative
            = CompositeStateType::Zero( compositeState.rows( ), compositeState.cols( ) );

    // Loop through the state-derivative model map and compute the elements of the composite
    // state derivative.
    for ( typename StateDerivativeModelMap::const_iterator iteratorStateDerivativeModels
          = stateDerivativeModelMap.begin( );
          iteratorStateDerivativeModels != stateDerivativeModelMap.end( );
          iteratorStateDerivativeModels++ )
    {
        // Set block indices for part state derivative.
        unsigned int startRow = boost::get< 0 >( iteratorStateDerivativeModels->first );
        unsigned int startColumn = boost::get< 1 >( iteratorStateDerivativeModels->first );
        unsigned int numberOfRows = boost::get< 2 >( iteratorStateDerivativeModels->first );
        unsigned int numberOfColumns = boost::get< 3 >( iteratorStateDerivativeModels->first );

        // Compute the associated part in the composite state derivative based on the mapping
        // provided.
        compositeStateDerivative.block( startRow, startColumn, numberOfRows, numberOfColumns )
                = iteratorStateDerivativeModels->second( independentVariable,
                                                         compositeState.block(
                                                             startRow, startColumn,
                                                             numberOfRows, numberOfColumns ) );
    }

    // Return the composite state derivative computed.
    return compositeStateDerivative;
}

//! Typedef for a composite state derivative model with independent-variable-type = double,
//! CompositeState(Derivative)Type = Eigen::MatrixXd, PartState(Derivative)Type = Eigen::MatrixXd.
typedef CompositeStateDerivativeModel< > CompositeStateDerivativeModelMatrixXd;

//! Typedef for shared-pointer to a CompositeStateDerivativeModelMatrixXd object.
typedef boost::shared_ptr< CompositeStateDerivativeModelMatrixXd >
CompositeStateDerivativeModelMatrixXdPointer;

//! Typedef for a composite state derivative model with independent-variable-type = double,
//! CompositeState(Derivative)Type = Eigen::VectorXd, PartState(Derivative)Type = Eigen::VectorXd.
typedef CompositeStateDerivativeModel< double, Eigen::VectorXd >
CompositeStateDerivativeModelVectorXd;

//! Typedef for shared-pointer to a CompositeStateDerivativeModelVectorXd object.
typedef boost::shared_ptr< CompositeStateDerivativeModelVectorXd >
CompositeStateDerivativeModelVectorXdPointer;

} // namespace state_derivative_models
} // namespace tudat

#endif // TUDAT_COMPOSITE_STATE_DERIVATIVE_MODEL_H
