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
 *      120823    K. Kumar          File created.
 *      120824    K. Kumar          Created from original StateDerivativeModel class.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_STATE_DERIVATIVE_MODEL_H
#define TUDAT_STATE_DERIVATIVE_MODEL_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace state_derivative_models
{

//! State derivative model class.
/*!
 * This class contains a state derivative model base class.
 * \tparam StateDerivativeType Data type of state derivative (default=Eigen::MatrixXd).
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::MatrixXd >
class StateDerivativeModel
{
private:

    //! Typedef for state derivative type.
    typedef StateType StateDerivativeType;

public:

    //! Virtual destructor.
    /*!
     * Virtual destructor that ensures that derived class destructors are called correctly.
     */
    virtual ~StateDerivativeModel( ) { }

    //! Compute state derivative.
    /*!
     * Computes state derivative. This is a pure virtual function, which should be implemented
     * in derived classes.
     * \return Computed state derivative.
     */
    virtual StateDerivativeType computeStateDerivative(
            const IndependentVariableType independentVariable, const StateType& state ) = 0;

protected:

private:
};

//! Typedef for shared-pointer to state derivative model (dynamically sized matrix).
typedef boost::shared_ptr< StateDerivativeModel< > > StateDerivativeModelXd;

//! Typedef for state derivative model (Vector6d).
typedef StateDerivativeModel< double, basic_mathematics::Vector6d > StateDerivativeModelVector6d;

//! Typedef for shared-pointer to state derivative model (Vector6d).
typedef boost::shared_ptr< StateDerivativeModelVector6d > StateDerivativeModelVector6dPointer;

} // namespace state_derivative_models
} // namespace tudat

#endif // TUDAT_STATE_DERIVATIVE_MODEL_H
