/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120222    B. Tong Minh      File created.
 *      121120    K. Kumar          Renamed file and removed code that belongs in derived classes.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H
#define TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Base class for reinitializable numerical integrators.
/*!
 * This base class enables derived classes to reinitialize the current state of the numerical
 * integrator.
 * \tparam IndependentVariableType The type of the independent variable.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType.
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd >
class ReinitializableNumericalIntegrator :
        public numerical_integrators::NumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType >
{
protected:

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::NumericalIntegrator<
    IndependentVariableType, StateType, StateDerivativeType > NumericalIntegratorBase;

public:

    //! Default constructor.
    ReinitializableNumericalIntegrator(
            const typename NumericalIntegratorBase::StateDerivativeFunction&
            aStateDerivativeFunction )
        : NumericalIntegratorBase( aStateDerivativeFunction )
    { }

    //! Default destructor.
    virtual ~ReinitializableNumericalIntegrator( ) { }

    //! Modify the state at the current value of the independent variable.
    /*!
     * Modify the state at the current value of the independent variable. This function is pure
     * virtual; hence it must be implemented in all derived classes.
     * \param newState The new state to set the current state to.
     */
    virtual void modifyCurrentState( const StateType& newState ) = 0;

protected:

private:
};

//! Typedef for shared-pointer to default, re-initializable numerical integrator.
/*!
 * Typedef for shared-pointer to a default, re-initializable numerical integrator
 * (IndependentVariableType = double, StateType = Eigen::VectorXd,
 * StateDerivativeType = Eigen::VectorXd).
 */
typedef boost::shared_ptr< ReinitializableNumericalIntegrator< > >
ReinitializableNumericalIntegratorXdPointer;

//! Typedef for a shared-pointer to a scalar, re-initializable numerical integrator.
/*!
 * Typedef for shared-pointer to a scalar numerical, re-initializable integrator
 * (IndependentVariableType = double, StateType = double, StateDerivativeType = double).
 */
typedef boost::shared_ptr< ReinitializableNumericalIntegrator< double, double, double > >
ReinitializableNumericalIntegratordPointer;


} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_REINITIALIZABLE_NUMERICAL_INTEGRATOR_H
