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
 *      110119    S. Billemont      File created.
 *      120813    P. Musegaas       Added some forgotten templates.
 *      120909    K. Kumar          Added Doxygen comments; updated variable-naming to be verbose;
 *                                  updated function names and ensured const-correctness.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      140124    H.P. Gijsen       fixed Doxygen warnings
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_FUNCTION_H
#define TUDAT_FUNCTION_H

#include <boost/shared_ptr.hpp>

namespace tudat
{
namespace basic_mathematics
{

//! Function interface to allow evaluation of a mathematical function.
/*!
 * This class adds the ability to a class to be evaluated as a mathematical function.
 * \tparam IndependentVariable Mathematical variable for independent input (e.g., x).
 * \tparam DependentVariable Mathematical variable for dependent output (e.g., y=f(x)).
 */
template< typename IndependentVariable = double, typename DependentVariable = double >
class Function
{
public:

    //! Default destructor.
    virtual ~Function( ) { }

    //! Compute mathematical function value. This implementation is a pure virtual function.
    /*!
     * Computes the value of the mathematical function used for the root-finder algorithm being
     * used.
     * \param inputValue Input value.
     * \return Computed mathematical function value.
     */
    virtual DependentVariable evaluate( const IndependentVariable inputValue ) = 0;
    
    //! Alias for evaluate( double ).
    inline DependentVariable operator( )( const double inputValue )
    {
        return evaluate( inputValue );
    }

    //! Evaluate the derivative of the function. This implementation is a pure virtual function.
    /*!
     * This evaluates the derivative of a given order of this function.
     * \param order               Order of the derivative to evaluate.
     * \param independentVariable Location where to evaluate the derivative.
     * \return the derivative of the function
     */
    virtual DependentVariable computeDerivative(
            const unsigned int order, const IndependentVariable independentVariable ) = 0;

    //! Evaluate the definite integral of the function.
    /*!
     * This evaluates the definite integral of a given order of this function. This implementation is a pure virtual function.
     * \param order      Order of the integral to evaluate.
     * \param lowerBound Integration lower bound (integrate from this point).
     * \param upperbound Integration upper bound (integrate to this point).
     * \return the definite numerical integral of this function
     */
    virtual DependentVariable computeDefiniteIntegral( const unsigned int order,
                                                       const IndependentVariable lowerBound,
                                                       const IndependentVariable upperbound ) = 0;
};

//! Typedef for shared-pointer to Function object.
/*!
 * Typedef for shared-pointer to Function object with IndependentVariable=double,
 * DependentVariable=double.
 */
typedef boost::shared_ptr< Function< > > FunctionPointer;

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_FUNCTION_H
