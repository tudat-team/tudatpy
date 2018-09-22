/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_FUNCTION_H
#define TUDAT_FUNCTION_H

#include <memory>

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
typedef std::shared_ptr< Function< > > FunctionPointer;

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_FUNCTION_H
