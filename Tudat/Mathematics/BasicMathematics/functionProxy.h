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

#ifndef TUDAT_FUNCTION_PROXY_H
#define TUDAT_FUNCTION_PROXY_H

#include <map>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"

namespace tudat
{
namespace basic_mathematics
{

//! Function interface to allow evaluation of a mathematical function.
/*!
 * Adds the ability for an already defined C++ function to be evaluated as a mathematical function.
 *
 * You need to set one C++ binding representing the mathematical function in the constructor.
 * Additional derivative and integral forms can be set using the addBinding function. If no 
 * explicit binding is found, numerical techniques are used to solve for derivatives and integrals.
 * \tparam IndependentVariable Mathematical variable for independent input (e.g., x).
 * \tparam DependentVariable Mathematical variable for dependent output (e.g., y=f(x)).
 */
template< typename IndependentVariable = double, typename DependentVariable = double >
class FunctionProxy : public BasicFunction< IndependentVariable, DependentVariable >
{
public:

    // Useful definitions.
    //! Typedef for a basic function.
    /*!
     * the independent and dependent variables are both doubles.
     */
    typedef BasicFunction< IndependentVariable, DependentVariable >           Parent;
    //! Typedef for the function that relates the dependent variable to the independent variable.
    /*!
     * the independent and dependent variables are both doubles.
     */
    typedef boost::function< DependentVariable( IndependentVariable ) >       FunctionSignature;
    //! Typedef for a shared pointer to the class FunctionProxy
    typedef boost::shared_ptr< FunctionProxy >                                FunctionProxyPointer;

    //! Create a Function object, using a specified function pointer.
    /*!
     * \param realFunction C++ function that will be used to evaluate this mathematical function.
     */
    FunctionProxy( FunctionSignature realFunction ) : realFunction_( realFunction ) { }

    //! Evaluate the mathematical function.
    /*!
     * Performs simple delegation of the evaluation to the user set function.
     * \param independentVariable Location where to evaluate the function.
     * \see Function::evaluate( IndependentVariable ).
     * \return the value of the dependent varaible.
     */
    virtual DependentVariable evaluate( const IndependentVariable independentVariable )
    {
        return realFunction_( independentVariable );
    }

    //! Evaluate the derivative of the function.
    /*!
     * This evaluates the derivative of a given order of this function.
     * \param order                   Order of the derivative to evaluate.
     * \param independentVariable     Location where to evaluate the derivative.
     * \return                        The function that gives the derivative.
     */
    DependentVariable computeDerivative( const unsigned int order,
                                         const IndependentVariable independentVariable )
    {
        // Find the corresponding function and evaluate it immediately at x.
        return findBinding( -static_cast< int >( order ) )( independentVariable );
    }

    //! Add a binding for an explicit form of a function derivative or integral.
    /*!
     * Using this method, you can bind additional C++ functions to either a specific derivative or
     * an integral of the base function.
     *
     * Don't add 0th order binding: it is bound separately (constructor) so no table lookup is
     * required.
     *
     * \param order Order of integration (negative for a derivative function).
     * \param function The functionsignature of the function set by the used.
     */
    void addBinding( const int order, FunctionSignature function )
    {
        functionCallTable_[ order ] = function;
    }
    
    //! Get a function representing the derivative or integral previously set.
    /*!
     * Retrieves a binding that was set using FunctionProxy::addBinding(int, FunctionSignature)
     * \param order Order of integration (negative for a derivative function).
     * \return the binding
     */
    FunctionSignature findBinding( int order = 0 )
    {
        // If 0th order, return the function itself.
        if ( order == 0 )
        {
            return realFunction_;
        }

        // Otherwise perform table lookup and return corresponding binding.
        return functionCallTable_.find( order )->second;
    }

protected:

private:

    //! Function representing the 0th order of this mathematical function.
    FunctionSignature realFunction_;

    //! Lookup table used to find explicit derivative or integral forms.
    std::map< int, FunctionSignature >  functionCallTable_;
};

// Type definitions for the commonly used double form.
typedef FunctionProxy< double, double >         UnivariateProxy;
typedef UnivariateProxy::FunctionSignature      UnivariateSignature;
typedef UnivariateProxy::FunctionProxyPointer   UnivariateProxyPointer;

//! Factory for creating a UnivariateProxyPointer from a C++ function UnivariateSignature.
/*!
 * \param function C++ function representing the mathematical function.
 */
inline UnivariateProxyPointer univariateProxy( UnivariateSignature function )
{
    // Return shared pointer to univariate function.
    return boost::make_shared< UnivariateProxy >( function );
}

//! Factory for creating a UnivaritateProxyPtr from a C++ function UnivariateSignature, with one
//! additional binding.
/*!
 * \param function  C++ function representing the mathematical function.
 * \param order1    Order of the additional derivative or integral function, i.e. order of
 *                  function1.
 * \param function1 Additional derivative or integral representation to bind.
 */
inline UnivariateProxyPointer univariateProxy( UnivariateSignature function,
                                               const int order1, UnivariateSignature function1 )
{
    // Create shared pointer to univariate function.
    UnivariateProxyPointer functionPointer = boost::make_shared< UnivariateProxy >( function );

    // Add a binding for an explicit form of the function derivative or integral.
    functionPointer->addBinding( order1, function1 );

    // Return function pointer.
    return functionPointer;
}

//! Factory for creating a UnivariateProxyPointer from a C++ function UnivariateSignature, with two
//! additional binding.
/*!
 * \param function  C++ function representing the mathematical function.
 * \param order1    Order of the first additional derivative or integral function, i.e. order of
 *                  function1.
 * \param function1 First additional derivative or integral representation to bind.
 * \param order2    Order of the second additional derivative or integral function, i.e. order of
 *                  function2.
 * \param function2 Second additional derivative or integral representation to bind.
 */
inline UnivariateProxyPointer univariateProxy( UnivariateSignature function,
                                               const int order1, UnivariateSignature function1,
                                               const int order2, UnivariateSignature function2 )
{
    // Create shared pointer to univariate function.
    UnivariateProxyPointer functionPointer = boost::make_shared< UnivariateProxy >( function );

    // Add bindings for the explicit form of the function derivative(s) or integral(s).
    functionPointer->addBinding( order1, function1 );
    functionPointer->addBinding( order2, function2 );

    // Return function pointer.
    return functionPointer;
}

//! Factory for creating a UnivariateProxyPointer from a C++ function UnivariateSignature, with three
//! additional binding.
/*!
 * \param function  C++ function representing the mathematical function.
 * \param order1    Order of the first additional derivative or integral function, i.e. order of
 *                  function1.
 * \param function1 First additional derivative or integral representation to bind.
 * \param order2    Order of the second additional derivative or integral function, i.e. order of
 *                  function2.
 * \param function2 Second additional derivative or integral representation to bind.
 * \param order3    Order of the third additional derivative or integral function, i.e. order of
 *                  function3.
 * \param function3 Third additional derivative or integral representation to bind.
 */
inline UnivariateProxyPointer univariateProxy( UnivariateSignature function,
                                               const int order1, UnivariateSignature function1,
                                               const int order2, UnivariateSignature function2,
                                               const int order3, UnivariateSignature function3 )
{
    // Create shared pointer to univariate function.
    UnivariateProxyPointer functionPointer = boost::make_shared< UnivariateProxy >( function );

    // Add bindings for the explicit form of the function derivative(s) or integral(s).
    functionPointer->addBinding( order1, function1 );
    functionPointer->addBinding( order2, function2 );
    functionPointer->addBinding( order3, function3 );

    // Return function pointer.
    return functionPointer;
}

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_FUNCTION_PROXY_H
