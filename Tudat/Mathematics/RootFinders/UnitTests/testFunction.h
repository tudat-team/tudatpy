/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Milea, A. Virtual Inheritance in C++, and solving the diamond problem,
 *          http://www.cprogramming.com/tutorial/virtual_inheritance.html, 1997-2011,
 *          last accessed: 9th September, 2012.
 *
 *    Notes
 *      You need to define what function implementation you use yourself to avoid the Diamond
 *      Problem (Milea, 2011).
 *
 */

#ifndef TUDAT_TEST_FUNCTION_H
#define TUDAT_TEST_FUNCTION_H

#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"

namespace tudat
{
namespace unit_tests
{

//! Simple definition of a test function, so that it can be used by all root-finder unit tests.
struct TestFunction
{
    //! Default destructor.
    virtual ~TestFunction( ) { }

    //! Expected true location of the root.
    virtual double getTrueRootLocation( ) = 0;
    
    //! Accuracy of the true value of the root.
    virtual double getTrueRootAccuracy( ) = 0;

    //! Get a reasonable initial guess of the root location.
    virtual double getInitialGuess( ) = 0;

    //! Get a reasonable lower boundary for the root location.
    virtual double getLowerBound( ) = 0;

    //! Get a reasonable upper boundary for the root location.
    virtual double getUpperBound( ) = 0;
};

} // namespace unit_tests
} // tudat

#endif // TUDAT_TEST_FUNCTION_H
