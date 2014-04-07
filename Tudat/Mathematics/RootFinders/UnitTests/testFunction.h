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
 *      120318    S. Billemont      File created.
 *      120402    T. Secretin       Code-check.
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
