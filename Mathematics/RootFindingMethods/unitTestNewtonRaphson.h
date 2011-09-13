/*! \file unitTestNewtonRaphson.h
 *    Header file of unit test file of Newton-Raphson method code.
 *    This unit test file will test the Newton-Raphson method code.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 11 January, 2011
 *    Last modified     : 20 January, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author              Comment
 *      110111    E. Iorfida          First creation of the code.
 *      110119    K. Kumar            Updated code to work with adaptor and
 *                                    abstract base implementation so that
 *                                    pointer-to-member functions are not
 *                                    required; filename changed.
 *      110120    E. Iorfida          Added necessary class that contains functions,
 *                                    to allow a directly test with adaptor.
 *      110120    K. Kumar            Added global mathematical function test.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef UNITTESTNEWTONRAPHSONMETHOD_H
#define UNITTESTNEWTONRAPHSONMETHOD_H

// Include statements.
#include <cmath>
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Class for NewtonRaphson unit test code.
/*!
 * This class contains functions, necessary to test NewtonRaphson method.
 */
class NewtonRaphsonTest
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NewtonRaphsonTest( )
    {
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~NewtonRaphsonTest( )
    {
    }

    //! Mathematical test function.
    /*!
     * Mathematical test function used by the Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeTestFunction( double& inputValue );

    //! First-derivative of mathematical test function.
    /*!
     * First-derivative of mathematical test function used by the
     * Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeFirstDerivativeTestFunction( double& inputValue );

protected:

private:
};

//! Global mathematical test function.
/*!
 * Global mathematical test function used by the Newton-Raphson algorithm.
 * \param inputValue Input value.
 */
double computeGlobalTestFunction( double& inputValue );

//! Global first-derivative mathematical test function.
/*!
 * Global first-derivative mathematical test function used by the
 * Newton-Raphson algorithm.
 * \param inputValue Input value.
 */
double computeGlobalFirstDerivativeTestFunction( double& inputValue );

//! Test of Newton-Raphson method code.
/*!
 * Test of Newton-Raphson method code.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testNewtonRaphsonMethod( );

}

#endif // UNITTESTNEWTONRAPHSONMETHOD_H

// End of file.
