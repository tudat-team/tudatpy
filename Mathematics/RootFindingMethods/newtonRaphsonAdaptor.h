/*! \file newtonRaphsonAdaptor.h
 *    This header file contains an adaptor template class for the NewtonRaphson
 *    class included in Tudat. This class is used to avoid the need for
 *    pointers-to-member-functions.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 19 January, 2011
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
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
 *      YYMMDD    Author            Comment
 *      110119    K. Kumar          First creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef NEWTONRAPHSONADAPTOR_H
#define NEWTONRAPHSONADAPTOR_H

// Include statements.
#include "Mathematics/RootFindingMethods/newtonRaphsonBase.h"

//! Template class for the NewtonRaphsonAdaptor.
/*!
 * This is a template class for an adaptor used by the NewtonRaphson class.
 */
template < class TClass >
class NewtonRaphsonAdaptor : public NewtonRaphsonBase
{
public:

    //! Definition of typedef.
    /*!
     * Functions to which the Newton-Raphson method is applied can be passed as
     * pointers-to-member-functions by making use of this adaptor class.
     */
    typedef double ( TClass::*pointerToTClassMemberFunction )( double& );

    //! Default constructor.
    /*!
     * Default constructor.
     */
    NewtonRaphsonAdaptor( ){ };

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~NewtonRaphsonAdaptor( ){ };

    //! Compute mathematical function.
    /*!
     * Computes the value of the mathematical function being used for the
     * Newton-Raphson method. The mathematical function is a member function
     * of the template class TClass.
     * \param inputValue Input value.
     */
    double computeFunction( double& inputValue );

    //! Compute first-derivative mathematical function.
    /*!
     * Computes the value of the first-derivative mathematical function being
     * used for the Newton-Raphson method. The first-derivative mathematical
     * function is a member function of the template class TClass.
     * \param inputValue Input value.
     */
    double computeFirstDerivativeFunction( double& inputValue );

    //! Set class that contains mathematical functions.
    /*!
     * Sets the class that contains the mathematical function and corresponding
     * first-derivative mathematical function used for the Newton-Raphson
     * method. These functions must be defined within the target class and must
     * be of the form: double functionName( double& )
     * \param pointerToClass Pointer to class.
     */
    void setClass( TClass* pointerToClass );

    //! Set pointer to mathematical function in target class.
    /*!
     * Sets a pointer to the mathematical function contained in the target
     * class, used for the Newton-Raphson method. This member function must be
     * defined within the target class and must be of the form:
     * double functionName( double& ).
     * \param function Pointer to mathematical ( member ) function.
     */
    void setPointerToFunction( pointerToTClassMemberFunction function );

    //! Set pointer to first-derivative mathematical function in target class.
    /*!
     * Sets a pointer to the first-derivative mathematical function contained
     * in the target class, used for the Newton-Raphson method. This member
     * function must be defined within the target class and must be of the
     * form: double functionName( double& ).
     * \param firstDerivativeFunction Pointer to first-derivative mathematical
     *          ( member ) function.
     */
    void setPointerToFirstDerivativeFunction( pointerToTClassMemberFunction
                                              firstDerivativeFunction );

protected:

private:

    //! Pointer to target class.
    /*!
     * Pointer to target class containing mathematical functions
     * necessary for the Newton-Raphson method.
     */
    TClass* pointerToTClass_;

    //! Pointer to mathematical ( member ) function.
    /*!
     * Pointer to mathematical function contained within target class, used
     * for the Newton-Rapshon method. This function must be of the form:
     * double functionName( double& ).
     */
    pointerToTClassMemberFunction pointerToFunction_;

    //! Pointer to first-derivative mathematical ( member ) function.
    /*!
     * Pointer to first-derivative mathematical function contained within
     * target class, used for the Newton-Rapshon method. This function must be
     *  of the form: double functionName( double& ).
     */
    pointerToTClassMemberFunction pointerToFirstDerivativeFunction_;
};

//! Compute mathematical function.
template < class TClass >
double NewtonRaphsonAdaptor< TClass >::computeFunction( double& inputValue )
{
    return ( pointerToTClass_->*pointerToFunction_ )( inputValue );
}

//! Compute first-derivative mathematical function.
template < class TClass >
double NewtonRaphsonAdaptor< TClass >::
        computeFirstDerivativeFunction( double& inputValue )
{
    return ( pointerToTClass_
             ->*pointerToFirstDerivativeFunction_ )( inputValue );
}

//! Set class that contains mathematical functions.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::setClass( TClass* pointerToClass )
{
    pointerToTClass_ = pointerToClass;
}

//! Set pointer to mathematical function contained in target class.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::
        setPointerToFunction( pointerToTClassMemberFunction function )
{
    pointerToFunction_ = function;
}

//! Set pointer to first-derivative mathematical function in target class.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::
        setPointerToFirstDerivativeFunction( pointerToTClassMemberFunction
                                             firstDerivativeFunction )
{
    pointerToFirstDerivativeFunction_ = firstDerivativeFunction;
}

#endif // NEWTONRAPHSONADAPTOR_H

// End of file.
