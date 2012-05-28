/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110119    K. Kumar          Creation of code.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *
 *    References
 *
 */

#ifndef TUDAT_NEWTON_RAPHSON_ADAPTOR_H
#define TUDAT_NEWTON_RAPHSON_ADAPTOR_H

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonBase.h"

namespace tudat
{

//! Template class for the NewtonRaphsonAdaptor.
/*!
 * A template class for an adaptor used by the NewtonRaphson class.
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
    NewtonRaphsonAdaptor( ) : pointerToTClass_( NULL ), pointerToFunction_( NULL ),
        pointerToFirstDerivativeFunction_( NULL ) { }

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
double NewtonRaphsonAdaptor< TClass >::computeFirstDerivativeFunction( double& inputValue )
{
    return ( pointerToTClass_->*pointerToFirstDerivativeFunction_ )( inputValue );
}

//! Set class that contains mathematical functions.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::setClass( TClass* pointerToClass )
{
    pointerToTClass_ = pointerToClass;
}

//! Set pointer to mathematical function contained in target class.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::setPointerToFunction( pointerToTClassMemberFunction function )
{
    pointerToFunction_ = function;
}

//! Set pointer to first-derivative mathematical function in target class.
template < class TClass >
void NewtonRaphsonAdaptor< TClass >::setPointerToFirstDerivativeFunction(
    pointerToTClassMemberFunction firstDerivativeFunction )
{
    pointerToFirstDerivativeFunction_ = firstDerivativeFunction;
}

} // namespace tudat

#endif // TUDAT_NEWTON_RAPHSON_ADAPTOR_H
