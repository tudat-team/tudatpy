/*! \file integratorAdaptor.h
 *    This header file contains an adaptor template class for the Integrator
 *    class included in Tudat. This class is used to avoid the need for
 *    specific pointers-to-member-functions.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 1 February, 2011
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110201    K. Kumar          First creation of code.
 *      110203    J. Melman         File checked.
 *      110207    K. Kumar          Path changed.
 */

#ifndef INTEGRATORADAPTOR_H
#define INTEGRATORADAPTOR_H

// Include statements.
#include "integratorBase.h"
#include "state.h"

//! Template class for the IntegratorAdaptor.
/*!
 * This is a template class for an adaptor used by the Integrator class.
 */
template < class TClass >
class IntegratorAdaptor : public IntegratorBase
{
public:

    // Definition of typedef.
    // State derivative function used by the Integrator class can be passed as
    // pointers-to-member-functions by making use of this adaptor class.
    typedef State* ( TClass::*pointerToTClassMemberFunction )( State* );

    //! Default constructor.
    /*!
     * Default constructor.
     */
    IntegratorAdaptor( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~IntegratorAdaptor( ){ }

    //! Compute state derivative.
    /*!
     * Computes the state derivative for the numerical integrator being used.
     * The state derivative function is a member function of the template class
     * TClass.
     * \param pointerToState State given as pointer to State object.
     * \return State derivative given as pointer to State object.
     */
    virtual State* computeStateDerivative( State* pointerToState );

    //! Set class that contains state derivative function.
    /*!
     * Sets the class that contains the state derivative function used for the
     * Integrator class. This function must be defined within the target
     * class and must be of the form: State* functionName( State* ).
     * \param pointerToClass Pointer to class.
     */
    void setClass( TClass* pointerToClass );

    //! Set pointer to state derivative function in target class.
    /*!
     * Sets a pointer to the state derivative function contained in the target
     * class, used for the Integrator class. This member function must be
     * defined within the target class and must be of the form:
     * State* functionName( State* ).
     * \param stateDerivativefunction Pointer to state derivative (member)
     *          function.
     */
    void setStateDerivativeFunction( pointerToTClassMemberFunction
                                     stateDerivativefunction );

protected:

private:

    //! Pointer to target class.
    /*!
     * Pointer to target class containing the state derivative function
     * necessary for the Integrator class.
     */
    TClass* pointerToTClass_;

    //! Pointer to state derivative (member) function.
    /*!
     * Pointer to state derivative function contained within target class, used
     * for the Integrator class. This function must be of the form:
     * State* functionName( State* ).
     */
    pointerToTClassMemberFunction pointerToStateDerivativeFunction_;
};

//! Compute state derivative.
template < class TClass >
State* IntegratorAdaptor< TClass >::computeStateDerivative(
        State* pointerToState )
{
    return ( pointerToTClass_
             ->*pointerToStateDerivativeFunction_ )( pointerToState );
}

//! Set class that contains state derivative function.
template < class TClass >
void IntegratorAdaptor< TClass >::setClass( TClass* pointerToClass )
{
    pointerToTClass_ = pointerToClass;
}

//! Set pointer to state derivative function in target class.
template < class TClass >
void IntegratorAdaptor< TClass >::
        setStateDerivativeFunction( pointerToTClassMemberFunction
                                    stateDerivativefunction )
{
    pointerToStateDerivativeFunction_ = stateDerivativefunction;
}

#endif // INTEGRATORADAPTOR_H

// End of file.
