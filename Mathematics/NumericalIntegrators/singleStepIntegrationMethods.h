/*! \file singleStepIntegrationMethods.h
 *    Header file that defines the base class for all single step integration
 *    methods included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 27 July, 2010
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
 *      100907    K. Kumar          File header and footer added.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110203    J. Melman         File checked.
 *      110207    K. Kumar          Path changed; moved integrate() function;
 *                                  added ostream overload.
 */

#ifndef SINGLESTEP_H
#define SINGLESTEP_H

// Include statements.
#include <iostream>
#include "integrator.h"

//! Single step integration methods class.
/*!
 * Base class for all single step integration methods included in Tudat.
 */
class SingleStepIntegrationMethods : public Integrator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SingleStepIntegrationMethods( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SingleStepIntegrationMethods( );

    //! Integrate.
    /*!
     * This function executes an integration.
     */
    void integrate( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param pointerToIntegrator Pointer to Integrator.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     Integrator* pointerToIntegrator );

protected:

    //! Compute next state.
    /*!
     * Computes the next state with respect to the state at the start of an
     * integration step.
     * \param stepsize Stepsize.
     */
    virtual void computeNextState_( const double& stepsize ) = 0;

private:
};

#endif // SINGLESTEP_H

// End of file.
