/*! \file euler.h
 *    Header file that defines the Euler integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 27 July, 2010
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
 *      100907    K. Kumar          File header and footer added.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         File checked.
 *      110207    K. Kumar          Path changed; moved integrat() function
 *                                  to SingleStepIntegrationMethods.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef EULER_H
#define EULER_H

// Include statements.
#include "singleStepIntegrationMethods.h"

//! Euler integrator class.
/*!
 * Implementation of Euler integrator.
 */
class Euler : public SingleStepIntegrationMethods
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    Euler( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~Euler( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param eulerIntegrator Euler integrator.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     Euler& eulerIntegrator );

protected:

private:

    //! Compute next state.
    /*!
     * Computes the next state with respect to the state at the start of an
     * integration step.
     * \param stepsize Stepsize.
     */
    void computeNextState_( const double& stepsize );
};

#endif // EULER_H

// End of file.
