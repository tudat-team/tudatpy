/*! \file euler.h
 *    Header file that defines the Euler integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dikx@student.tudelft.nl
 *
 *    Date created      : 27 July, 2010
 *    Last modified     : 29 September, 2010
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
 *      YYMMDD    author              comment
 *      100907    K. Kumar            File header and footer added.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor comment modifications.
 */

#ifndef EULER_H
#define EULER_H

// Include statements.
#include "singleStepIntegrationMethods.h"

//! Euler Integrator class.
/*!
 * Implementation of Euler integrator in Tudat.
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

    //! Perform Euler integration.
    /*!
     * This function performs Euler integration.
     */
    void integrate( );

protected:

private:

    //! Compute next state vector.
    /*!
     * This function computes the next state vector with respect to the state
     * vector at the start of an integration step.
     * \param stepsize Stepsize
     */
    void computeNextStateVector( const double& stepsize );
};

#endif // EULER_H

// End of file.
