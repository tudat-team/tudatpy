/*! \file numericalPropagator.cpp
 *    Source file that defines the numerical propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 15 September, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The code only handles 1st-order Ordinary Differential Equations (ODEs)
 *      at the moment. In future, the code will be changed to accommodate
 *      2nd-order ODE's.
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
 *      100915    K. Kumar          File created.
 *      100928    K. Kumar          Modified propagate function; added fixed
 *                                  output interval option.
 *      100929    J. Melman         Corrected several alignments; added some
 *                                  comments.
 *      100929    K. Kumar          Minor comment modifications.
 *      110202    K. Kumar          Updated code to use Integrator adaptor
 *                                  instead of pointers-to-member functions;
 *                                  Update code to use State class.
 *      110203    J. Melman         Possibly encountered a mistake regarding
 *                                  the addition of the forces.
 *      110207    K. Kumar          Updated code with corrections.
 */

// Include statements.
#include "numericalPropagator.h"

// Using declarations.
using std::endl;

//! Default constructor.
NumericalPropagator::NumericalPropagator( )
{
}

//! Default destructor.
NumericalPropagator::~NumericalPropagator( )
{
}

//! Set integrator for propagation.
void NumericalPropagator::setIntegrator( Integrator* pointerToIntegrator )
{
    // Set pointer to object of Integrator class.
    pointerToIntegrator_ = pointerToIntegrator;
}

//! Add force model for propagation of a body.
void NumericalPropagator::addForceModel( Body* pointerToBody,
                                         ForceModel* pointerToForceModel )
{
    // Set force model in vector container for given body.
    bodiesToPropagate_[ pointerToBody ]
            .vectorContainerOfPointersToForceModels_
                     .push_back( pointerToForceModel );
}

// End of file.
