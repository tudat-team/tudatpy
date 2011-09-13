/*! \file state.cpp
 *    This source file contains a base class for all state classes in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 26 October, 2010
 *    Last modified     : 4 February, 2011
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
 *      101022    K. Kumar          First creation of code.
 *      101202    J. Melman         Changed path and dot added.
 *      110204    K. Kumar          Changed path; added missing comment.
 *      110207    K. Kumar          Added ostream overload.
 */

// Include statements.
#include "state.h"

// Using declarations.
using std::endl;

//! Default constructor.
State::State( )
{
}

//! Default destructor.
State::~State( )
{
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, State& stateObject )
{
    stream << "The state is set to: " << endl;
    stream << stateObject.state << endl;

    // Return stream.
    return stream;
}

// End of file.
