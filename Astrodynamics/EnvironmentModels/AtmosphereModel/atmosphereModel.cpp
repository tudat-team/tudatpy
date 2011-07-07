/*! \file atmosphereModel.cpp
*    Source file that defines the baseclass atmosphere model included in Tudat.
*
*    Path              : /Astrodynamics/EnvironmentModels/AtmosphereModel/
*    Version           : 3
*    Check status      : Unchecked
*
*    Author            : F.M. Engelen
*    Affiliation       : Delft University of Technology
*    E-mail address    : F.M.Engelen@student.tudelft.nl
*
*    Checker           : K. Kumar
*    Affiliation       : Delft University of Technology
*    E-mail address    : K.Kumar@tudelft.nl
*
*    Checker           : J. Melman
*    Affiliation       : Delft University of Technology
*    E-mail address    : J.C.P.Melman@tudelft.nl
*
*    Date created      : 10 March, 2011
*    Last modified     : 27 April, 2011
*
*    References
*
*    Notes
*       Base class for standard and reference atmospheres.
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
*      110310    F.M. Engelen      File created.
*      110324    J. Melman         Set time to zero by default.
*      110412    F.M. Engelen      Structure to delete default altitude model.
*      110427    F.M. Engelen      Changed input arguments to altitude, longitude, latitude.
*/

// Include statements.
#include "atmosphereModel.h"

//! Default constructor.
AtmosphereModel::AtmosphereModel( )
{
}

//! Default destructor.
AtmosphereModel::~AtmosphereModel( )
{
}

// End of file.
