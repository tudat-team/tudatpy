/*! \file convertMeanAnomalyBase.cpp
 *    This source file contains a base class to convert mean anomaly to
 *    eccentric and hyperbolic eccentric anomalies.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
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
 *    Date created      : 10 February, 2011
 *    Last modified     : 10 February, 2011
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
 *      110210    K. Kumar          First creation of code.
 */

// Include statements.
#include "convertMeanAnomalyBase.h"

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Default constructor.
ConvertMeanAnomalyBase::ConvertMeanAnomalyBase( )
    : eccentricity_( -1.0 ),
      meanAnomaly_( -1.0 ),
      pointerToNewtonRaphson_( NULL )
{
}

//! Default destructor.
ConvertMeanAnomalyBase::~ConvertMeanAnomalyBase( )
{
}

//! Set Newton-Raphson method.
void ConvertMeanAnomalyBase::
        setNewtonRaphson( NewtonRaphson* pointerToNewtonRaphson )
{
    // Set pointer to Newton-Raphson method.
    pointerToNewtonRaphson_ = pointerToNewtonRaphson;
}

//! Set eccentricity.
void ConvertMeanAnomalyBase::setEccentricity( const double& eccentricity )
{
    eccentricity_ = eccentricity;
}

//! Set mean anomaly.
void ConvertMeanAnomalyBase::setMeanAnomaly( const double& meanAnomaly )
{
    meanAnomaly_ = meanAnomaly;
}

}

// End of file.
