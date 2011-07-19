/*! \file aerodynamicCoefficientInterface.cpp
 *    This file contains the source definition of the AerodynamicCoefficientInterface
 *    base class included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 2
 *    Check status      : not checked
 *
 *    Author           :  F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@tudelft.nl
 *
 *    Date created      : 08 June 2011
 *    Last modified     : 14 July 2011
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
 *      110608    F.M. Engelen      First creation of code.
 *      110714    D. Dirkx          Class name change and other minor changes during code check.

 */

#include "aerodynamicCoefficientInterface.h"

//! Default constructor.
AerodynamicCoefficientInterface::AerodynamicCoefficientInterface( )
{
}

//! Default destructor.
AerodynamicCoefficientInterface::~AerodynamicCoefficientInterface( )
{
}

//! Set reference area.
void AerodynamicCoefficientInterface::setReferenceArea(
        const double& referenceArea )
{
    referenceArea_ = referenceArea;
}

//! Gets reference area.
double AerodynamicCoefficientInterface::getReferenceArea( )
{
    return referenceArea_ ;
}

//! Sets reference length.
void AerodynamicCoefficientInterface::setReferenceLength(
        const double& referenceLength )
{
    referenceLength_ = referenceLength;
}

//! Gets reference length.
double AerodynamicCoefficientInterface::getReferenceLength( )
{
    return referenceLength_ ;
}

//! Sets lateral reference length.
void AerodynamicCoefficientInterface::setLateralReferenceLength(
        const double& lateralReferenceLength )
{
    lateralReferenceLength_ = lateralReferenceLength;
}

//! Gets lateral reference length.
double AerodynamicCoefficientInterface::getLateralReferenceLength( )
{
    return lateralReferenceLength_;
}

//! Sets moment reference point.
void AerodynamicCoefficientInterface::setMomentReferencePoint(
        const Vector3d& momentReferencePoint)
{
    momentReferencePoint_ = momentReferencePoint;
}

//! Gets moment reference point.
VectorXd AerodynamicCoefficientInterface::getMomentReferencePoint( )
{
    return momentReferencePoint_ ;
}

//! Retreive the current force coefficients.
Vector3d AerodynamicCoefficientInterface::getCurrentForceCoefficients( )
{
    return currentForceCoefficients_;
}

//! Retreive the current moment coefficients.
Vector3d AerodynamicCoefficientInterface::getCurrentMomentCoefficients( )
{
    return currentMomentCoefficients_;
}

//! Function to set the force coefficients.
void AerodynamicCoefficientInterface::setCurrentForceCoefficients(
        const Vector3d& currentForceCoefficients )
{
    currentForceCoefficients_ = currentForceCoefficients;
}

//! Function to set the moment coefficients.
 void AerodynamicCoefficientInterface::setCurrentMomentCoefficients(
         const Vector3d& currentMomentCoefficients )
 {
     currentMomentCoefficients_ = currentMomentCoefficients;
 }

  // End of file.
