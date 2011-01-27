/*! \file sphericalHarmonicsGravityField.cpp
 *    Source file that defines the spherical harmonics gravity field model
 *    included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModel/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 17 November, 2010
 *    Last modified     : 6 January, 2011
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
 *      101117    K. Kumar            File created.
 *      101214    K. Kumar            Updated getGradientOfPotential() and
 *                                    getLaplacianOfPotential().
 *      101215    K. Kumar            Simplified getGradientOfPotential() and
 *                                    getLaplacianOfPotential().
 *      101216    K. Kumar            Updated functions to use position vector
 *                                    of origin for relative position.
 *      110106    K. Kumar            Added set/get functions for degree and
 *                                    order of expansion.
 */

// Include statements.
#include "sphericalHarmonicsGravityField.h"

// Using directives.
using mathematics::raiseToIntegerPower;

//! Default constructor.
SphericalHarmonicsGravityField::SphericalHarmonicsGravityField( )
{
    // Initialize variables.
    referenceRadius_ = -0.0;
    degreeOfExpansion_ = -0;
    orderOfExpansion_ = -0;
}

//! Default destructor.
SphericalHarmonicsGravityField::~SphericalHarmonicsGravityField( )
{
}

//! Set the reference radius.
void SphericalHarmonicsGravityField::setReferenceRadius(
        const double& referenceRadius )
{
    referenceRadius_ = referenceRadius;
}

//! Set degree of spherical harmonics gravity field expansion.
void SphericalHarmonicsGravityField::setDegreeOfExpansion(
        const int& degreeOfExpansion )
{
    degreeOfExpansion_ = degreeOfExpansion;
}

//! Set order of spherical harmonics gravity field expansion.
void SphericalHarmonicsGravityField::setOrderOfExpansion(
        const int& orderOfExpansion )
{
    orderOfExpansion_ = orderOfExpansion;
}

//! Get the reference radius.
double SphericalHarmonicsGravityField::getReferenceRadius( )
{
    return referenceRadius_;
}

//! Get degree of spherical harmonics gravity field expansion.
double SphericalHarmonicsGravityField::getDegreeOfExpansion( )
{
    return degreeOfExpansion_;
}

//! Get order of spherical harmonics gravity field expansion.
double SphericalHarmonicsGravityField::getOrderOfExpansion( )
{
    return orderOfExpansion_;
}

//! Get the gravitational potential.
double SphericalHarmonicsGravityField::getPotential(
        const Vector3d& positionVector )
{
    // Compute relative position vector.
    relativePositionVector_ = positionVector - positionVectorOfOrigin_;

    // Return the gravitational potential for a point mass.
    return gravitationalParameter_ / relativePositionVector_.norm();
}

//! Get the gradient of the gravitational potential.
Vector3d SphericalHarmonicsGravityField::getGradientOfPotential(
        const Vector3d& positionVector )
{
    // Compute relative position vector.
    relativePositionVector_ = positionVector - positionVectorOfOrigin_;

    // Compute and return gradient of potential.
    return -gravitationalParameter_
            * relativePositionVector_
            / raiseToIntegerPower(
                    relativePositionVector_.norm(), 3 );
}

//! Get the Laplacian of the gravitational potential.
Matrix3d SphericalHarmonicsGravityField::
        getLaplacianOfPotential( const Vector3d& positionVector )
{
    // Declare local variables.
    // Declare identity matrix for computations.
    Matrix3d identityMatrix_;

    // Compute relative position vector.
    relativePositionVector_ = positionVector - positionVectorOfOrigin_;

    // Set identity matrix to square size of state vector.
    identityMatrix_.setIdentity( 3, 3 );

    // Compute and return Laplacian of potential.
    return gravitationalParameter_
            / raiseToIntegerPower(
                    relativePositionVector_.norm( ), 5 )
            * ( ( 3.0 * relativePositionVector_
                  * relativePositionVector_.transpose( ) )
                - ( relativePositionVector_.squaredNorm( )
                    * identityMatrix_ ) );
}

// End of file.
