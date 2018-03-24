/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Function to retrieve the current values of the cosine coefficients that are to be estimated.
Eigen::VectorXd SphericalHarmonicsCosineCoefficients::getParameterValue( )
{
    Eigen::VectorXd parameterVector = Eigen::VectorXd::Zero( parameterSize_ );

    Eigen::MatrixXd coefficientBlock = getCosineCoefficients_( );

    for(  unsigned int i = 0; i < blockIndices_.size( ); i++ )
    {
        parameterVector( i ) = coefficientBlock( blockIndices_.at( i ).first, blockIndices_.at( i ).second );
    }
    return parameterVector;

}

//! Function to reset the cosine coefficients that are to be estimated.
void SphericalHarmonicsCosineCoefficients::setParameterValue( const Eigen::VectorXd parameterValue )
{
    Eigen::MatrixXd coefficients = getCosineCoefficients_( );

    for(  unsigned int i = 0; i < blockIndices_.size( ); i++ )
    {
        coefficients( blockIndices_.at( i ).first, blockIndices_.at( i ).second ) = parameterValue( i );
    }
    setCosineCoefficients_( coefficients );
}

}

}


