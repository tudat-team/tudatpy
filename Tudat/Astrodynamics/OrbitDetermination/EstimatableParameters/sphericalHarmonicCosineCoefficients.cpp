/*    Copyright (c) 2010-2019, Delft University of Technology
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

    for( unsigned int i = 0; i < blockIndices_.size( ); i++ )
    {
        coefficients( blockIndices_.at( i ).first, blockIndices_.at( i ).second ) = parameterValue( i );
    }
    setCosineCoefficients_( coefficients );
}

//! Function to get a list of Kaula constraint values for gravity field coefficients at given degrees and indices
Eigen::VectorXd getKaulaConstraintVector(
        const std::vector< std::pair< int, int > > blockIndices,
        const double constraintMultiplier )
{
    Eigen::VectorXd constraints = Eigen::VectorXd::Zero( blockIndices.size( ) );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        constraints( i ) = constraintMultiplier / ( blockIndices.at( i ).first * blockIndices.at( i ).first );
    }
    return constraints;
}

//! Function to get a list of Kaula constraint values for gravity field coefficients for given parameter
Eigen::VectorXd getKaulaConstraintVector(
        const std::shared_ptr< SphericalHarmonicsCosineCoefficients > parameter,
        const double constraintMultiplier )
{
    return getKaulaConstraintVector( parameter->getBlockIndices( ), constraintMultiplier );
}


}

}


