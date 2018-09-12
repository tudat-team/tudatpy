/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Schaub, H. and Junkins, J., Analytical Mechanics of Space Systems, 2nd ed., ser. Education Series.
 *          American Institute of Aeronautics and Astronautics, 2002.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/attitudeElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert quaternions to modified Rodrigues parameters.
Eigen::Vector4d convertQuaternionsToModifiedRodriguesParameterElements( const Eigen::Vector4d& quaternionElements )
{
    // Declare eventual output vector
    Eigen::Vector4d convertedModifiedRodriguesParameterElements;

    // Extract quaternions
    double etaQuaternionParameter = quaternionElements( etaQuaternionIndex );

    // Convert quaternions to modified Rodrigues parameters (or SMPR)
    bool shadowFlag = etaQuaternionParameter < 0;
    double conversionSign = shadowFlag ? - 1.0 : 1.0; // conversion is slightly different for SMRP and MRP
    convertedModifiedRodriguesParameterElements.segment( sigma1ModifiedRodriguesParametersIndex, 3 ) = conversionSign *
            quaternionElements.segment( epsilon1QuaternionIndex, 3 ) / ( 1 + conversionSign * etaQuaternionParameter );
    convertedModifiedRodriguesParameterElements( shadowFlagModifiedRodriguesParametersIndex ) = shadowFlag ? 1.0 : 0.0;

    // Give output
    return convertedModifiedRodriguesParameterElements;
}

//! Convert modified Rodrigues parameters to quaternions.
Eigen::Vector4d convertModifiedRodriguesParametersToQuaternionElements( const Eigen::Vector4d& modifiedRodriguesParameterElements )
{
    // Declare eventual output vector
    Eigen::Vector4d convertedQuaternionElements;

    // Extract modified Rodrigues parameters
    Eigen::Vector3d modifiedRodriguesParametersVector =
            modifiedRodriguesParameterElements.segment( sigma1ModifiedRodriguesParametersIndex, 3 );

    // Precompute often used variables
    double modifiedRodriguesParametersMagnitudeSquared = modifiedRodriguesParametersVector.squaredNorm( );

    // Convert modified Rodrigues parameters to quaternions
    double conversionSign = ( int( modifiedRodriguesParameterElements( shadowFlagModifiedRodriguesParametersIndex ) ) == 1 ) ?
                - 1.0 : 1.0; // converion is slightly different for SMRP and MRP
    convertedQuaternionElements( etaQuaternionIndex ) = conversionSign *
            ( 1.0 - modifiedRodriguesParametersMagnitudeSquared ) /
            ( 1.0 + modifiedRodriguesParametersMagnitudeSquared );
    convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) = conversionSign *
            2.0 / ( 1.0 + modifiedRodriguesParametersMagnitudeSquared ) *
            modifiedRodriguesParametersVector;

    // Give output
    return convertedQuaternionElements;
}

//! Convert quaternions to exponential map.
Eigen::Vector4d convertQuaternionsToExponentialMapElements( const Eigen::Vector4d& quaternionElements )
{
    using mathematical_constants::PI;

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector4d convertedExponentialMapElements;

    // Convert quaternions to exponential map (or SEM)
    double exponentialMapMagnitude = 2.0 * std::acos( quaternionElements( etaQuaternionIndex ) );
    bool shadowFlag = std::fabs( exponentialMapMagnitude ) > PI;
    Eigen::Vector3d exponentialMapVector = quaternionElements.segment( epsilon1QuaternionIndex, 3 );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        exponentialMapVector *= 48.0 / ( 24.0 - exponentialMapMagnitude * exponentialMapMagnitude );
    }
    else
    {
        exponentialMapVector *= shadowFlag ?
                    - ( 2.0 * PI - exponentialMapMagnitude ) / std::sin( 0.5 * exponentialMapMagnitude ) : // shadow exponential map
                    exponentialMapMagnitude / std::sin( 0.5 * exponentialMapMagnitude ); // exponential map
    }
    convertedExponentialMapElements.segment( e1ExponentialMapIndex, 3 ) = exponentialMapVector;
    convertedExponentialMapElements( shadowFlagExponentialMapIndex ) = shadowFlag ? 1.0 : 0.0; // set flag

    // Give output
    return convertedExponentialMapElements;
}

//! Convert exponential map to quaternions.
Eigen::Vector4d convertExponentialMapToQuaternionElements( const Eigen::Vector4d& exponentialMapElements )
{
    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector4d convertedQuaternionElements;

    // Extract exponential map
    Eigen::Vector3d exponentialMapVector = exponentialMapElements.segment( e1ExponentialMapIndex, 3 );
    double exponentialMapMagnitude = exponentialMapVector.norm( );

    // Convert exponential map to quaternions
    convertedQuaternionElements( etaQuaternionIndex ) = std::cos( 0.5 * exponentialMapMagnitude );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) =
                exponentialMapVector * ( 0.5 - exponentialMapMagnitude * exponentialMapMagnitude / 48.0 );
    }
    else
    {
        convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) =
                exponentialMapVector / exponentialMapMagnitude * std::sin( 0.5 * exponentialMapMagnitude );
    }

    // Change sign based on shadow flag
    double conversionSign = ( int( exponentialMapElements( shadowFlagExponentialMapIndex ) ) == 1 ) ?
                - 1.0 : 1.0; // converion is slightly different for SEM and EM

    // Give output
    return conversionSign * convertedQuaternionElements;
}

} // namespace orbital_element_conversions

} // namespace tudat
