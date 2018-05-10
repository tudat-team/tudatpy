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
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      Wertz, J. R. Mission geometry; orbit and constellation design and management.
 *      Mengali, G., Quarta, A.A. Fondamenti di meccanica del volo spaziale.
 *      Wertz, J.R. Mission Geometry; Orbit and Constellation Design and Management, Spacecraft
 *          Orbit and Attitude Systems, Microcosm Press, Kluwer Academic Publishers, 2001.
 *      Advanced Concepts Team, ESA. Keplerian Toolbox, http://sourceforge.net/projects/keptoolbox,
 *          last accessed: 21st April, 2012.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
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
    convertedQuaternionElements( etaQuaternionIndex ) = conversionSign * ( 1.0 - modifiedRodriguesParametersMagnitudeSquared ) /
            ( 1.0 + modifiedRodriguesParametersMagnitudeSquared );
    convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) = conversionSign *
            2.0 / ( 1.0 + modifiedRodriguesParametersMagnitudeSquared ) *
            modifiedRodriguesParametersVector;

    // Give output
    return convertedQuaternionElements;
}

//! Convert quaternions to exponential map.
Eigen::Vector3d convertQuaternionsToExponentialMapElements( const Eigen::Vector4d& quaternionElements )
{
    using mathematical_constants::PI;

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector3d convertedExponentialMapElements = quaternionElements.segment( epsilon1QuaternionIndex, 3 );

    // Convert quaternions to exponential map (or SEM)
    double exponentialMapMagnitude = 2.0 * std::acos( quaternionElements( etaQuaternionIndex ) );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        convertedExponentialMapElements *= 48.0 / ( 24.0 - exponentialMapMagnitude * exponentialMapMagnitude );
    }
    else
    {
        convertedExponentialMapElements *= ( std::fabs( exponentialMapMagnitude ) > PI ) ?
                    - ( 2.0 * PI - exponentialMapMagnitude ) / std::sin( 0.5 * exponentialMapMagnitude ) : // shadow exponential map
                    exponentialMapMagnitude / std::sin( 0.5 * exponentialMapMagnitude ); // exponential map
    }

    // Give output
    return convertedExponentialMapElements;
}

//! Convert exponential map to quaternions.
Eigen::Vector4d convertExponentialMapToQuaternionElements( const Eigen::Vector3d& exponentialMapElements )
{
    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Declare eventual output vector
    Eigen::Vector4d convertedQuaternionElements;

    // Extract exponential map
    double exponentialMapMagnitude = exponentialMapElements.norm( );

    // Convert exponential map to quaternions
    convertedQuaternionElements( etaQuaternionIndex ) = std::cos( 0.5 * exponentialMapMagnitude );
    if ( std::fabs( exponentialMapMagnitude ) < singularityTolerance )
    {
        convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) =
                exponentialMapElements * ( 0.5 - exponentialMapMagnitude * exponentialMapMagnitude / 48.0 );
    }
    else
    {
        convertedQuaternionElements.segment( epsilon1QuaternionIndex, 3 ) =
                exponentialMapElements / exponentialMapMagnitude * std::sin( 0.5 * exponentialMapMagnitude );
    }

    // Give output
    return convertedQuaternionElements;
}

} // namespace orbital_element_conversions

} // namespace tudat
