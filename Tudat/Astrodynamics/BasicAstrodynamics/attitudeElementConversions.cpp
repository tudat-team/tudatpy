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
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

#include "Tudat/Basics/utilities.h"

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

//! Transform quaternion to opposite rotation.
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::Vector4d >& quaternionHistoryMap )
{
    // Get total number of rows
    unsigned int numberOfRows = quaternionHistoryMap.size( );

    // Create concatenation of vectors
    Eigen::VectorXd timeHistoryVector;
    Eigen::MatrixXd quaternionHistoryMatrix;
    timeHistoryVector.resize( numberOfRows, 1 );
    quaternionHistoryMatrix.resize( 4, numberOfRows );
    int i = 0;
    for ( std::map< double, Eigen::Vector4d >::const_iterator mapIterator = quaternionHistoryMap.begin( );
          mapIterator != quaternionHistoryMap.end( ); mapIterator++ )
    {
        timeHistoryVector[ i ] = mapIterator->first;
        quaternionHistoryMatrix.col( i ) = mapIterator->second;
        i++;
    }

    // Compute numerical derivative of first quaternion elements
    std::cout << "Rows: " << quaternionHistoryMatrix.rows( ) << ", cols: " << quaternionHistoryMatrix.cols( ) << std::endl;
    Eigen::VectorXd vectorOfEtas = quaternionHistoryMatrix.row( etaQuaternionIndex ).transpose( );
    Eigen::VectorXd timeDerivativeOfEtas;
    timeDerivativeOfEtas.resize( numberOfRows - 1, 1 );
    for ( unsigned int i = 0; i < numberOfRows - 1; i++ )
    {
        timeDerivativeOfEtas[ i ] =
                ( vectorOfEtas[ i + 1 ] - vectorOfEtas[ i ] ) /
                ( timeHistoryVector[ i + 1 ] - timeHistoryVector[ i ] );
    }

    // Get standard deviation of first elements
    double conversionThreshold = statistics::computeStandardDeviationOfVectorComponents( timeDerivativeOfEtas );
}

//! Transform quaternion in translational or rotational state to opposite rotation.
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::VectorXd >& stateHistoryMap,
                                           const propagators::IntegratedStateType stateType )
{
    // Select index based on input state type
    int quaternionStartIndex;
    switch ( stateType )
    {
    case propagators::translational_state:
        quaternionStartIndex = static_cast< int >( etaUSM7Index );
        break;
    case propagators::rotational_state:
        quaternionStartIndex = static_cast< int >( etaQuaternionIndex );
        break;
    default:
        throw std::runtime_error( "Error in conversion of quaternion history."
                                  "Only translational and rotational propagators are supported." );
    }

    // Extract exponential map from state and convert history
    std::map< double, Eigen::Vector4d > quaternionHistoryMap;
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = stateHistoryMap.begin( );
          mapIterator != stateHistoryMap.end( ); mapIterator++ )
    {
        quaternionHistoryMap[ mapIterator->first ] = mapIterator->second.segment( quaternionStartIndex, 4 );
    }
    convertQuaternionHistoryToMatchSigns( quaternionHistoryMap );

    // Replace state elements with new quaterions
    for ( std::map< double, Eigen::VectorXd >::iterator mapIterator = stateHistoryMap.begin( );
          mapIterator != stateHistoryMap.end( ); mapIterator++ )
    {
        mapIterator->second.segment( quaternionStartIndex, 4 ) = quaternionHistoryMap[ mapIterator->first ];
    }
}

} // namespace orbital_element_conversions

} // namespace tudat
