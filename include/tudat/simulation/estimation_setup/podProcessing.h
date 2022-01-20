/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PODPROCESSING_H
#define TUDAT_PODPROCESSING_H

#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

////! Function to create a single vector of times from all observation times.
///*!
// *  Function to create a single vector of times from all observation times, created by concatenating all observation times,
// *  in the order of first observable type and the link ends, as they are stored in the input data type (PodInputType).
// *  \param measurementData Data structure containing all measurement data, first by observable type, then by link ends.
// *  \return Concatenated vector of times.
// */
//template< typename ObservationScalarType = double, typename TimeType = double >
//std::vector< TimeType > getConcatenatedTimeVector(
//        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
//{
//    // Iterate over all observations and concatenate the time vectors.
//    std::vector< TimeType > concatenatedTimes;
//    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator
//         observablesIterator = measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
//    {
//        // Get data type size.
//        int currentObservableSize = getObservableSize( observablesIterator->first );

//        // Iterate over all link ends
//        for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::SingleObservablePodInputType::const_iterator
//             dataIterator = observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( );
//             dataIterator++  )
//        {
//            // Copy data in the case that observable size is 1
//            if( currentObservableSize == 1 )
//            {
//                concatenatedTimes.insert( concatenatedTimes.end( ), dataIterator->second.second.first.begin( ),
//                                          dataIterator->second.second.first.end( ) );
//            }
//            // For observable size N>1, coopy time tag N times for each observation
//            else
//            {
//                std::vector< TimeType > originalTimesList = dataIterator->second.second.first;
//                std::vector< TimeType > currentTimesList;
//                for( unsigned int i = 0; i < originalTimesList.size( ); i++ )
//                {
//                    for( int j = 0; j < currentObservableSize; j++ )
//                    {
//                        currentTimesList.push_back( originalTimesList.at( i ) );
//                    }
//                }

//                concatenatedTimes.insert( concatenatedTimes.end( ), currentTimesList.begin( ),
//                                          currentTimesList.end( ) );

//            }
//        }
//    }

//    return concatenatedTimes;
//}

////! Function to create a single vector of link end indices for all observations
///*!
// *  Function to create a single vector of link end indices for all observations, created by assigning an integer to each set of
// *  link ends, and these indices for all observations, in the order of first observable type and the link ends,
// *  as they are stored in the input data type (PodInputType). Output vector has same size as vector of observables, so
// *  for an observation of size N, N values are added to the output vector
// *  \param measurementData Data structure containing all measurement data, first by observable type, then by link ends.
// *  \return Concatenated vector of link end indices (first), mapping of LinkEnds to index (second)
// */
//template< typename ObservationScalarType = double, typename TimeType = double >
//std::pair< std::vector< int >, std::map< observation_models::LinkEnds, int > > getConcatenatedGroundStationIndex(
//        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
//{
//    std::vector< int > concatenatedIds;
//    std::map< observation_models::LinkEnds, int > stationIds;

//    int maximumStationId = 0;
//    int currentStationId;

//    // Iterate over all observations and concatenate the link end indices.
//    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
//         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
//    {
//        // Get current observable size
//        int currentObservableSize = getObservableSize( observablesIterator->first );

//        for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::SingleObservablePodInputType::const_iterator
//             dataIterator =
//             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
//        {
//            // Check if current LinkEnds have an index assign to them, and create it if not
//            if( stationIds.count( dataIterator->first ) == 0 )
//            {
//                stationIds[ dataIterator->first  ] = maximumStationId;
//                currentStationId = maximumStationId;
//                maximumStationId++;
//            }
//            else
//            {
//                currentStationId = stationIds[ dataIterator->first ];
//            }

//            // Add indices for current observables
//            for( unsigned int i = 0; i < dataIterator->second.second.first.size( ); i++ )
//            {
//                for( int j = 0; j < currentObservableSize; j++ )
//                {
//                    concatenatedIds.push_back( currentStationId );
//                }
//            }
//        }
//    }

//    return std::make_pair( concatenatedIds, stationIds );
//}

////! Function to create a single vector of observable types for all observations
///*!
// *  Function to create a single vector of observable types for all observations, in the order of first observable type and the
// *  link ends, as they are stored in the input data type (PodInputType). Output vector has same size as vector of observables, so
// *  for an observation of size N, N values are added to the output vector
// *  \param measurementData Data structure containing all measurement data, first by observable type, then by link ends.
// *  \return Concatenated vector of observable types.
// */
//template< typename ObservationScalarType = double, typename TimeType = double >
//std::vector< int > getConcatenatedObservableTypes(
//        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
//{
//    std::vector< int > concatenatedIds;

//    // Iterate over all observations and concatenate the time vectors.
//    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
//         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
//    {
//        int currentObservable = static_cast< int >( observablesIterator->first );
//        int currentObservableSize = getObservableSize( observablesIterator->first );

//        for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::SingleObservablePodInputType::const_iterator
//             dataIterator =
//             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
//        {
//            // Add indices for current observables
//            for( unsigned int i = 0; i < dataIterator->second.second.first.size( ); i++ )
//            {
//                for( int j = 0; j < currentObservableSize; j++ )
//                {
//                    concatenatedIds.push_back( currentObservable );
//                }
//            }
//        }
//    }

//    return concatenatedIds;
//}

////! Function to create a single vector of observations from full observation input data.
///*!
// *  Function to create a single vector of observations from full observation input data, created by concatenating all observations
// *  in the order of first observable type and the link ends, as they are stored in the input data type (PodInputType).
// *  \param measurementData Data structure containing all measurement data, first by observable type, then by link ends.
// *  \return Concatenated vector of observations.
// */
//template< typename ObservationScalarType = double, typename TimeType = double >
//Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedMeasurementVector(
//        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
//{
//    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
//    typedef OrbitDeterminationManager< ObservationScalarType, TimeType > OdManager;

//    // Get total number of observations
//    int totalNumberOfObservations = OdManager::getNumberOfObservationsPerObservable( measurementData ).second;

//    // Iterate over all observations and concatenate the time vectors.
//    ObservationVectorType concatenatedObservations = ObservationVectorType::Zero( totalNumberOfObservations, 1 );

//    // Iterate over all link ends
//    int currentIndex = 0;
//    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
//         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
//    {
//        for( typename std::map< observation_models::LinkEnds, std::pair< ObservationVectorType,
//             std::pair< std::vector< TimeType >, observation_models::LinkEndType > > >::const_iterator dataIterator =
//             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
//        {
//            concatenatedObservations.segment( currentIndex, dataIterator->second.first.rows( ) ) = dataIterator->second.first;
//            currentIndex += dataIterator->second.first.rows( );
//        }
//    }

//    return concatenatedObservations;
//}


//! Function to sort the information matrix by the time of the associated observations.
/*!
 *  Function to sort the information matrix (e.g. matrix with partial derivatives of observables w.r.t. estimated parameters) by
 *  the time of the associated observations, in ascending order. That is, the time associated withthe observation wrt which the
 *  partial is stored in a given line of the information matrix is used to determine the new row position of that
 *  row in the ordered matrix. The lowest time is the first row; the highest time is the last row.
 *  \param measurementData Set of all measurement data, ordered by observable type and link end set
 *  \param typeAndLinkSortedInformationMatrix information matrix, in original ordering (observable type and link ends)
 *  as produced by OrbitDeterminationManager.
 *  \param timeOrder Order in which the observation time vector was sorted to obtain the time-sorted information matrix (returned
 *  by reference_
 *  \return Information matrix, ordered by associated observation time
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::MatrixXd, std::vector< TimeType > > getTimeOrderedInformationMatrix(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedInformationMatrix,
        std::vector< int >& timeOrder )
{
    // Retrieve unordered vector of times
    std::vector< TimeType > concatenatedTimes = measurementData->getConcatenatedTimeVector( );

    // Sort the concatesnated time vector, and get the order of the sorting.
    std::pair< std::vector< int >, std::vector< TimeType > > sortOutput = utilities::getSortOrderOfVectorAndSortedVector(
                concatenatedTimes );
    timeOrder = sortOutput.first;

    // Retrieve sort order and check consistency of its size.
    std::vector< int > timeVectorSortOrder = sortOutput.first;
    if( static_cast< int >( timeVectorSortOrder.size( ) ) != typeAndLinkSortedInformationMatrix.rows( ) )
    {
        std::cout<<timeVectorSortOrder.size( )<<std::endl<<typeAndLinkSortedInformationMatrix.rows( )<<std::endl;
        throw std::runtime_error( "Error when sorting information matrix by time, sizes incompatible" );
    }

    // Create information matrix with entries in order of time.
    int numberOfColumns = typeAndLinkSortedInformationMatrix.cols( );
    Eigen::MatrixXd sortedMatrix = Eigen::MatrixXd::Zero( typeAndLinkSortedInformationMatrix.rows( ),
                                                          typeAndLinkSortedInformationMatrix.cols( ) );
    for( unsigned int i = 0; i < timeVectorSortOrder.size( ); i++ )
    {
        sortedMatrix.block( i, 0, 1, numberOfColumns ) =
                typeAndLinkSortedInformationMatrix.block( timeVectorSortOrder[ i ], 0, 1, numberOfColumns );
    }

    // Return time-sorted information matrix and time order.
    return std::make_pair( sortedMatrix, sortOutput.second );
}

//! Function to create a map of the estimation covariance as a function of time
/*!
 *  Function to create a map of the estimation covariance as a function of time
 *  \param measurementData Data structure containing all observable values, as well as associated times and reference link ends
 *  \param typeAndLinkSortedNormalizedInformationMatrix Information matrix, normalized by the normalizationFactors, and
 *  sorted as in the normalizationFactors: first by observable type, then by link ends
 *  \param normalizationFactors Values by which the parameters (and partials) have been normalized, in order to stabilize
 *  the solution of the normal equations
 *  \param outputTimes Times at which the covariance is to be computed for the output map
 *  \param diagonalOfWeightMatrix Vector containing the diagonal of the weights matrix used in the estimation
 *  \param unnormalizedInverseAPrioriCovariance Inverse a priori covariance matrix, with parameters not normalized by
 *  normalizationFactors
 *  \return Covariance (map values) as a function of time (map keys) for the given estimation input settings and output times.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::MatrixXd > calculateCovarianceUsingDataUpToEpoch(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedNormalizedInformationMatrix,
        const Eigen::VectorXd& normalizationFactors,
        const std::vector< double >& outputTimes,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& unnormalizedInverseAPrioriCovariance )
{
    int totalNumberOfParameters = unnormalizedInverseAPrioriCovariance.cols( );

    Eigen::MatrixXd normalizedInverseAPrioriCovariance = Eigen::MatrixXd(  totalNumberOfParameters, totalNumberOfParameters );
    for( int j = 0; j < totalNumberOfParameters; j++ )
    {
        for( int k = 0; k < totalNumberOfParameters; k++ )
        {
            normalizedInverseAPrioriCovariance( j, k ) = unnormalizedInverseAPrioriCovariance( j, k ) /
                    ( normalizationFactors( j ) * normalizationFactors( k ) );
        }
    }

    // Check consistency of input data
    if( normalizedInverseAPrioriCovariance.cols( ) != normalizedInverseAPrioriCovariance.rows( ) )
    {
        throw std::runtime_error( "Error when calculating covariance as function of time, a priori covariance is not square" );
    }

    if( typeAndLinkSortedNormalizedInformationMatrix.cols( ) != totalNumberOfParameters )
    {
        throw std::runtime_error(
                    "Error when calculating covariance as function of time, number of parameters is inconsistent with information matrix" );
    }

    if( normalizationFactors.rows( ) != totalNumberOfParameters )
    {
        throw std::runtime_error(
                    "Error when calculating covariance as function of time, number of parameters is inconsistent with normalization factors" );
    }
    if( typeAndLinkSortedNormalizedInformationMatrix.rows( ) != diagonalOfWeightMatrix.rows( ) )
    {
        throw std::runtime_error(
                    "Error when calculating covariance as function of time, weights are inconsistent with partials" );
    }

    // Order information matrix by time of observations
    std::vector< int > timeOrder;
    std::pair< Eigen::MatrixXd, std::vector< TimeType > > timeOrderedMatrixOutput =
            getTimeOrderedInformationMatrix< ObservationScalarType, TimeType >(
                measurementData, typeAndLinkSortedNormalizedInformationMatrix, timeOrder );
    std::vector< TimeType > orderedTimeVector = timeOrderedMatrixOutput.second;

    Eigen::VectorXd timeOrderedDiagonalOfWeightMatrix = Eigen::VectorXd::Zero( diagonalOfWeightMatrix.rows( ) );
    for( unsigned int i = 0; i < timeOrder.size( ); i++ )
    {
        timeOrderedDiagonalOfWeightMatrix( i ) = diagonalOfWeightMatrix( timeOrder.at( i ) );
    }

    // Create lookupn scheme for time value
    interpolators::BinarySearchLookupScheme< TimeType > timeLookup =
            interpolators::BinarySearchLookupScheme< TimeType >( orderedTimeVector );

    // Declare return map.
    std::map< TimeType, Eigen::MatrixXd > covarianceMatrixHistory;

    // Initialize loop variables.
    TimeType currentTime = orderedTimeVector[ 0 ];
    unsigned int currentIndex;
    Eigen::MatrixXd currentInverseNormalizedCovarianceMatrix;

    // Loop over matrix at given time interval.
    for( unsigned int i = 0; i < outputTimes.size( ); i++ )
    {
        // Increment time (no covariance computed for t=t0)
        currentTime = outputTimes.at( i );

        // Find index in list of times.
        currentIndex = timeLookup.findNearestLowerNeighbour( currentTime );

        if( currentIndex != orderedTimeVector.size( ) - 1 )
        {
            while( orderedTimeVector.at( currentIndex) == orderedTimeVector.at( currentIndex + 1 ) )
            {
                currentIndex++;
                if( currentIndex == orderedTimeVector.size( ) - 1 )
                {
                    break;
                }
            }
        }

        if( orderedTimeVector.size( ) <= currentIndex )
        {
            throw std::runtime_error( "Error when getting covariance as a function of time, output time not found" );
        }

        // Create information matrix up to current time.
        Eigen::MatrixXd currentInformationMatrix = timeOrderedMatrixOutput.first.block(
                    0, 0, currentIndex + 1, timeOrderedMatrixOutput.first.cols( ) );

        // Create inverse of covariance matrix
        currentInverseNormalizedCovarianceMatrix = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
                    currentInformationMatrix, timeOrderedDiagonalOfWeightMatrix.segment( 0, currentIndex + 1 ), normalizedInverseAPrioriCovariance );
        Eigen::MatrixXd covarianceMatrix = currentInverseNormalizedCovarianceMatrix.inverse( );

        for( int i = 0; i < covarianceMatrix.rows( ); i++ )
        {
            for( int j = 0; j < covarianceMatrix.rows( ); j++ )
            {
                covarianceMatrix( i, j ) /= normalizationFactors( i ) * normalizationFactors( j );
            }
        }

        covarianceMatrixHistory[ currentTime ] = covarianceMatrix;
    }


    return covarianceMatrixHistory;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::MatrixXd > calculateCovarianceUsingDataUpToEpoch(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedNormalizedInformationMatrix,
        const Eigen::VectorXd& normalizationFactors,
        const double outputTimeStep,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& unnormalizedInverseAPrioriCovariance )
{
    Eigen::VectorXd timeVector =
            utilities::convertStlVectorToEigenVector( measurementData->getConcatenatedTimeVector( ) );
    double minimumTime = timeVector.minCoeff( );
    double maximumTime = timeVector.maxCoeff( );
    double currentTime = minimumTime;
    std::vector< double > outputTimes;

    // Loop over matrix at given time interval.
    while( currentTime < maximumTime )
    {
        // Increment time (no covariance computed for t=t0)
        currentTime += outputTimeStep;
        outputTimes.push_back( currentTime );
    }

    return calculateCovarianceUsingDataUpToEpoch(
                measurementData, typeAndLinkSortedNormalizedInformationMatrix, normalizationFactors,
                outputTimes, diagonalOfWeightMatrix, unnormalizedInverseAPrioriCovariance );
}

//! Function to create a map of the estimation covariance as a function of time
/*!
 *  Function to create a map of the estimation covariance as a function of time
 *  \param podInputData Data structure containing all input to the orbit determination process.
 *  \param podOutputData Data structure containing all output of the orbit determination process.
 *  \param outputTimes Times at which the covariance is to be computed for the output map.
 *  \return Covariance (map values) as a function of time (map keys) for the given estimation input settings and output times.
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          typename ParameterScalarType = double >
std::map< TimeType, Eigen::MatrixXd >  calculateCovarianceUsingDataUpToEpoch(
        const std::shared_ptr< PodInput< ObservationScalarType, TimeType > >& podInputData,
        const std::shared_ptr< PodOutput< ParameterScalarType > >& podOutputData,
        const std::vector< double >& outputTimes )
{
    return calculateCovarianceUsingDataUpToEpoch< ObservationScalarType, TimeType >(
                podInputData->getObservationsAndTimes( ), podOutputData->normalizedInformationMatrix_,
                podOutputData->informationMatrixTransformationDiagonal_, outputTimes,
                podOutputData->weightsMatrixDiagonal_, podInputData->getInverseOfAprioriCovariance( ) );
}

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          typename ParameterScalarType = double >
std::map< TimeType, Eigen::MatrixXd >  calculateCovarianceUsingDataUpToEpoch(
        const std::shared_ptr< PodInput< ObservationScalarType, TimeType > >& podInputData,
        const std::shared_ptr< PodOutput< ParameterScalarType > >& podOutputData,
        const double outputTimeStep )
{
    return calculateCovarianceUsingDataUpToEpoch< ObservationScalarType, TimeType >(
                podInputData->getObservationCollection( ), podOutputData->normalizedInformationMatrix_,
                podOutputData->informationMatrixTransformationDiagonal_, outputTimeStep,
                podOutputData->weightsMatrixDiagonal_, podInputData->getInverseOfAprioriCovariance( ) );
}


}

}
#endif // TUDAT_PODPROCESSING_H
