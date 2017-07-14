#ifndef TUDAT_PODPROCESSING_H
#define TUDAT_PODPROCESSING_H

#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a single vector of times from all observation times.
/*!
 *  Function to create a single vector of times from all observation times, created by concatenating all observation times,
 *  in th eorder of observable type and link ends, respectively, as they are stored in the input data type.
 *  \param measurementData Data tructure containing all measurement data, first by observable type, then by link ends.
 *  \return Concatenated vector of times.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< TimeType > getConcatenatedTimeVector(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
{
    // Iterate over all observations and concatenate the time vectors.
    std::vector< TimeType > concatenatedTimes;
    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
    {
        int currentObservableSize = getObservableSize( observablesIterator->first );

        if( currentObservableSize < 0 )
        {
            std::cerr<<"Error when getting osbervable size in pod output; size is dynamics"<<std::endl;
        }

        for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::SingleObservablePodInputType::const_iterator dataIterator =
             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
        {
            if( currentObservableSize == 1 )
            {
                concatenatedTimes.insert( concatenatedTimes.end( ), dataIterator->second.second.first.begin( ),
                                          dataIterator->second.second.first.end( ) );
            }
            else
            {
                std::vector< TimeType > originalTimesList = dataIterator->second.second.first;
                std::vector< TimeType > currentTimesList;
                for( unsigned int i = 0; i < originalTimesList.size( ); i++ )
                {
                    for( unsigned j = 0; j < currentObservableSize; j++ )
                    {
                        currentTimesList.push_back( originalTimesList.at( i ) );
                    }
                }

                concatenatedTimes.insert( concatenatedTimes.end( ), currentTimesList.begin( ),
                                          currentTimesList.end( ) );

            }
        }
    }

    return concatenatedTimes;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::pair< std::vector< int >, std::map< observation_models::LinkEnds, int > > getConcatenatedGroundStationIndex(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
{
    // Iterate over all observations and concatenate the time vectors.
    std::vector< int > concatenatedIds;
    std::map< observation_models::LinkEnds, int > stationIds;

    int maximumStationId = 0;
    int currentStationId;

    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
    {
        for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::SingleObservablePodInputType::const_iterator dataIterator =
             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
        {
            if( stationIds.count( dataIterator->first ) == 0 )
            {
                stationIds[ dataIterator->first  ] = maximumStationId;
                currentStationId = maximumStationId;
                maximumStationId++;
            }
            else
            {
                currentStationId = stationIds[ dataIterator->first ];
            }
            for( unsigned int i = 0; i < dataIterator->second.second.first.size( ); i++ )
            {
                concatenatedIds.push_back( currentStationId );
            }
        }
    }

    return std::make_pair( concatenatedIds, stationIds );
}

template< typename ObservationScalarType = double, typename TimeType = double >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedMeasurementVector(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData )
{
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef OrbitDeterminationManager< ObservationScalarType, TimeType > OdManager;

    int totalNumberOfObservations = OdManager::getNumberOfObservationsPerObservable( measurementData ).second;
    // Iterate over all observations and concatenate the time vectors.
    ObservationVectorType concatenatedObservations = ObservationVectorType::Zero( totalNumberOfObservations, 1 );

    int currentIndex = 0;
    for( typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType::const_iterator observablesIterator =
         measurementData.begin( ); observablesIterator != measurementData.end( ); observablesIterator++ )
    {
        for( typename std::map< observation_models::LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, observation_models::LinkEndType > > >::const_iterator dataIterator =
             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
        {
            concatenatedObservations.segment( currentIndex, dataIterator->second.first.rows( ) ) = dataIterator->second.first;
            currentIndex += dataIterator->second.first.rows( );
        }
    }

    return concatenatedObservations;
}


//! Function to sort the observation matrix by the time of the associated observations.
/*!
 *  Function to sort the observation matrix by the time of the associated observations, in ascending order. That is, the time associated with
 *  the observation wrt which the partial is stored in a given line of the observation matrix is used to determine the new row position of that
 *  row in the ordered matrix.
 *  \param measurementData Set of all measurement data, ordered by observable type and link end set
 *  \param typeAndLinkSortedObservationMatrix Observation matrix, in original ordering (observable type and link ends) as produced by
 *  OrbitDeterminationManager.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::pair< Eigen::MatrixXd, std::vector< TimeType > > getTimeOrderedObservationMatrix(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedObservationMatrix )
{
    //typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    std::vector< TimeType > concatenatedTimes = getConcatenatedTimeVector< ObservationScalarType, TimeType >( measurementData ) ;

    // Sort the concatesnated time vector, and get the order of the sorting.
    std::pair< std::vector< int >, std::vector< TimeType > > sortOutput = utilities::getSortOrderOfVectorAndSortedVector(
                concatenatedTimes );

    // Retrieve sort order and check consistency of its size.
    std::vector< int > timeVectorSortOrder = sortOutput.first;
    if( static_cast< int >( timeVectorSortOrder.size( ) ) != typeAndLinkSortedObservationMatrix.rows( ) )
    {
        std::cerr<<"Error when sorting information matrix by time, sizes incompatible"<<std::endl;
    }

    // Create observation matrix with entries in order of time.
    int numberOfColumns = typeAndLinkSortedObservationMatrix.cols( );
    Eigen::MatrixXd sortedMatrix = Eigen::MatrixXd::Zero( typeAndLinkSortedObservationMatrix.rows( ), typeAndLinkSortedObservationMatrix.cols( ) );
    for( unsigned int i = 0; i < timeVectorSortOrder.size( ); i++ )
    {
        sortedMatrix.block( i, 0, 1, numberOfColumns ) =
                typeAndLinkSortedObservationMatrix.block( timeVectorSortOrder[ i ], 0, 1, numberOfColumns );
    }

    // Return time-sorted observation matrix and time order.
    return std::make_pair( sortedMatrix, sortOutput.second );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::MatrixXd > calculateCovarianceMatrixAsFunctionOfTime(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedNormalizedObservationMatrix,
        const Eigen::VectorXd& normalizationFactors,
        const double outputTimeStep,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& normalizedInverseAPrioriCovariance )
{
    // Check consistency of input data
    if( normalizedInverseAPrioriCovariance.cols( ) != normalizedInverseAPrioriCovariance.rows( ) )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 1"<<std::endl;
    }

    int totalNumberOfParameters = normalizedInverseAPrioriCovariance.cols( );

    if( typeAndLinkSortedNormalizedObservationMatrix.cols( ) != totalNumberOfParameters )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 3"<<std::endl;
    }

    if( normalizationFactors.rows( ) != totalNumberOfParameters )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 4"<<std::endl;
    }
    if( typeAndLinkSortedNormalizedObservationMatrix.rows( ) != diagonalOfWeightMatrix.rows( ) )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 5"<<std::endl;
    }

    // Order observation matrix by time of observations
    std::pair< Eigen::MatrixXd, std::vector< TimeType > > timeOrderedMatrixOutput =
            getTimeOrderedObservationMatrix< ObservationScalarType, TimeType >(
                measurementData, typeAndLinkSortedNormalizedObservationMatrix );
    std::vector< TimeType > timeOrder = timeOrderedMatrixOutput.second;

    // Create lookupn scheme for time value
    interpolators::BinarySearchLookupScheme< TimeType > timeLookup =
            interpolators::BinarySearchLookupScheme< TimeType >( timeOrder );

    // Create matrix to unnormalize parameters.
    Eigen::MatrixXd unnormalizationMatrix = normalizationFactors.asDiagonal( );

    // Declare return map.
    std::map< TimeType, Eigen::MatrixXd > covarianceMatrixHistory;

    // Initialize loop variables.
    TimeType currentTime = timeOrder[ 0 ];
    int currentIndex;
    Eigen::MatrixXd currentInverseCovarianceMatrix;

    // Loop over matrix at given time interval.
    while( currentTime < timeOrder[ timeOrder.size( ) - 1 ] )
    {
        // Increment time (done at beginning of loop so that one entry is made for entire matrix, i.e. when current time > largest time)
        currentTime += outputTimeStep;

        // Find index in list of times.
        currentIndex = timeLookup.findNearestLowerNeighbour( currentTime );

        // Create information matrix up to current time.
        Eigen::MatrixXd currentInformationMatrix = timeOrderedMatrixOutput.first.block(
                    0, 0, currentIndex, timeOrderedMatrixOutput.first.cols( ) );

        // Create inverse of covariance matrix
        currentInverseCovarianceMatrix = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
                    currentInformationMatrix, diagonalOfWeightMatrix.segment( 0, currentIndex ),
                    unnormalizationMatrix.inverse( ) * normalizedInverseAPrioriCovariance * unnormalizationMatrix.inverse( ) );
        covarianceMatrixHistory[ currentTime ] = unnormalizationMatrix.inverse( ) *
                currentInverseCovarianceMatrix.inverse( ) * unnormalizationMatrix.inverse( );
    }


    return covarianceMatrixHistory;
}

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          typename ParameterScalarType = double >
std::map< TimeType, Eigen::MatrixXd >  calculateCovarianceMatrixAsFunctionOfTime(
        const boost::shared_ptr< PodInput< ObservationScalarType, TimeType > >& podInputData,
        const boost::shared_ptr< PodOutput< ParameterScalarType > >& podOutputData,
        const double outputTimeStep )
{
    return calculateCovarianceMatrixAsFunctionOfTime< ObservationScalarType, TimeType >(
                podInputData->getObservationsAndTimes( ), podOutputData->normalizedInformationMatrix_,
                podOutputData->informationMatrixTransformationDiagonal_, outputTimeStep,
                podOutputData->weightsMatrixDiagonal_, podInputData->getInverseOfAprioriCovariance( ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void printPodData(
        const typename OrbitDeterminationManager< ObservationScalarType, TimeType >::PodInputType& measurementData,
        const Eigen::MatrixXd& typeAndLinkSortedNormalizedObservationMatrix,
        const Eigen::VectorXd& normalizationFactors,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& normalizedInverseAPrioriCovariance,
        const std::string folderName,
        const std::string fileNameSuffix,
        const bool sortInformationMatrixByObservableFirst = 0 )
{
    // Check consistency of input data
    if( normalizedInverseAPrioriCovariance.cols( ) != normalizedInverseAPrioriCovariance.rows( ) )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 1"<<std::endl;
    }

    int totalNumberOfParameters = normalizedInverseAPrioriCovariance.cols( );

    if( typeAndLinkSortedNormalizedObservationMatrix.cols( ) != totalNumberOfParameters )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 3"<<std::endl;
    }

    if( normalizationFactors.rows( ) != totalNumberOfParameters )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 4"<<std::endl;
    }
    if( typeAndLinkSortedNormalizedObservationMatrix.rows( ) != diagonalOfWeightMatrix.rows( ) )
    {
        std::cerr<<"Error when calculating covariance as function of time, error 5"<<std::endl;
    }

    // Order observation matrix by time of observations
    std::pair< Eigen::MatrixXd, std::vector< TimeType > > informationMatrix;
    if( !sortInformationMatrixByObservableFirst )
    {
        informationMatrix =
                getTimeOrderedObservationMatrix< ObservationScalarType, TimeType >( measurementData, typeAndLinkSortedNormalizedObservationMatrix );
        input_output::writeMatrixToFile(
                    informationMatrix.first,
                    folderName + "timeOrderedInformationMatrix"+fileNameSuffix+".dat" );
    }
    else
    {
        input_output::writeMatrixToFile(
                    typeAndLinkSortedNormalizedObservationMatrix,
                    folderName + "informationMatrix"+fileNameSuffix+".dat" );
    }


    input_output::writeMatrixToFile(
                normalizationFactors,
                folderName + "partialNormalizationFactors"+fileNameSuffix+".dat" );

}


template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType,
          typename ParameterScalarType = double >
void printPodData(
        const boost::shared_ptr< PodInput< ObservationScalarType, TimeType > >& podInputData,
        const boost::shared_ptr< PodOutput< ParameterScalarType > >& podOutputData,
        const std::string folderName,
        const std::string fileNameSuffix,
        const bool sortInformationMatrixByObservableFirst = 0 )
{
    return printPodData< ObservationScalarType, TimeType >(
                podInputData->observationsAndTimes_, podOutputData->normalizedInformationMatrix_,
                podOutputData->informationMatrixTransformationDiagonal_,
                podOutputData->weightsMatrixDiagonal_, podInputData->inverseOfAprioriCovariance_,
                folderName, fileNameSuffix, sortInformationMatrixByObservableFirst );
}

}

}
#endif // TUDAT_PODPROCESSING_H
