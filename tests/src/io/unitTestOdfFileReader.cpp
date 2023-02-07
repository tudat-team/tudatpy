
#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>

#include <boost/make_shared.hpp>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/math/interpolators/lookupScheme.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/astro/orbit_determination/parseOdfFile.h"

int main( )
{
    using namespace tudat;
    //   boost::shared_ptr< orbit_determination::ProcessedOdfFileContents > odfContents =
    //           orbit_determination::parseOdfFileContents(
    //               input_output::readOdfFile( "/home/dominic/Downloads/mromagr2017_117_0745xmmmv1.odf" ) );

    //   std::map< observation_models::ObservableType,
    //           std::vector< boost::shared_ptr< orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocks =
    //           odfContents->dataBlocks;

    //   for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
    //   {
    //       int counter = 0;
    //       for( unsigned int i = 0; i < it->second.size( ); i++ )
    //       {
    //           boost::shared_ptr< orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
    //                   boost::dynamic_pointer_cast< orbit_determination::ProcessdOdfFileDopplerData >(
    //                       it->second.at( i ) );
    //           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );

    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getCompressionTimes( ),
    //                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getObservationData( ),
    //                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getReferenceFrequencies( ),
    //                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           input_output::writeDataMapToTextFile(
    //                       currentDopplerData->getRampFlags( ),
    //                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
    //           counter++;
    //       }
    //   }

    std::vector< boost::filesystem::path > files = input_output::listAllFilesInDirectory(
                "/home/dominic/Software/MercuryData/odf/", false );

    std::vector< boost::shared_ptr< orbit_determination::ProcessedOdfFileContents > > odfContentsList;

    //for( unsigned int i = 0; i < files.size( ); i++ )
    for( unsigned int i = 0; i < files.size( ); i++ )
    {
        if( i % 100 == 0 )
        {
            //std::cout<<i<<std::endl;
        }
        std::string fileString = files.at( i ).string( );
        int stringSize = fileString.size( );

        if( fileString.substr( stringSize - 3, stringSize -1 ) == "dat" )
        {
            odfContentsList.push_back( orbit_determination::parseOdfFileContents(
                                           input_output::readOdfFile( "/home/dominic/Software/MercuryData/odf/" + fileString ) ) );
        }
    }


    boost::shared_ptr< orbit_determination::ProcessedOdfFileContents > mergedData =
            orbit_determination::mergeOdfFileContents(
                odfContentsList );

    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
            boost::shared_ptr< orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocksMerged =
            mergedData->processedDataBlocks;

    double filterStartTime = 58.0 * 365.25 * 86400.0;
    double filterEndTime = 62.0 * 365.25 * 86400.0;

    double filterFrequency = 200.0E3;

    for( auto it = dataBlocksMerged.begin( ); it != dataBlocksMerged.end( ); it++ )
    {
        int counter = 0;

        std::vector< double > observationTimes, observables, referenceFrequencies;
        std::vector< std::string > receivingStation, transmittingStation;
        std::vector< double > receivingRampFrequency, transmittingRampFrequency;
        std::vector< std::string > originFiles;
        std::vector< bool > rampFlags;


        for( auto linkIt = it->second.begin( ); linkIt != it->second.end( ); linkIt++ )
        {
            boost::shared_ptr< orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
                    boost::dynamic_pointer_cast< orbit_determination::ProcessdOdfFileDopplerData >( linkIt->second );

            std::vector< double > currentObservationTimes = currentDopplerData->observationTimes,
                    currentObservables = currentDopplerData->observableValues,
                    currentReferenceFrequencies = currentDopplerData->referenceFrequency;
            std::vector< std::string > currentOriginFiles = currentDopplerData->originFile;
            std::vector< bool > currentRampFlags = currentDopplerData->rampingFlag;
            std::cout<<it->first<<" "<<linkIt->first.first<<" "<<linkIt->first.second<<" "<<currentObservationTimes.size( )<<std::endl;

            boost::shared_ptr< orbit_determination::RampedReferenceFrequencyInterpolator > transmitterRampInterpolator;
            if( mergedData->rampInterpolators.count(
                        boost::lexical_cast< int >( linkIt->first.first ) ) != 0 )
            {
                transmitterRampInterpolator = mergedData->rampInterpolators.at(
                            boost::lexical_cast< int >( linkIt->first.first ) );
            }
            else
            {
                transmitterRampInterpolator = NULL;
            }

            boost::shared_ptr< orbit_determination::RampedReferenceFrequencyInterpolator > receiverRampInterpolator;
            if( mergedData->rampInterpolators.count(
                        boost::lexical_cast< int >( linkIt->first.second ) ) != 0 )
            {
                receiverRampInterpolator = mergedData->rampInterpolators.at(
                            boost::lexical_cast< int >( linkIt->first.second ) );
            }
            else
            {
                receiverRampInterpolator = NULL;
            }

            for( unsigned int j = 0; j < currentObservationTimes.size( ); j++ )
            {
                if( currentObservationTimes.at( j ) > filterStartTime &&
                        currentObservationTimes.at( j ) < filterEndTime &&
                        std::fabs( currentObservables.at( j ) ) < filterFrequency
                        && std::stoi( linkIt->first.second ) > 0 && std::stoi( linkIt->first.second ) < 100
                         && std::stoi( linkIt->first.first ) > 0 && std::stoi( linkIt->first.first ) < 100 )
                {

                    observationTimes.push_back( currentObservationTimes.at( j ) );
                    observables.push_back( currentObservables.at( j ) );
                    //originFiles.push_back( currentOriginFiles.at( j ) );
                    rampFlags.push_back( currentRampFlags.at( j ) );
                    referenceFrequencies.push_back( currentReferenceFrequencies.at( j ) );

                    receivingStation.push_back( linkIt->first.second );
                    transmittingStation.push_back( linkIt->first.first );

                    bool testBoolean;

                    if( receiverRampInterpolator != NULL )
                    {
                       receivingRampFrequency.push_back( receiverRampInterpolator->getCurrentReferenceFrequency(
                                                             currentObservationTimes.at( j ), testBoolean ) );
                    }
                    else
                    {
                        receivingRampFrequency.push_back( TUDAT_NAN );
                    }

                    if( transmitterRampInterpolator != NULL )
                    {
                       transmittingRampFrequency.push_back( transmitterRampInterpolator->getCurrentReferenceFrequency(
                                                                currentObservationTimes.at( j ), testBoolean ) );
                    }
                    else
                    {
                        transmittingRampFrequency.push_back( TUDAT_NAN );
                    }



                    //           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );

                    //           input_output::writeDataMapToTextFile(
                    //                       currentDopplerData->getCompressionTimes( ),
                    //                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
                    //           input_output::writeDataMapToTextFile(
                    //                       currentDopplerData->getObservationData( ),
                    //                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
                    //           input_output::writeDataMapToTextFile(
                    //                       currentDopplerData->getReferenceFrequencies( ),
                    //                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
                    //           input_output::writeDataMapToTextFile(
                    //                       currentDopplerData->getRampFlags( ),
                    //                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
                    counter++;
                }
            }
        }

        Eigen::MatrixXd dataMatrix = Eigen::MatrixXd( rampFlags.size( ), 8 );
        for( unsigned int j = 0; j < rampFlags.size( ); j++ )
        {
            dataMatrix( j, 0 ) = observationTimes.at( j );
            dataMatrix( j, 1 ) = observables.at( j );
            dataMatrix( j, 2 ) = boost::lexical_cast< double >( receivingStation.at( j ) );
            dataMatrix( j, 3 ) = boost::lexical_cast< double >( transmittingStation.at( j ) );
            dataMatrix( j, 4 ) = referenceFrequencies.at( j );
            dataMatrix( j, 5 ) = static_cast< double >( rampFlags.at( j ) );
            dataMatrix( j, 6 ) = static_cast< double >( receivingRampFrequency.at( j ) );
            dataMatrix( j, 7 ) = static_cast< double >( transmittingRampFrequency.at( j ) );
        }


        input_output::writeMatrixToFile(
                    dataMatrix, "odfFileSummary_" + std::to_string( it->first ), 16, "/home/dominic/Documents/" );
    }


}

