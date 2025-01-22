/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E, 2008, JPL/DSN
 */

#ifndef TUDAT_PROCESSTRACKINGTXTFILE_H
#define TUDAT_PROCESSTRACKINGTXTFILE_H

#include "tudat/basics/utilities.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{
namespace observation_models
{

/*!
 * Map containing all the tracking data types that correspond to an observable. If a specific set of TrackingDataTypes (keys)
 * Is present in the file, the corresponding ObservableType (value) can later be added to the  ObservationCollection
 */
static const std::map<ObservableType, std::vector<input_output::TrackingDataType>> observableRequiredDataTypesMap = {
    {n_way_range, {input_output::TrackingDataType::n_way_light_time}},
    {doppler_measured_frequency, {input_output::TrackingDataType::doppler_measured_frequency}},
    {dsn_n_way_averaged_doppler, {input_output::TrackingDataType::doppler_averaged_frequency}},
};


/*!
 * Function to extract a vector of all the observableTypes that can be deduced from a list of DataTypes
 * @param availableDataTypes Vector of available TrackingDataTypes
 * @return vector of observable types that can be derived from the tracking data types
 */
std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes);

//! Class containing the processed file contents for a Tracking txt data file
template< typename ObservationScalarType = double, typename TimeType = double >
class ProcessedTrackingTxtFileContents
{
public:
  /*!
   * Constructor
   * @param rawTrackingTxtContents Raw contents object
   * @param spacecraftName Name used to refer to the spacecraft of interest
   * @param earthFixedGroundStationPositions map of earth fixed positions of the ground stations.
   */
  ProcessedTrackingTxtFileContents(const std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents,
                                   const std::string spacecraftName,
                                   const std::map<std::string, Eigen::Vector3d>& earthFixedGroundStationPositions)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName),
        earthFixedGroundStationPositions_(earthFixedGroundStationPositions), initialised_(false)
  {
        initialise();
  }

  //! Main initialisation sequence to create the processed file contents.
  void initialise()
  {
    updateLinkEnds();
    updateObservationTimes(); // updateLinkEnds must come first
    updateAvailableObservableTypes();
    updateObservations();
    initialised_ = true;
  }

// Update  information from the raw file data
public:

  //! Method to extract the first and last time of the observations
  std::pair< TimeType, TimeType> getStartAndEndTime()
  {
    const std::vector< TimeType >& times = getObservationTimes();
    std::pair<TimeType, TimeType> startEndTime({ times.front(), times.back() });
    return startEndTime;
  }

  //! Method to create a vector of the linkEnds for each of the lines in the file (different ground stations might be used)
  void updateLinkEnds()
  {
      // Clear any previous values
      linkEndsVector_.clear();

      // Get information from raw data file
      const auto& metaDataStrMap = rawTrackingTxtFileContents_->getMetaDataStrMap();
      const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();

      // Deduce linkends representation
      LinkEndsRepresentation linkEndsRepresentation = getLinkEndsRepresentation();

      // Create a vector of LinkEnds based on how they are represented
      // This currently only implements the DSN transmitter and receiver
      switch (linkEndsRepresentation)
      {
      case dsn_transmitting_receiving_station_nr:
      {
          const auto& dsnTransmitterIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_transmitting_station_nr);
          const auto& dsnReceiverIds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::dsn_receiving_station_nr);

          for (size_t i = 0; i < numDataRows; ++i)
          {
              std::string transmitterName = getStationNameFromStationId(dsnTransmitterIds.at(i));
              std::string receiverName = getStationNameFromStationId(dsnReceiverIds.at(i));
              LinkEnds currentLinkEnds{
                  {transmitter, LinkEndId("Earth", transmitterName)},
                  {reflector, LinkEndId(spacecraftName_, "")},
                  {receiver, LinkEndId("Earth", receiverName)} };
              linkEndsVector_.push_back(currentLinkEnds);
          }
          break;
      }

      case transmitting_receiving_station_name:
      {
          std::string receivingStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::receiving_station_name);
          std::string transmittingStationName = rawTrackingTxtFileContents_->getMetaDataStrMap().at(input_output::TrackingDataType::transmitting_station_name);

          for (size_t i = 0; i < numDataRows; ++i)
          {
              LinkEnds currentLinkEnds{
                  {transmitter, LinkEndId("Earth", transmittingStationName)},
                  {reflector, LinkEndId(spacecraftName_, "")},
                  {receiver, LinkEndId("Earth", receivingStationName)},
              };
              linkEndsVector_.push_back(currentLinkEnds);
          }
          break;
      }

          // Throw error if representation not implemented
      default: {
          throw std::runtime_error("Error while processing tracking txt file: LinkEnds representation not recognised or implemented.");
      }
      }

      // Creating a set with all the distinct LinkEnds
      linkEndsSet_ = utilities::vectorToSet(linkEndsVector_);
  }

  //! Method to (re)calculate the times from the Tracking file
  void updateObservationTimes()
  {
      // Clear any previous values
      observationTimes_.clear();

      // Get data map and time representation
      const auto& numDataRows = rawTrackingTxtFileContents_->getNumRows();
      TimeRepresentation timeRepresentation = getTimeRepresentation();

      // Depending on the time representation, convert further to tdb seconds since j2000
      switch (timeRepresentation)
      {
      case tdb_seconds_j2000:
      {
          observationTimes_ =
              utilities::staticCastVector< TimeType, double >(
                  rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::tdb_reception_time_j2000) );
          break;
      }
      case utc_seconds_j2000:
      {
          observationTimesUtc_ =
              utilities::staticCastVector< TimeType, double >(
                  rawTrackingTxtFileContents_->getDoubleDataColumn( input_output::TrackingDataType::utc_reception_time_j2000 ) );
          observationTimes_ = convertTimesTdbFromJ2000( observation_models::receiver );
          break;
      }
      case calendar_day_time:
      {
          // Convert dates to Julian days since J2000
          auto years = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::year);
          auto months = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::month);
          auto days = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::day);
          auto hours = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::hour);
          auto minutes = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::minute);
          auto seconds = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::second);
          // Convert to seconds and add to utc times
          observationTimesUtc_.clear( );
          for( unsigned int i = 0; i < years.size( ); i++ )
          {
              observationTimesUtc_.push_back(basic_astrodynamics::DateTime(
                  years.at( i ), months.at( i ), days.at( i ), hours.at( i ), minutes.at( i ), seconds.at( i ) ).epoch< TimeType >() );
          }
          // Convert to TDB
          observationTimes_ = convertTimesTdbFromJ2000( observation_models::receiver );

          break;
      }
          // Throw error if representation not implemented
      default:
      {
          throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
      }
      }

      // Get the delays in the time tag (or set to 0.0 if not specified)
      std::vector<double> timeTagDelays = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::time_tag_delay, 0.0);
      for (size_t idx = 0; idx < observationTimes_.size(); ++idx) {
          observationTimes_[idx] -= timeTagDelays[idx];
      }
  }

  //! Method to update the observations
  void updateObservations()
  {
      // Update the observableTypes that one can expect to process
      observationMap_.clear();

      for (const ObservableType observableType : observableTypes_)
      {
          std::vector< ObservationScalarType > observableValues;

          // Convert the raw data to required observables
          switch (observableType)
          {
          case n_way_range:
          {
              // Conversion function for n-way range
              auto lightTimeRangeConversion =
                  [](ObservationScalarType lightTime, ObservationScalarType lightTimeDelay) {
                  return (lightTime - lightTimeDelay) * physical_constants::getSpeedOfLight< ObservationScalarType >(); };

              // Extract columns from raw data and convert to observable values
              std::vector< ObservationScalarType > lightTimes =
                  utilities::staticCastVector< ObservationScalarType, double >( rawTrackingTxtFileContents_->getDoubleDataColumn(
                  input_output::TrackingDataType::n_way_light_time) );
              std::vector< ObservationScalarType > lightTimeDelays =
                  utilities::staticCastVector< ObservationScalarType, double >( rawTrackingTxtFileContents_->getDoubleDataColumn(
                  input_output::TrackingDataType::light_time_measurement_delay, 0.0 ) );
              observableValues = utilities::convertVectors< ObservationScalarType >(lightTimeRangeConversion, lightTimes, lightTimeDelays);
              break;
          }
          case doppler_measured_frequency:
          {
              // Conversion function for doppler measured frequency
              auto dopplerFrequencyConversion = [](double dopplerFrequency, double dopplerBaseFrequency) {
                  return dopplerFrequency + dopplerBaseFrequency; };

              // Extract columns from raw data and convert to observable values
              std::vector<double> dopplerFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_measured_frequency);
              std::vector<double> dopplerBaseFrequencies = rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_base_frequency);

              // Check if any of the base frequencies are zero
              if (std::any_of(dopplerBaseFrequencies.begin(), dopplerBaseFrequencies.end(), [](double baseFrequency) { return baseFrequency == 0.0; }))
              {
                  std::cerr << "Warning when processing doppler_measured_frequency. Doppler base frequency is zero." << std::endl;
              }

              observableValues = utilities::convertVectors< ObservationScalarType >(dopplerFrequencyConversion, dopplerFrequencies, dopplerBaseFrequencies);
              break;
          }
          case dsn_n_way_averaged_doppler:
          {
              observableValues = utilities::staticCastVector< ObservationScalarType, double >(
                  rawTrackingTxtFileContents_->getDoubleDataColumn(input_output::TrackingDataType::doppler_averaged_frequency ) );
              break;
          }
          default: {
              throw std::runtime_error("Error while processing tracking txt file. ObservableType conversion not implemented");
          }
          }

          // Store observables
          observationMap_[observableType] = observableValues;
      }
  }


    //! Update the available observableTypes described by the processed data file
  void updateAvailableObservableTypes()
  {
        observableTypes_ = findAvailableObservableTypes(rawTrackingTxtFileContents_->getDataColumnTypes());
  }

  /*!
   * Convert the UTC observation times to TDB
   * @return vector of TDB observation times
   */
  std::vector< TimeType > convertTimesTdbFromJ2000( const LinkEndType referenceLinkEnd )
  {
      // Get the timescale converter

      // Check if there is one LinkEnds per observation time
      if (linkEndsVector_.size() != observationTimesUtc_.size())
      {
          throw std::runtime_error("Error while processing tracking data: vector of linkEnds and observationTimes not of equal size");
      }

      // Ge a vector of ground station positions
      std::vector<Eigen::Vector3d> groundStationPositions;
      for (const auto& linkEnds : linkEndsVector_)
      {
          try
          {
              std::string currentGroundStation = linkEnds.at( referenceLinkEnd ).getStationName( );
              groundStationPositions.push_back( earthFixedGroundStationPositions_.at( currentGroundStation ));
          }
          catch( const std::runtime_error& error )
          {
              throw std::runtime_error( "Error when creating getting link ends for time conversion in tracking file processing: "
                                        + std::string( error.what( ) ) + ". For reference link end " + std::to_string( static_cast< int >( referenceLinkEnd ) ) );
          }
      }

      // Convert to TDB using the GS positions
      std::vector< TimeType > observationTimesTdb =
          timeScaleConverter_->getCurrentTimes< TimeType >(
          basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, observationTimesUtc_, groundStationPositions );
      return observationTimesTdb;
  }

  //! Utility function to get the ground station id
  //! TODO: This is temporary until the functionality to read from file is implemented
  static std::string getStationNameFromStationId(const int stationId, const std::string networkPrefix="DSS-")
  {
    return networkPrefix + std::to_string(stationId);
  }

  void updateAncilliarySettings( const ObservableType observableType,
                                 std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySettings )
  {
      if (!utilities::containsAll(observableTypes_, std::vector< ObservableType >( { observableType } ) ) )
      {
          throw std::runtime_error( "Error when getting ancilliary settings from processed file contents, could not find " +
          getObservableName( observableType ) );
      }

      switch( observableType )
      {
      case dsn_n_way_averaged_doppler:
      {
          ancilliarySettings->setAncilliaryDoubleData( doppler_integration_time, getObservationTimeStep( ) );
      }
      default:
          break;
      }
  }

// Settings interpreting the file format.
private:

  //! Enum of implemented time representations. Combination of columns in the processed file contents
  enum TimeRepresentation
  {
    calendar_day_time,
    tdb_seconds_j2000,
    utc_seconds_j2000
  };

  //! Enum of implemented link ends representations. Combination of columns in the processed file contents
  enum LinkEndsRepresentation
  {
    dsn_transmitting_receiving_station_nr,
    transmitting_receiving_station_name
  };

  //! Deduce the time representation from available columns
  TimeRepresentation getTimeRepresentation()
  {

      // Get all the data types from the raw file contents
      auto const& availableDataTypes = rawTrackingTxtFileContents_->getDataColumnTypes();

      // Return representation based on available data types
      if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::utc_reception_time_j2000}))
      {
          return utc_seconds_j2000;
      }
      if (utilities::containsAll(availableDataTypes,
                                 std::vector<input_output::TrackingDataType>{
                                     input_output::TrackingDataType::year,
                                     input_output::TrackingDataType::month,
                                     input_output::TrackingDataType::day,
                                     input_output::TrackingDataType::hour,
                                     input_output::TrackingDataType::minute,
                                     input_output::TrackingDataType::second }))
      {
          return calendar_day_time;
      }
      if (utilities::containsAll(availableDataTypes, std::vector<input_output::TrackingDataType>{input_output::TrackingDataType::tdb_reception_time_j2000}))
      {
          return tdb_seconds_j2000;
      }

      // Throw an error if no match is found
      throw std::runtime_error("Error while processing tracking txt file: Time representation not recognised or implemented.");
  }

  //! Deduce the link ends representation from available columns
  LinkEndsRepresentation getLinkEndsRepresentation()
  {
      // Get all the available data columns
      auto const& availableDataTypes = rawTrackingTxtFileContents_->getAllAvailableDataTypes();

      // Porvide link Ends representation based on available columns

      if (utilities::containsAll(availableDataTypes,
                                 std::vector<input_output::TrackingDataType>{
                                     input_output::TrackingDataType::dsn_transmitting_station_nr,
                                     input_output::TrackingDataType::dsn_receiving_station_nr } ) )
      {
          return dsn_transmitting_receiving_station_nr;
      }

      if (utilities::containsAll(availableDataTypes,
                                 std::vector<input_output::TrackingDataType>{
                                     input_output::TrackingDataType::transmitting_station_name,
                                     input_output::TrackingDataType::receiving_station_name } ) )
      {
          return transmitting_receiving_station_name;
      }

      // Throw error if no match is found
      throw std::runtime_error("Error while processing tracking txt file: Link Ends representation not recognised or implemented.");
  }


// Getters
public:
  bool is_initialised() const
  {
    return initialised_;
  }

  const std::vector<ObservableType>& getObservableTypes() const
  {
    return observableTypes_;
  }

  const std::vector< TimeType >& getObservationTimes() const
  {
    return observationTimes_;
  }

  const std::vector<LinkEnds>& getLinkEndsVector() const
  {
    return linkEndsVector_;
  }

  const std::set<LinkEnds>& getLinkEndsSet() const
  {
    return linkEndsSet_;
  }

  int getNumRows() const
  {
    return rawTrackingTxtFileContents_->getNumRows();
  }

  const std::map<ObservableType, std::vector<ObservationScalarType>> getObservationMap() const
  {
    return observationMap_;
  }

    std::shared_ptr<input_output::TrackingTxtFileContents> getRawTrackingTxtFileContents( )
    {
      return rawTrackingTxtFileContents_;
    }


private:

  double getObservationTimeStep( )
  {
      if( observationTimesUtc_.size( ) < 2 )
      {
          throw std::runtime_error( "Error when getting integration time for processed file contents, size is < 2" );
      }
      double observationTimeStep = observationTimesUtc_.at( 1 ) - observationTimesUtc_.at( 0 );
      for( unsigned int i = 1; i < observationTimesUtc_.size( ); i++ )
      {
          double testObservationTimeStep = observationTimesUtc_.at( i ) - observationTimesUtc_.at( i - 1 );
          if( std::fabs( observationTimeStep - testObservationTimeStep ) > 0.01 )// 50.0 * std::numeric_limits< double >::epsilon( ) * observationTimesUtc_.at( i - 1 )  )
          {
              std::cout<<std::setprecision( 19 )<<i<<" "<<observationTimesUtc_.at( i )<<" "<<observationTimesUtc_.at( i - 1 )<<" "<<observationTimesUtc_.at( i - 2 )<<" "<<testObservationTimeStep<<" "<<observationTimeStep<<" "<<testObservationTimeStep - observationTimeStep<<" "<<
                       50.0 * std::numeric_limits< double >::epsilon( ) * observationTimesUtc_.at( i - 1 )<<std::endl;
              throw std::runtime_error( "Error when getting integration time for processed file contents, step is not equal" );
          }
      }

      return observationTimeStep;
  }

  //! TrackingTxtFileContents raw file contents
  std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents_;

  //! Name of the spacecraft (as chosen by the user)
  std::string spacecraftName_;

  //! Vector of TDB observation times
  std::vector<TimeType> observationTimes_;

  std::vector<TimeType> observationTimesUtc_;

    //! Map of observations
  std::map<ObservableType, std::vector<ObservationScalarType>> observationMap_;

  //! Vector of observable types that can be derived from the tracking data
  std::vector<ObservableType> observableTypes_;

  //! Vector of linkEnds for every data type
  std::vector<LinkEnds> linkEndsVector_;

  //! Set of distinct linkEnds
  std::set<LinkEnds> linkEndsSet_;

  //! Map with earth fixed positions for all ground stations
  std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions_;

  //! Boolean verifying that the initialisation is complete
  bool initialised_ = false;

  //! Create link ends representation
  std::map<LinkEndType, input_output::TrackingDataType> linkEndsRepresentation_;

  std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter_ = earth_orientation::createDefaultTimeConverter( );

};


/*!
 * Function to create an observation collection from the processed Tracking file data
 * @param processedTrackingTxtFileContents
 * @param observableTypesToProcess
 * @param ancillarySettings
 * @return observation collection
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
    createTrackingTxtFileObservationSets(
    const std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings() )
{

    // Make sure processing the tracking file was successful
    if (!processedTrackingTxtFileContents->is_initialised())
    {
        throw std::runtime_error("Error while processing tracking txt file: processedTrackingTxtFileContents was never initialised.");
    }

    // Get the double map and set of distinct LinkEnds
    const auto& allObservationsMap = processedTrackingTxtFileContents->getObservationMap();
    std::set<LinkEnds> linkEndsToProcess = processedTrackingTxtFileContents->getLinkEndsSet();

    // Check observable types to process. If empty, process all available, if impossible, throw error
    std::vector<ObservableType> availableObservableTypes = processedTrackingTxtFileContents->getObservableTypes();
    if (observableTypesToProcess.empty())
    {
        observableTypesToProcess = availableObservableTypes;
    }
    if (!utilities::containsAll(availableObservableTypes, observableTypesToProcess))
    {
        throw std::runtime_error("Error while processing Tracking txt file. Not enough information to extract requested observables");
    }

    if( observableTypesToProcess.size( ) > 1 )
    {
        throw std::runtime_error( "Error, can only process one observable type at a time from text file." );
    }


    // Initialise necessary maps
    std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > > observationSets;
    std::map<ObservableType, std::map<LinkEnds, std::vector<TimeType>>> observationTimesMap;
    std::map<ObservableType, std::map<LinkEnds, std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >>> observablesMap;

    // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
    std::vector<TimeType> allObservationTimes = processedTrackingTxtFileContents->getObservationTimes() ;
    std::vector<LinkEnds> linkEndsVector = processedTrackingTxtFileContents->getLinkEndsVector();
    std::set<LinkEnds> linkEndsSet = processedTrackingTxtFileContents->getLinkEndsSet();

    std::shared_ptr< ObservationAncilliarySimulationSettings > updatedAncilliarySettings =
        std::make_shared<ObservationAncilliarySimulationSettings>( ancillarySettings );
    processedTrackingTxtFileContents->updateAncilliarySettings( availableObservableTypes.at( 0 ), updatedAncilliarySettings );

    // Prepare maps that order all observations per observable type and link ends
    // This is necessary for files where the linkends are not always the same
    // Loop over all observation times (rows in the file)
    for (size_t i = 0; i < allObservationTimes.size(); ++i)
    {
        // Get link ends of this observation
        LinkEnds& currentLinkEnds = linkEndsVector[i];

        // Loop over the available observable types
        for (const ObservableType& currentObservableType : observableTypesToProcess)
        {

            // Prepare a container to store the data columns
            Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> currentObservable(1);
            currentObservable(0) = allObservationsMap.at(currentObservableType).at( i );

            // Store in correct maps
            observationTimesMap[currentObservableType][currentLinkEnds].push_back(allObservationTimes.at( i ));
            observablesMap[currentObservableType][currentLinkEnds].push_back(currentObservable);
        }
    }

    // Fill the observation collection
    for (ObservableType& currentObservableType : observableTypesToProcess)
    {
        for (const LinkEnds& currentLinkEnds : linkEndsSet)
        {
            observationSets[currentObservableType][currentLinkEnds].push_back(
                std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
                    currentObservableType,
                    currentLinkEnds,
                    observablesMap[currentObservableType][currentLinkEnds],
                    observationTimesMap[currentObservableType][currentLinkEnds],
                    receiver, // TODO: make more flexible to allow for other reference link ends
                    std::vector<Eigen::VectorXd>(),
                    nullptr, updatedAncilliarySettings ));
        }
    }

    return observationSets;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFilesObservationCollection(
    const std::vector< std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > > processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings())
{
    std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > > observationSets;

    for( unsigned int i = 0; i < processedTrackingTxtFileContents.size( ); i++ )
    {
        std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
            processedObervationSet = createTrackingTxtFileObservationSets< ObservationScalarType, TimeType >(
            processedTrackingTxtFileContents.at( i ), observableTypesToProcess, ancillarySettings );
        for( auto it : processedObervationSet )
        {
            for( auto it2 : it.second )
            {
                for( unsigned int j = 0; j < it2.second.size( ); j++ )
                {
                    observationSets[ it.first ][ it2.first ].push_back( it2.second.at( j ) );
                }
            }
        }
    }
    return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >( observationSets );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    const std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings())
{
    return createTrackingTxtFilesObservationCollection< ObservationScalarType, TimeType >(
        { processedTrackingTxtFileContents }, observableTypesToProcess, ancillarySettings );
}

/*!
 * Function to create an observation collection from the raw Tracking file data
 * @param rawTrackingTxtFileContents The raw tracking file contents
 * @param spacecraftName Name of the spacecraft
 * @param observableTypesToProcess Vector of observable types that need to be in the collection
 * @param earthFixedGroundStationPositions map of ground station positions
 * @param ancillarySettings
 * @return observation collection
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    const std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents,
    const std::string spacecraftName,
    const std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const std::map<std::string, Eigen::Vector3d> earthFixedGroundStationPositions = simulation_setup::getCombinedApproximateGroundStationPositions(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings())
{
    // Create processed tracking file contents
    auto processedTrackingTxtFileContents = std::make_shared<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >(
        rawTrackingTxtFileContents, spacecraftName, earthFixedGroundStationPositions );

    // Create observation collection and return
    return createTrackingTxtFileObservationCollection(processedTrackingTxtFileContents,
                                                    observableTypesToProcess,
                                                    ancillarySettings);
}

void setStationFrequenciesFromTrackingData(
    const std::map< std::string, std::vector< std::tuple< std::vector< double >, std::vector< double >, std::vector< double > > > >& rampInformation,
    simulation_setup::SystemOfBodies& bodies );

template< typename ObservationScalarType = double, typename TimeType = double >
void setStationFrequenciesFromTrackingData(
    const std::vector< std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > > processedTrackingTxtFileContents,
    simulation_setup::SystemOfBodies& bodies )
{
    std::map< std::string, std::vector< std::tuple< std::vector< double >, std::vector< double >, std::vector< double > > > > rampInformation;
    for( unsigned int i = 0; i < processedTrackingTxtFileContents.size( ); i++ )
    {
        std::shared_ptr<input_output::TrackingTxtFileContents> fileContents = processedTrackingTxtFileContents.at( i )->getRawTrackingTxtFileContents( );
        std::vector< double > rampUtcTimes = fileContents->getDoubleDataColumn( input_output::TrackingDataType::utc_ramp_referencee_j2000 );
        std::vector< double > frequencyRampRates = fileContents->getDoubleDataColumn( input_output::TrackingDataType::transmission_frequency_linear_term );
        std::vector< double > frequencyValues = fileContents->getDoubleDataColumn( input_output::TrackingDataType::transmission_frequency_constant_term );

        std::string transmitterName;
        if( processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).size( ) != 1 )
        {
            throw std::runtime_error( "Error when getting link ends from IFMS file, found multiple link ends sets." +
                                      std::to_string( processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).size( ) ) );
        }
        else
        {
            LinkEnds currentLinkEnds = *(processedTrackingTxtFileContents.at( i )->getLinkEndsSet( ).begin( ) );
            transmitterName = currentLinkEnds.at( transmitter ).stationName_;
        }
        rampInformation[ transmitterName ].push_back( std::make_tuple( rampUtcTimes, frequencyValues, frequencyRampRates ) );
    }
    setStationFrequenciesFromTrackingData( rampInformation, bodies );
}

template< typename ObservationScalarType = double, typename TimeType = double >
void setTrackingDataInformationInBodies(
    const std::vector< std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > > processedTrackingTxtFileContents,
    simulation_setup::SystemOfBodies& bodies,
    ObservableType observableType )
{
    for( unsigned int i = 0; i < processedTrackingTxtFileContents.size( ); i++ )
    {
        std::vector<ObservableType> availableObservableTypes = processedTrackingTxtFileContents.at( i )->getObservableTypes( );
        if ( !utilities::containsAll( availableObservableTypes, std::vector< ObservableType >( { observableType } ) ) )
        {
            throw std::runtime_error(
                "Error while processing Tracking txt file for body properties. Observable not found: " + getObservableName( observableType ) );
        }
    }


    switch( observableType )
    {
    case dsn_n_way_averaged_doppler:
    {
        setStationFrequenciesFromTrackingData( processedTrackingTxtFileContents, bodies );
    }
    default:
        break;
    }

}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createMultiStationIfmsObservedObservationCollectionFromFiles(
    const std::vector< std::string >& ifmsFileNames,
    simulation_setup::SystemOfBodies& bodies,
    const std::string& targetName,
    const std::vector< std::string >& groundStationNames,
    const FrequencyBands& receptionBand,
    const FrequencyBands& transmissionBand,
    const bool applyTroposphereCorrection = true,
    const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
    simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    if( groundStationNames.size( ) != ifmsFileNames.size( ) )
    {
        throw std::runtime_error( "Error when loading IFMS files, list of files has different size than list of stations" );
    }
    std::vector< std::shared_ptr< input_output::TrackingTxtFileContents > > rawIfmsDataList;

    for( std::string ifmsFileName : ifmsFileNames )
    {
        rawIfmsDataList.push_back( input_output::readIfmsFile( ifmsFileName, applyTroposphereCorrection ) );
    }

    std::vector< std::shared_ptr< ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > > processedIfmsFiles;
    for( unsigned int i = 0; i < rawIfmsDataList.size( ); i++ )
    {
        rawIfmsDataList.at( i )->addMetaData( input_output::TrackingDataType::receiving_station_name, groundStationNames.at( i ) );
        rawIfmsDataList.at( i )->addMetaData( input_output::TrackingDataType::transmitting_station_name, groundStationNames.at( i ) );
        processedIfmsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >(
            rawIfmsDataList.at( i ), targetName, earthFixedGroundStationPositions ) );
    }

    setTrackingDataInformationInBodies(
        processedIfmsFiles, bodies, dsn_n_way_averaged_doppler );

    ObservationAncilliarySimulationSettings ancilliarySettings;
    ancilliarySettings.setAncilliaryDoubleVectorData(frequency_bands, { static_cast< double >( transmissionBand ), static_cast< double >( receptionBand ) });
    ancilliarySettings.setAncilliaryDoubleData( doppler_reference_frequency, 0.0 );
    ancilliarySettings.setAncilliaryDoubleData( reception_reference_frequency_band, convertFrequencyBandToDouble( receptionBand ) );

    return observation_models::createTrackingTxtFilesObservationCollection< ObservationScalarType, TimeType >(
        processedIfmsFiles, std::vector<ObservableType>(), ancilliarySettings );
}


template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createIfmsObservedObservationCollectionFromFiles(
    const std::vector< std::string >& ifmsFileNames,
    simulation_setup::SystemOfBodies& bodies,
    const std::string& targetName,
    const std::string& groundStationName,
    const FrequencyBands& receptionBand,
    const FrequencyBands& transmissionBand,
    const bool applyTroposphereCorrection = true,
    const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
    simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    std::vector< std::string > groundStationNameList( ifmsFileNames.size( ), groundStationName);

    return createMultiStationIfmsObservedObservationCollectionFromFiles< ObservationScalarType, TimeType >(
        ifmsFileNames, bodies, targetName, groundStationNameList, receptionBand, transmissionBand, applyTroposphereCorrection, earthFixedGroundStationPositions );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createFdetsObservedObservationCollectionFromFile(
    const std::string& fdetsFileName,
    const double& baseFrequency,
    const std::vector<std::string>& columnTypes,
    const std::string& targetName,
    const std::string& transmittingStationName,
    const std::string& receivingStationName,
    const FrequencyBands& receptionBand,
    const FrequencyBands& transmissionBand,
    const std::map< std::string, Eigen::Vector3d >& earthFixedGroundStationPositions =
    simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    using namespace input_output;
    std::shared_ptr<TrackingTxtFileContents> fdetsFileContents = readFdetsFile( fdetsFileName, columnTypes );
    fdetsFileContents->addMetaData( TrackingDataType::receiving_station_name, receivingStationName );
    fdetsFileContents->addMetaData( TrackingDataType::transmitting_station_name, transmittingStationName );
    fdetsFileContents->addMetaData( TrackingDataType::doppler_base_frequency, baseFrequency);

    std::vector< std::shared_ptr< ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > > > processedFdetsFiles;
    processedFdetsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents< ObservationScalarType, TimeType > >(
        fdetsFileContents, targetName, earthFixedGroundStationPositions ) );

    // Define ancilliary settings
    ObservationAncilliarySimulationSettings ancilliarySettings;
    ancilliarySettings.setAncilliaryDoubleVectorData(frequency_bands, { static_cast< double >( transmissionBand ), static_cast< double >( receptionBand ) });

    return observation_models::createTrackingTxtFilesObservationCollection< ObservationScalarType, TimeType >(
        processedFdetsFiles, std::vector<ObservableType>( ), ancilliarySettings );
}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
