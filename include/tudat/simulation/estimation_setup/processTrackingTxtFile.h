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
#include "tudat/simulation/estimation_setup/observations.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"

// TODO:
//  Ancillary settings (frequency, time scale, ...)
//  Format checks
//  Unit Checks
//  Add time conversions



namespace tudat
{
namespace observation_models
{

//! Utility function that might belong elsewhere
template< typename TimeScalarType = double >
TimeScalarType convertCalendarDateToJulianDaySinceJ2000(const int calendarYear,
                                                        const int calendarMonth,
                                                        const int calendarDay,
                                                        const int calendarHour,
                                                        const int calendarMinutes,
                                                        const TimeScalarType calendarSeconds)
{
  return basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch<TimeScalarType>
      (calendarYear, calendarMonth, calendarDay, calendarHour, calendarMinutes,
       calendarSeconds, basic_astrodynamics::JULIAN_DAY_ON_J2000);
}

//! Utility function to check if a container contains all elements of another container
template< typename T, typename U >
bool containsAll(const T& referenceArray, const U searchArray)
{
  return std::includes(referenceArray.begin(), referenceArray.end(), searchArray.begin(), searchArray.end());
}

//! Utility function to apply the convert function to vectors of arguments
template< typename ConvertFunc, typename FirstArg, typename... Args >
std::vector<double> convertVectors(ConvertFunc convertFunc, const std::vector<FirstArg>& firstArg, const std::vector<Args>& ... args)
{
  std::vector<double> result;
  size_t N = firstArg.size(); // Assuming all vectors have the same size
  for (size_t i = 0; i < N; ++i) {
    result.push_back(convertFunc(firstArg[i], args[i]...));
  }
  return result;
}

//! Utility function to convert a vector to a set
template< typename T >
std::set<T> vectorToSet(std::vector<T> vector)
{
  std::set<T> set;
  for (T& elem : vector) {
    set.insert(elem);
  }
  return set;
}

//! Map containing all the tracking data types that correspond to an observable. If a specific set of TrackingDataTypes (keys)
//! Is present in the file, the corresponding ObservableType (value) can later be added to the  ObservationCollection
static const std::map<std::vector<input_output::TrackingDataType>, ObservableType> dataTypeToObservableMap = {
    {{input_output::TrackingDataType::n_way_light_time}, n_way_range},
    // Todo: More to be implemented
};

//! Function to extract a vector of all the observableTypes that can be deduced from a list of DataTypes
std::vector<ObservableType> findAvailableObservableTypes(std::vector<input_output::TrackingDataType> availableDataTypes);

//! Function to create ancillary settings
template< typename TimeType = double >
ObservationAncilliarySimulationSettings createTrackingAncillarySettings(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
{
  // TODO: Currently not used
  return ObservationAncilliarySimulationSettings();
}

//! Class containing the processed file contents for a Tracking txt data file
class ProcessedTrackingTxtFileContents
{
//  Constructors
public:
  ProcessedTrackingTxtFileContents(std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtContents, std::string spacecraftName)
      : rawTrackingTxtFileContents_(rawTrackingTxtContents), spacecraftName_(spacecraftName)
  {
    initialise();
  }

  void initialise()
  {
    updateObservationTimes();
    updateObservations();
    updateLinkEnds();
    updateObservableTypes();
    initialised_ = true;
  }


// Update  information from the raw file data
public:

  //! Method to extract the first and last time of the observations
  std::pair<double, double> getStartAndEndTime()
  {
    const std::vector<double>& times = getObservationTimes();
    std::pair<double, double> startEndTime({times.front(), times.back()});
    return startEndTime;
  }

  //! Method to (re)calculate the times from the Tracking file
  void updateObservationTimes();

  //! Method to create a vector of the linkEnds for each of the lines in the file (different ground stations might be used)
  void updateLinkEnds();

  // TODO: It would be better if the observables are already extracted here, instead of in the separate function
  void updateObservations();

  void updateObservableTypes()
  {
    observableTypes_ = findAvailableObservableTypes(rawTrackingTxtFileContents_->getDataColumnTypes());
  }

  //! Utility function to get the ground station id
  std::string getStationNameFromStationId(const std::string networkPrefix, const int stationId)
  {
    return networkPrefix + std::to_string(stationId);
  }

// Settings interpreting the file format.
private:

  enum TimeRepresentation
  {
    day_time,
  };

  enum LinkEndsRepresentation
  {
    dsn_transmitting_receiving_station_nr,
  };

  TimeRepresentation getTimeRepresentation();

  LinkEndsRepresentation getLinkEndsRepresentation();


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

  const std::vector<double>& getObservationTimes() const
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

  // TODO:
  //  Instead of this, process the data in the updateObservations functions!
  const std::map<input_output::TrackingDataType, std::vector<double>>& getDoubleDataMap()
  {
    return rawTrackingTxtFileContents_->getDoubleDataMap();
  }

private:
  std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents_;
  std::string spacecraftName_;
  std::vector<double> observationTimes_;
  std::vector<ObservableType> observableTypes_;
  std::vector<LinkEnds> linkEndsVector_;
  std::set<LinkEnds> linkEndsSet_;
  bool initialised_;
};

//! Function to create an observation collection from the processed Tracking file data
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<observation_models::ProcessedTrackingTxtFileContents> processedTrackingTxtFileContents,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{

  if (!processedTrackingTxtFileContents->is_initialised()) {
    throw std::runtime_error("Error while processing tracking txt file: processedTrackingTxtFileContents was never initialised.");
  }

  const auto& dataMap = processedTrackingTxtFileContents->getDoubleDataMap();

  // Check expected sizes
  if (observableTypesToProcess.empty()) {
    observableTypesToProcess = processedTrackingTxtFileContents->getObservableTypes();
  }
  std::set<LinkEnds> linkEndsToProcess = processedTrackingTxtFileContents->getLinkEndsSet();

  std::map<ObservableType, std::map<LinkEnds, std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
      observationSets;

  std::map<ObservableType, std::map<LinkEnds, std::vector<TimeType>>> observationTimesMap;
  std::map<ObservableType, std::map<LinkEnds, std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> >>> observablesMap;

  // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
  std::vector<TimeType> allObservationTimes = processedTrackingTxtFileContents->getObservationTimes();
  std::vector<LinkEnds> linkEndsVector = processedTrackingTxtFileContents->getLinkEndsVector();
  std::set<LinkEnds> linkEndsSet = processedTrackingTxtFileContents->getLinkEndsSet();

  for (size_t i = 0; i < allObservationTimes.size(); ++i) {
    LinkEnds& currentLinkEnds = linkEndsVector[i];

    for (auto& pair : dataTypeToObservableMap) {
      std::vector<input_output::TrackingDataType> currentTrackingDataTypes = pair.first;
      const ObservableType& currentObservableType = pair.second;

      if (containsAll(observableTypesToProcess, std::vector<ObservableType>{pair.second})) {
        size_t currentSize = currentTrackingDataTypes.size();
        Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> currentObservable(currentSize);

        for (size_t j = 0; j < currentSize; ++j) {
          currentObservable[j] = dataMap.at(currentTrackingDataTypes[j])[i];
        }
        observationTimesMap[currentObservableType][currentLinkEnds].push_back(allObservationTimes[i]);
        observablesMap[currentObservableType][currentLinkEnds].push_back(currentObservable);
      }
    }
  }

  for (ObservableType& currentObservableType : observableTypesToProcess) {
    for (const LinkEnds& currentLinkEnds : linkEndsSet) {
      observationSets[currentObservableType][currentLinkEnds].push_back(std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
          currentObservableType,
          currentLinkEnds,
          observablesMap[currentObservableType][currentLinkEnds],
          observationTimesMap[currentObservableType][currentLinkEnds],
          receiver, // Todo: Not sure how this is used (copied from odfFiles)
          std::vector<Eigen::VectorXd>(),
          nullptr,
          std::make_shared<ObservationAncilliarySimulationSettings>(ancillarySettings)));
    }
  }

  return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
}

//! Function to create an observation collection from the raw Tracking file data
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<observation_models::ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<input_output::TrackingTxtFileContents> rawTrackingTxtFileContents,
    std::string spacecraftName,
    std::vector<ObservableType> observableTypesToProcess = std::vector<ObservableType>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{
  auto processedTrackingTxtFileContents = std::make_shared<observation_models::ProcessedTrackingTxtFileContents>(rawTrackingTxtFileContents, spacecraftName);
  return createTrackingTxtFileObservationCollection(processedTrackingTxtFileContents, observableTypesToProcess, ancillarySettings, startAndEndTimesToProcess);
}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSTRACKINGTXTFILE_H
