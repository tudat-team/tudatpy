/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
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
//#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/basic_astro/timeConversions.h"
//#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/simulation/estimation_setup/observations.h"
//#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
//#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/interpolators/lookupScheme.h"
//#include "tudat/math/quadrature/trapezoidQuadrature.h"

namespace tudat
{
namespace observation_models
{

// This map specifies what fields are expected in the TrackingTxtFileContents to extract a specific observable
std::map<ObservableType, std::vector<input_output::TrackingDataType>> requiredDataTypesForObservable = {
    {one_way_range, {}},
    {angular_position, {}},
    {position_observable, {}},
    {one_way_doppler, {}},
    {one_way_differenced_range, {}},
    {n_way_range, {}},
    {two_way_doppler, {}},
    {euler_angle_313_observable, {}},
    {velocity_observable, {}},
    {relative_angular_position, {}},
    {n_way_differenced_range, {}},
    {relative_position_observable, {}},
    {dsn_one_way_averaged_doppler, {}},
    {dsn_n_way_averaged_doppler, {}},
};

// Function to check if a container contains all elements of another container
template< typename T, typename U >
bool includesAll(const T& referenceArray, const U& searchArray)
{
  return std::includes(referenceArray.begin(), referenceArray.end(), searchArray.begin(), searchArray.end());
}

// Create ancillary settings
template< typename TimeType = double >
ObservationAncilliarySimulationSettings createOdfAncillarySettings(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents)
{
  // Currently not implemented
  return ObservationAncilliarySimulationSettings();
}



// Create observation collection
// TODO: This assumes the ancillary settings are the same for all the observables. For txt files this will usually be the case
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr<ObservationCollection<ObservationScalarType, TimeType> >
createTrackingTxtFileObservationCollection(
    std::shared_ptr<input_output::TrackingTxtFileContents> trackingTxtFileContents,
    std::map<ObservableType, LinkEnds> observableTypesToProcess = std::map<ObservableType, LinkEnds>(),
    const ObservationAncilliarySimulationSettings& ancillarySettings = ObservationAncilliarySimulationSettings(),
    std::pair<TimeType, TimeType> startAndEndTimesToProcess = std::make_pair<TimeType, TimeType>(TUDAT_NAN, TUDAT_NAN))
{
  const auto& availableDataTypes = trackingTxtFileContents->getDataColumnTypes();
  const auto& dataMap = trackingTxtFileContents->getDoubleDataMap();
  const auto& numDataRows = trackingTxtFileContents->getNumRows();

  // Check expected sizes
  if (observableTypesToProcess.empty()) {
    throw std::runtime_error("Error while processing tracking txt file: Empty list of observableTypes to process");
  }

  // Create and fill single observation sets
  std::map<ObservableType,
           std::map<LinkEnds,
                    std::vector<std::shared_ptr<SingleObservationSet<ObservationScalarType, TimeType> > > > >
      observationSets;

  // C++ 17 allows looping over two arrays at once.
  for (auto& pair : observableTypesToProcess) {
    ObservableType currentObservableType = pair.first;
    LinkEnds currentLinkEnds = pair.second;

    // Check if the file contains the necessary data types.
    const auto& requiredDataTypes = requiredDataTypesForObservable[currentObservableType];
    if (!includesAll(availableDataTypes, requiredDataTypes)) {
      throw std::runtime_error(
          "Error while processing tracking txt file: Not all requiredDataTypes are present in the provided Tracking file contents");
    }

    // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
    std::vector<std::vector<TimeType> > observationTimes;
    std::vector<std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > > observables;

    // TODO: Do this conversion in a separate function (including more options for how the time is stored)
    auto years = dataMap.at(input_output::TrackingDataType::year);
    auto months = dataMap.at(input_output::TrackingDataType::month);
    auto days = dataMap.at(input_output::TrackingDataType::day);
    auto hours = dataMap.at(input_output::TrackingDataType::hour);
    auto minutes = dataMap.at(input_output::TrackingDataType::minute);
    auto seconds = dataMap.at(input_output::TrackingDataType::second);
    for (size_t i; i < numDataRows; ++i) {
      observationTimes.push_back(tudat::basic_astrodynamics::convertCalendarDateToJulianDay(years[i],
                                                                                            months[i],
                                                                                            days[i],
                                                                                            hours[i],
                                                                                            minutes[i],
                                                                                            seconds[i]));
    }

    // FIXME: Fill the values and convert the times


    // Look at separateSingleLinkOdfData;


    // Todo: Maybe, truncating the observations might be useful. Currently not implemented, because most files until now relatively small


    // Create the single observation sets and save them
    for (unsigned int i = 0; i < observationTimes.size(); ++i) {
      observationSets[currentObservableType][currentLinkEnds].push_back(
          std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
              currentObservableType, currentLinkEnds, observables.at(i), observationTimes.at(i),
              receiver,
              std::vector<Eigen::VectorXd>(),
              nullptr, std::make_shared<ObservationAncilliarySimulationSettings>(
                  ancillarySettings)));
    }

    return std::make_shared<observation_models::ObservationCollection<ObservationScalarType, TimeType> >(observationSets);
  }
}

// TODO: Add a variation to the function where two vectors are provided instead of a map!





//
//  for (auto observableTypeIterator = trackingTxtFileContents->getProcessedDataBlocks().begin();
//       observableTypeIterator != trackingTxtFileContents->getProcessedDataBlocks().end(); ++observableTypeIterator) {
//    ObservableType currentObservableType = observableTypeIterator->first;
//
//
//
//    // Check if an observation set should be created for the current observable type
//    if (std::count(observableTypesToProcess.begin(), observableTypesToProcess.end(), currentObservableType) == 0) {
//      continue;
//    }
//
//    for (auto linkEndsIterator = observableTypeIterator->second.begin();
//         linkEndsIterator != observableTypeIterator->second.end(); ++linkEndsIterator) {
//      LinkEnds currentLinkEnds = linkEndsIterator->first;
//      std::shared_ptr<ProcessedOdfFileSingleLinkData> currentOdfSingleLinkData = linkEndsIterator->second;
//
//      // Get vectors of times, observations, and ancillary settings for the current observable type and link ends
//      std::vector<std::vector<TimeType> > observationTimes;
//      std::vector<std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > > observables;
//      std::vector<ObservationAncilliarySimulationSettings> ancillarySettings;
//
//      // Fill vectors
//      separateSingleLinkOdfData(
//          currentObservableType, currentOdfSingleLinkData, observationTimes, observables, ancillarySettings);
//
//      std::vector<std::vector<TimeType> > truncatedObservationTimes;
//      std::vector<std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > > truncatedObservables;
//
//      if (std::isnan(static_cast< double >( startAndEndTimesToProcess.first )) &&
//          std::isnan(static_cast< double >( startAndEndTimesToProcess.second ))) {
//        truncatedObservationTimes = observationTimes;
//        truncatedObservables = observables;
//      } else {
//        for (unsigned int i = 0; i < observationTimes.size(); ++i) {
//          std::vector<TimeType> singleTruncatedObservationTimes;
//          std::vector<Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1> > singleTruncatedObservables;
//
//          for (unsigned int j = 0; j < observationTimes.at(i).size(); ++j) {
//            if (((!std::isnan(static_cast< double >( startAndEndTimesToProcess.first ))
//                && observationTimes.at(i).at(j) >= startAndEndTimesToProcess.first) ||
//                std::isnan(static_cast< double >( startAndEndTimesToProcess.first ))) &&
//                ((!std::isnan(static_cast< double >( startAndEndTimesToProcess.second ))
//                    && observationTimes.at(i).at(j) <= startAndEndTimesToProcess.second) ||
//                    std::isnan(static_cast< double >( startAndEndTimesToProcess.second )))) {
//              singleTruncatedObservationTimes.push_back(observationTimes.at(i).at(j));
//              singleTruncatedObservables.push_back(observables.at(i).at(j));
//            }
//          }
//          truncatedObservationTimes.push_back(singleTruncatedObservationTimes);
//          truncatedObservables.push_back(singleTruncatedObservables);
//        }
//      }
//
//      // Create the single observation sets and save them
//      for (unsigned int i = 0; i < observationTimes.size(); ++i) {
//        sortedObservationSets[currentObservableType][currentLinkEnds].push_back(
//            std::make_shared<SingleObservationSet<ObservationScalarType, TimeType> >(
//                currentObservableType, currentLinkEnds, truncatedObservables.at(i), truncatedObservationTimes.at(i),
//                receiver,
//                std::vector<Eigen::VectorXd>(),
//                nullptr, std::make_shared<ObservationAncilliarySimulationSettings>(
//                    ancillarySettings.at(i))));
//      }
//    }
//  }
//
//  return std::make_shared<ObservationCollection<ObservationScalarType, TimeType> >(
//      sortedObservationSets);
//}

} // namespace observation_models
} // namespace tudat

#endif // TUDAT_PROCESSODFFILE_H
