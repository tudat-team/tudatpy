/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONVIABILITY_H
#define TUDAT_CREATEOBSERVATIONVIABILITY_H

#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"

namespace tudat
{

namespace observation_models
{

//! Class to define settings for observation viability calculator creation
/*!
 *  Class to define settings for observation viability calculator creation. Settinsg are defined by type of viability check,
 *  link end at which the check is to be performed, as well as double and string specifiers, which provide any needed
 *  specific parameters of the vaibility check (meaning of both parameters is given in constructor doxygen, as well as
 *  in-code comments of ObservationViabilityType
 */
class ObservationViabilitySettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param observationViabilityType Type of viability that is to be checked
     * \param associatedLinkEnd Link end at which viability is to be checked
     * \param stringParameter String parameter providing additional specification for viability check, meaning is type-dependent:
     * Elevation angle: none; Avoidance angle: body for which avoidance angle is to be calculated; Occultation: body for which
     * occultation is to be calculated
     * \param doubleParameter Double parameter providing additional specification for viability check, meaning is type-dependent:
     * Elevation angle: minimum allowed elevation angle; Avoidance angle: minimum allowed avoidance angle; Occultation: none
     */
    ObservationViabilitySettings( const ObservationViabilityType observationViabilityType,
                                  const std::pair< std::string, std::string > associatedLinkEnd,
                                  const std::string stringParameter = "",
                                  const double doubleParameter = TUDAT_NAN ):
        observationViabilityType_( observationViabilityType ), associatedLinkEnd_( associatedLinkEnd ), stringParameter_( stringParameter ),
        doubleParameter_( doubleParameter )
    {}

    virtual ~ObservationViabilitySettings( ) = default;

    //! Type of viability that is to be checked
    ObservationViabilityType observationViabilityType_;

    std::pair< std::string, std::string > getAssociatedLinkEnd( )
    {
        return associatedLinkEnd_;
    }

    std::string getStringParameter( )
    {
        return stringParameter_;
    }

    double getDoubleParameter( )
    {
        return doubleParameter_;
    }

protected:
    //! Link end at which viability is to be checked
    std::pair< std::string, std::string > associatedLinkEnd_;

    //! String parameter providing additional specification for viability check.
    /*!
     * String parameter providing additional specification for viability check, meaning is type-dependent:
     * Elevation angle: none; Avoidance angle: body for which avoidance angle is to be calculated; Occultation: body for which
     * occultation is to be calculated
     */
    std::string stringParameter_;

    //! Double parameter providing additional specification for viability check
    /*!
     *  Double parameter providing additional specification for viability check, meaning is type-dependent:
     *  Elevation angle: minimum allowed elevation angle; Avoidance angle: minimum allowed avoidance angle; Occultation: none
     */
    double doubleParameter_;
};

class ObservationBoundariesViabilitySettings : public ObservationViabilitySettings
{
public:
    ObservationBoundariesViabilitySettings( const std::pair< std::string, std::string > associatedLinkEnd,
                                            const std::vector< std::pair< double, double > > boundaries ):
        ObservationViabilitySettings( observation_boundaries, associatedLinkEnd ), boundaries_( boundaries )
    {
        for( unsigned int i = 0; i < boundaries.size( ); i++ )
        {
            if( boundaries.at( i ).first > boundaries.at( i ).second )
            {
                throw std::runtime_error( "Error when making observation boundaries viability settings, boundary pair is inconsistent: " +
                                          std::to_string( boundaries.at( i ).first ) + " is larger than " +
                                          std::to_string( boundaries.at( i ).second ) );
            }
        }
    }

    const std::vector< std::pair< double, double > >& getBoundaries( ) const
    {
        return boundaries_;
    }

private:
    std::vector< std::pair< double, double > > boundaries_;
};

class BodyInSunlightViabilitySettings : public ObservationViabilitySettings
{
public:
    BodyInSunlightViabilitySettings( const std::pair< std::string, std::string > associatedLinkEnd,
                                     const std::vector< std::string > occultingBodies ):
        ObservationViabilitySettings( body_in_sunlight, associatedLinkEnd ), occultingBodies_( occultingBodies )
    {}

    const std::vector< std::string >& getOccultingBodies( ) const
    {
        return occultingBodies_;
    }

private:
    std::vector< std::string > occultingBodies_;
};

inline std::shared_ptr< ObservationViabilitySettings > observationBoundariesViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const std::vector< std::pair< double, double > > boundaries )
{
    return std::make_shared< ObservationBoundariesViabilitySettings >( associatedLinkEnd, boundaries );
}

inline std::vector< std::shared_ptr< ObservationViabilitySettings > > observationBoundariesViabilitySettings(
        const std::vector< std::pair< std::string, std::string > > associatedLinkEnds,
        const std::vector< std::pair< double, double > > boundaries )
{
    std::vector< std::shared_ptr< ObservationViabilitySettings > > viabilitySettingsList;
    for( unsigned int i = 0; i < associatedLinkEnds.size( ); i++ )
    {
        viabilitySettingsList.push_back(
                std::make_shared< ObservationBoundariesViabilitySettings >( associatedLinkEnds.at( i ), boundaries ) );
    }
    return viabilitySettingsList;
}

inline std::vector< std::shared_ptr< ObservationViabilitySettings > > elevationAngleViabilitySettings(
        const std::vector< std::pair< std::string, std::string > > associatedLinkEnds,
        const double elevationAngle )
{
    std::vector< std::shared_ptr< ObservationViabilitySettings > > viabilitySettingsList;
    for( unsigned int i = 0; i < associatedLinkEnds.size( ); i++ )
    {
        viabilitySettingsList.push_back( std::make_shared< ObservationViabilitySettings >(
                minimum_elevation_angle, associatedLinkEnds.at( i ), "", elevationAngle ) );
    }
    return viabilitySettingsList;
}

inline std::shared_ptr< ObservationViabilitySettings > elevationAngleViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const double elevationAngle )
{
    return std::make_shared< ObservationViabilitySettings >( minimum_elevation_angle, associatedLinkEnd, "", elevationAngle );
}

inline std::shared_ptr< ObservationViabilitySettings > groundStationDarknessViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const double maximumSunElevationAngle = -12.0 * mathematical_constants::PI / 180.0 )
{
    return std::make_shared< ObservationViabilitySettings >( ground_station_darkness, associatedLinkEnd, "", maximumSunElevationAngle );
}

inline std::vector< std::shared_ptr< ObservationViabilitySettings > > bodyAvoidanceAngleViabilitySettings(
        const std::vector< std::pair< std::string, std::string > > associatedLinkEnds,
        const std::string bodyToAvoid,
        const double avoidanceAngle )
{
    std::vector< std::shared_ptr< ObservationViabilitySettings > > viabilitySettingsList;
    for( unsigned int i = 0; i < associatedLinkEnds.size( ); i++ )
    {
        viabilitySettingsList.push_back( std::make_shared< ObservationViabilitySettings >(
                body_avoidance_angle, associatedLinkEnds.at( i ), bodyToAvoid, avoidanceAngle ) );
    }
    return viabilitySettingsList;
}

inline std::shared_ptr< ObservationViabilitySettings > bodyAvoidanceAngleViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const std::string bodyToAvoid,
        const double avoidanceAngle )
{
    return std::make_shared< ObservationViabilitySettings >( body_avoidance_angle, associatedLinkEnd, bodyToAvoid, avoidanceAngle );
}

inline std::vector< std::shared_ptr< ObservationViabilitySettings > > bodyOccultationViabilitySettings(
        const std::vector< std::pair< std::string, std::string > > associatedLinkEnds,
        const std::string occultingBody )
{
    std::vector< std::shared_ptr< ObservationViabilitySettings > > viabilitySettingsList;
    for( unsigned int i = 0; i < associatedLinkEnds.size( ); i++ )
    {
        viabilitySettingsList.push_back(
                std::make_shared< ObservationViabilitySettings >( body_occultation, associatedLinkEnds.at( i ), occultingBody ) );
    }
    return viabilitySettingsList;
}

inline std::shared_ptr< ObservationViabilitySettings > bodyOccultationViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const std::string occultingBody )
{
    return std::make_shared< ObservationViabilitySettings >( body_occultation, associatedLinkEnd, occultingBody );
}

inline std::shared_ptr< ObservationViabilitySettings > bodyInSunlightViabilitySettings(
        const std::pair< std::string, std::string > associatedLinkEnd,
        const std::vector< std::string > occultingBodies )
{
    return std::make_shared< BodyInSunlightViabilitySettings >( associatedLinkEnd, occultingBodies );
}

//! Typedef for vector of ObservationViabilitySettings pointers
typedef std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > ObservationViabilitySettingsList;

//! Function to filter list of observationViabilitySettings, so that only those relevant for single set of link ends are retained
/*!
 * Function to filter list of observationViabilitySettings, so that only those relevant for single set of link ends are retained
 * \param observationViabilitySettings Full, unfiltered, list of observation viability settings
 * \param linkEnds Link ends for which the relevant observation vaibilityies are to be retrieved
 * \return List of observationViabilitySettings that are relevant for linkEnds
 */
ObservationViabilitySettingsList filterObservationViabilitySettings( const ObservationViabilitySettingsList& observationViabilitySettings,
                                                                     const LinkEnds& linkEnds );

//! Function to create an object to check if observation value is within given boundaries
/*!
 * Function to create an object to check if observation value is within given boundaries
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds Link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings Object that defines the settings for the creation of the viability check creation
 * (settings must be compatible with observation boundaries check). If ground station is not specified (by associatedLinkEnd_.second in
 * observationViabilitySettings), check is performed for all ground stations on (or c.o.m. of) body (defined by associatedLinkEnd_.first)
 * automatically. \return Object to check if observation value is within given boundaries
 */
std::shared_ptr< ObservationBoundariesViabilityCalculator > createObservationBoundariesCalculator(
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings );

//! Function to create an object to check if a minimum elevation angle condition is met for an observation
/*!
 * Function to create an object to check if a minimum elevation angle condition is met for an observation
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds Link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings Object that defines the settings for the creation of the viability check creation
 * (settings must be compatible with minimum elevation angle check).  Ground station must ve specified by
 * associatedLinkEnd_.second in observationViabilitySettings.
 * \param stationName Name of the ground station for which calculator is to be computed (if no station is explicitly given in
 * observationViabilitySettings).
 * \return Object to check if a minimum elevation angle condition is met for an observation
 */
std::shared_ptr< MinimumElevationAngleCalculator > createMinimumElevationAngleCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings,
        const std::string& stationName );

std::shared_ptr< GroundStationDarknessCalculator > createGroundStationDarknessCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings,
        const std::string& stationName );

//! Function to create an object to check if a body avoidance angle condition is met for an observation
/*!
 * Function to create an object to check if a body avoidance angle condition is met for an observation
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds Link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings Object that defines the settings for the creation of the viability check creation
 * (settings must be compatible with body avoidance angle check). If ground station is not specified (by
 * associatedLinkEnd_.second in observationViabilitySettings), check is performed for all ground stations on (or c.o.m. of) body
 * (defined by associatedLinkEnd_.first) automatically.
 * \return Object to check if a body avoidance angle condition is met for an observation
 */
std::shared_ptr< BodyAvoidanceAngleCalculator > createBodyAvoidanceAngleCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings );

//! Function to create an object to check if a body occultation condition is met for an observation
/*!
 * Function to create an object to check if a body occultation condition is met for an observation
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds Link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings Object that defines the settings for the creation of the viability check creation
 * (settings must be compatible with body occultation check).  If ground station is not specified (by
 * associatedLinkEnd_.second in observationViabilitySettings), check is performed for all ground stations on (or c.o.m. of) body
 * (defined by associatedLinkEnd_.first) automatically, or fo
 * \return Object to check if a body occultation condition is met for an observation
 */
std::shared_ptr< OccultationCalculator > createOccultationCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings );

std::shared_ptr< BodyInSunlightCalculator > createBodyInSunlightCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings );

//! Function to create an list of obervation viability conditions for a single set of link ends
/*!
 * Function to create an list of obervation viability conditions for a single set of link ends
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds Link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings List of viability settings from which viability check objects are to be created
 * \return List of obervation viability conditions for a single set of link ends
 */
std::vector< std::shared_ptr< ObservationViabilityCalculator > > createObservationViabilityCalculators(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings );

//! Function to create an list of obervation viability conditions for a number of sets of link ends, for a single observable type
/*!
 * Function to create an list of obervation viability conditions for a number of sets of link ends, for a single observable type
 * \param bodies Map of body objects that constitutes the environment
 * \param linkEnds List of link ends for which viability check object is to be made
 * \param observationType Type of observable for which viability check object is to be made
 * \param observationViabilitySettings List of viability settings from which viability check objects are to be created
 * \return List of obervation viability conditions for a number of sets of link ends, for a single observable type
 */
std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > > createObservationViabilityCalculators(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< LinkEnds > linkEnds,
        const ObservableType observationType,
        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings );

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CREATEOBSERVATIONVIABILITY_H
