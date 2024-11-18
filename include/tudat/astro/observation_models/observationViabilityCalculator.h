/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONVIABILITYCALCULATOR_H
#define TUDAT_OBSERVATIONVIABILITYCALCULATOR_H

#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/linearAlgebra.h"

#include "tudat/astro/basic_astro/missionGeometry.h"

#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{

namespace observation_models
{


//! Enum defining possible checks which can be performed for observation viability,
/*!
 *  Enum defining possible checks which can be performed for observation viability, the string and double parameter shown in the
 *  comments are used to define it when making an ObservationViabilitySettings object
 */
enum ObservationViabilityType
{
    minimum_elevation_angle,    //properties: no string, double = elevation angle
    body_avoidance_angle,       //properties: string = body to avoid, double = avoidance angle
    body_occultation            //properties: string = occulting body, no double
};


//! Base class for determining whether an observation is possible or not
/*!
 *  Base class for determining whether an observation is possible or not. Derived classes implement specific checks, such as
 *  minimum elevation angle, body occultation, etc. The input from which the viability of an observation is calculated is a vector
 *  of times and states of the link ends involved in the observation for which the viability is checked, in the order as provided
 *  by the computeObservationsAndLinkEndData of the associated ObservationModel.
 */
class ObservationViabilityCalculator
{
public:

    //! Base class constructor
    ObservationViabilityCalculator( ){ }

    //! Base class destructor
    virtual ~ObservationViabilityCalculator( ){ }

    //! Pure virtual base class function for determining whether an observation is viable.
    /*!
     *  Pure virtual base class function for determining whether an observation is viable. The input from which the viability of
     *  an observation is calculated is a vector of times and states of the link ends involved in the observation for which the
     *  viability is checked.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the  of
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */
    virtual bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                      const std::vector< double >& linkEndTimes ) = 0;
};

double getEvaluationEpochOfViabilityBody(
    const std::vector< Eigen::Vector6d >& linkEndStates,
    const std::vector< double >& linkEndTimes,
    const std::pair< int, int > linkEndIndexPair,
    const std::function< Eigen::Vector6d( const double ) > viabilityBodyStateFunction  );

//! Function to check whether an observation is viable
/*!
 * Function to check whether an observation is viable. The input from which the viability of an observation is calculated are a
 * vector  of times and states of the link ends involved in the observation.
 * \param times Vector of times of the link ends involved in the observation, in the order as provided by the
 * function computeObservationsAndLinkEndData of the associated ObservationModel.
 * \param states Vector of states of the link ends involved in the observation, in the order as provided by the
 * function computeObservationsAndLinkEndData of the associated ObservationModel.
 * \param linkEnds Link ends for current observation
 * \param viabilityCalculators List of viability calculators, for each set of link ends (function retrieves vector for linkEnds
 * input.
 * \return True if observation is viable, false if not.
 */
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times, const LinkEnds& linkEnds,
        const std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > >& viabilityCalculators );

//! Function to check whether an observation is viable
/*!
 * Function to check whether an observation is viable. The input from which the viability of an observation is calculated are a
 * vector  of times and states of the link ends involved in the observation.
 * \param times Vector of times of the link ends involved in the observation, in the order as provided by the
 * function computeObservationsAndLinkEndData of the associated ObservationModel.
 * \param states Vector of states of the link ends involved in the observation, in the order as provided by the
 * function computeObservationsAndLinkEndData of the associated ObservationModel.
 * \param viabilityCalculators List of viability calculators
 * \return True if observation is viable, false if not.
 */
bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times,
        const std::vector< std::shared_ptr< ObservationViabilityCalculator > >& viabilityCalculators );


//! Function to check whether an observation is possible based on minimum elevation angle criterion at one link end.
class MinimumElevationAngleCalculator: public ObservationViabilityCalculator
{
public:

    //! Constructor.
    /*!
     * Constructor, takes a variable defining the geometry and the minimmum elevation angle and current angle calculator.
     * \param linkEndIndices Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes
     * vectors to use for elevation angle calculator when isObservationViable is called. The second entry of the pair is the index
     * of the target that is being observed from the ground station at which the elevation angle is check, the first entry is
     * the index of the ground station at which the check  is performed. From each entry of this vector, a vector is created for
     * which the elevation angle is checked.
     * \param minimumElevationAngle Minimum elevation angle that is allowed.
     * \param pointingAngleCalculator Object to calculate pointing angles (elevation angle) at ground station
     */
    MinimumElevationAngleCalculator(
            const std::vector< std::pair< int, int > > linkEndIndices,
            const double minimumElevationAngle,
            const std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator ):
        linkEndIndices_( linkEndIndices ), minimumElevationAngle_( minimumElevationAngle ),
        pointingAngleCalculator_( pointingAngleCalculator ){ }

    //! Destructor
    ~MinimumElevationAngleCalculator( ){ }

    //! Function for determining whether the elevation angle at station is sufficient to allow observation.
    /*!
     *  Function for determining whether the elevation angle at station is sufficient to allow observation. The input from which
     *  the viability of an observation is calculated are a vector of times and states of link ends involved in the observation.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the  of
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the  of
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */
    bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                              const std::vector< double >& linkEndTimes );
private:

    //! Vector of indices denoting which combinations of entries of vectors are to be used in isObservationViable  function
    /*!
     *  Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes vectors to use for
     *  elevation angle calculator when isObservationViable is called. The second entry of the pair is the index of the target
     *  that is being observed from the ground station at which the elevation angle is check, the first entry is the index of the
     *  ground station at which the  check is performed. From each entry of this vector, a vector is created for which the
     *  elevation angle is checked.
     */
    std::vector< std::pair< int, int > > linkEndIndices_;

    //! Minimum elevation angle that is allowed.
    double minimumElevationAngle_;

    //! Object to calculate pointing angles (elevation angle) at ground station
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator_;
};

double computeMinimumLinkDistanceToPoint( const Eigen::Vector3d& observingBody,
                                          const Eigen::Vector3d& transmittingBody,
                                          const Eigen::Vector3d& relativePoint );

double computeCosineBodyAvoidanceAngle( const Eigen::Vector3d& observingBody,
                                        const Eigen::Vector3d& transmittingBody,
                                        const Eigen::Vector3d& bodyToAvoid );

double computeCosineBodyAvoidanceAngle( const std::vector< Eigen::Vector6d >& linkEndStates,
                                        const std::pair< int, int > observingAndTransmittingIndex,
                                        const Eigen::Vector3d& bodyToAvoid );

//! Function to check whether an observation is possible, based on body avoidance angle, as vied from single link end
/*!
 *  Function to check whether an observation is possible, based on body avoidance angle, as vied from single link end. A typical
 *  example of this is the sun avoidance (or solar separation) angle, on which a limit is typically imposed. Calculations in this
 *  class computes the angle between the line-of-sight vector(s) from a single link end, and the vecotr from this link end to the
 *  center of mass of the body that is to be 'avoided'.
 *  NOTE: This class computes the position of the body that is to be avoided in at the time halfway between transmission and
 *  reception time of the signal
 */
class BodyAvoidanceAngleCalculator: public ObservationViabilityCalculator
{
public:

    //! Constructor
    /*!
     *  \brief BodyAvoidanceAngleCalculator
     *  \param linkEndIndices Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes
     *  vectors to use for elevation angle calculator when isObservationViable is called. The second entry of the pair is the
     *  index of the target that is being observed from the ground station at which the avoidance angle is checked, the first
     *  entry is the index of the ground station at which the check is performed. From each entry of this vector, a vector is
     *  created for which the avoidance angle to bodyToAvoid_ is checked
     *  \param bodyAvoidanceAngle Minimum angle between line-of-sight vector and vector to avoided body that is allowed
     *  \param stateFunctionOfBodyToAvoid Function that returns the inertial state of the body that is to be avoided as a functiom
     *  of time
     *  \param bodyToAvoid Name of the body that is to be avoided.
     */
    BodyAvoidanceAngleCalculator( const std::vector< std::pair< int, int > > linkEndIndices,
                                  const double bodyAvoidanceAngle,
                                  const std::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid,
                                  const std::string bodyToAvoid ):
        linkEndIndices_( linkEndIndices ),
        bodyAvoidanceAngle_( bodyAvoidanceAngle ),
        stateFunctionOfBodyToAvoid_( stateFunctionOfBodyToAvoid ),
        bodyToAvoid_( bodyToAvoid ){ }

    //! Destructor
    ~BodyAvoidanceAngleCalculator( ){ }

    //! Function for determining whether the avoidance angle to a given body at station is sufficient to allow observation.
    /*!
     *  Function for determining whether the avoidance angle to a given body at station is sufficient to allow observation.
     *  The input from which the viability of an observation is calculated are a vector of times and states of link ends involved
     *  in the observation.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */
    bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                              const std::vector< double >& linkEndTimes );

private:

    //! Vector of indices denoting which combinations of entries of vectors to isObservationViable are to be used.
    /*!
     *  Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes vectors to use for
     *  elevation angle calculator when isObservationViable is called. The second entry of the pair is the index of the target that
     *  is being observed from the ground station at which the avoidance angle is checked, the first entry is the index of the
     *  ground station at which the check is performed. From each entry of this vector, a vector is created for which the
     *  avoidance angle to bodyToAvoid_ is checked.
     */
    std::vector< std::pair< int, int > > linkEndIndices_;

    //! Minimum angle between line-of-sight vector and vector to avoided body that is allowed
    double bodyAvoidanceAngle_;

    //! Function that returns the inertial state of the body that is to be avoided as a function of time
    std::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid_;

    //! Name of the body that is to be avoided.
    std::string bodyToAvoid_;
};

bool computeOccultation(
    const Eigen::Vector3d observer1Position,
    const Eigen::Vector3d observer2Position,
    const Eigen::Vector3d occulterPosition,
    const double radius );

//! Function to check whether an observation is possible, based on whether a given body is causing an occultation of the link
/*!
 *  Function to check whether an observation is possible, based on whether a given body is causing an occultation of the link.
 *  NOTE: This class computes the position of the occulting body in at the time halfway between transmission and reception time of
 *  the signal
 */
class OccultationCalculator: public ObservationViabilityCalculator
{
public:

    OccultationCalculator(
            const std::vector< std::pair< int, int > > linkEndIndices,
            const std::function< Eigen::Vector6d( const double ) > stateFunctionOfOccultingBody,
            const double radiusOfOccultingBody ):
        linkEndIndices_( linkEndIndices ),
        stateFunctionOfOccultingBody_( stateFunctionOfOccultingBody ),
        radiusOfOccultingBody_( radiusOfOccultingBody ){ }

    //! Function for determining whether the link is occulted during the observataion.
    /*!
     *  Function for determining whether the link is occulted during the observataion.
     *  The input from which the viability of an observation is calculated are a vector of times and states of link ends involved
     *  in the observation.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the
     *  function computeObservationsAndLinkEndData of the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */
    bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                              const std::vector< double >& linkEndTimes );

private:

    //! Vector of indices denoting which combinations of entries of vectors to isObservationViable are to be used.
    /*!
     *  Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes vectors to use for
     *  occultation calculation when isObservationViable is called. The second entry of the pair is the index of the target that
     *  is being observed from the ground station at which the signal is transmitted, the first entry is the index of the
     *  ground station at which it is received. From each entry of this vector, a vector is created for which the
     *  it if checked if it is occulted.
     */
    std::vector< std::pair< int, int > > linkEndIndices_;

    //! Function that returns the inertial state of the body that for which it is checked whether it occults the link.
    std::function< Eigen::Vector6d( const double ) > stateFunctionOfOccultingBody_;

    //! Radius of body causing occultation.
    double radiusOfOccultingBody_;
};

}

}

#endif // TUDAT_OBSERVATIONVIABILITYCALCULATOR_H
