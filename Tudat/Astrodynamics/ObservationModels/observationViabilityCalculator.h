#ifndef OBSERVATIONVIABILITYCALCULATOR_H
#define OBSERVATIONVIABILITYCALCULATOR_H

#include <vector>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"

#include "Tudat/Astrodynamics/GroundStations/pointingAnglesCalculator.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{

//! Base class for determining whether an observation is possible or not
/*!
 *  Base class for determining whether an observation is possible or not. Derived classes implement specific checks, such as minimum elevation
 *  angle, body occultation, etc. The input from which the viability of an observation is calculated is a vector od times and states of the link
 *  ends involved in the observation for which the viability is checked, in the order as provided by the computeObservationsAndLinkEndData of
 *  the associated ObservationModel.
 */
class ObservationViabilityCalculator
{
public:

    //! Base class constructor
    /*!
     *  Base class constructor (empty)
     */
    ObservationViabilityCalculator( ){ }

    //! Base class destructor
    /*!
     *  Base class destructor (empty).
     */
    virtual ~ObservationViabilityCalculator( ){ }

    //! Pure virtual base class function for determining whether an observation is viable.
    /*!
     *  Pure virtual base class function for determining whether an observation is viable. The input from which the viability of
     *  an observation is calculated is a vector of times and states of the link ends involved in the observation for which the viability is checked.
     *  Function to be implemented in derived class.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the  of
     *  computeObservationsAndLinkEndData the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the  of
     *  computeObservationsAndLinkEndData the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */
    virtual bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                      const std::vector< double >& linkEndTimes ) = 0;
};

bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times, const LinkEnds& linkEnds,
        const std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > >& viabilityCalculators );

bool isObservationViable(
        const std::vector< Eigen::Vector6d >& states, const std::vector< double >& times,
        const std::vector< boost::shared_ptr< ObservationViabilityCalculator > >& viabilityCalculators );


//! Function to check whether an observation is possible based on minimum elevation angle criterion.
/*!
 *  Function to check whether an observation is possible based on minimum elevation angle criterion at one link end.
 */
class MinimumElevationAngleCalculator: public ObservationViabilityCalculator
{
public:

    //! Constructor.
    /*!
     *  Constructor, takes a variable defining the geometry and the minimmum elevation angle and current angle calculator.
     *  \param linkEndIndices Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes vectors
     *  to use for elevation angle calculator when isObservationViable is called. The first entry of the pair is the index of the target that is being
     *  observed from the ground station at which the elevation angle is check, the second entry is the index of the ground station at which the check
     *  is performed. From each entry of this vector, a vector is created for which the elevation angle is checked.
     *  \param minimumElevationAngle Minimum elevation angle that is allowed.
     *  \param pointingAngleCalculator Object to calculate pointing angles (elevation angle) at ground station
     */
    MinimumElevationAngleCalculator( const std::vector< std::pair< int, int > > linkEndIndices,
                                     const double minimumElevationAngle,
                                     const boost::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator ):
        linkEndIndices_( linkEndIndices ), minimumElevationAngle_( minimumElevationAngle ),
        pointingAngleCalculator_( pointingAngleCalculator ){ }

    ~MinimumElevationAngleCalculator( ){ }

    //! Function for determining whether the elevation angle at station is sufficient to allow observation.
    /*!
     *  Function for determining whether the elevation angle at station is sufficient to allow observation. The input from which the viability of
     *  an observation is calculated is a vector of times and states of the link ends involved in the observation for which the viability is checked.
     *  \param linkEndStates Vector of states of the link ends involved in the observation, in the order as provided by the  of
     *  computeObservationsAndLinkEndData the associated ObservationModel.
     *  \param linkEndTimes Vector of times of the link ends involved in the observation, in the order as provided by the  of
     *  computeObservationsAndLinkEndData the associated ObservationModel.
     *  \return True if observation is viable, false if not.
     */

    bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                              const std::vector< double >& linkEndTimes );
private:

    //! Vector of indices denoting which combinations of entries of vectors to isObservationViable are to be used.
    /*!
     *  Vector of indices denoting which combinations of entries from the linkEndIndices and linkEndTimes vectors to use for
     *  elevation angle calculator when isObservationViable is called. The first entry of the pair is the index of the target that is being
     *  observed from the ground station at which the elevation angle is check, the second entry is the index of the ground station at which the
     *  check is performed. From each entry of this vector, a vector is created for which the elevation angle is checked.
     */
    std::vector< std::pair< int, int > > linkEndIndices_;

    //! Minimum elevation angle that is allowed.
    /*!
     *  Minimum elevation angle that is allowed.
     */
    double minimumElevationAngle_;

    //! Object to calculate pointing angles (elevation angle) at ground station
    /*!
     *  Object to calculate pointing angles (elevation angle) at ground station
     */
    boost::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator_;
};

class BodyAvoidanceAngleCalculator: public ObservationViabilityCalculator
{
public:

    BodyAvoidanceAngleCalculator( const std::vector< std::pair< int, int > > linkEndIndices,
                                  const double bodyAvoidanceAngle,
                                  const boost::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid,
                                  const std::string bodyToAvoid ):
        linkEndIndices_( linkEndIndices ),
        bodyAvoidanceAngle_( bodyAvoidanceAngle ),
        stateFunctionOfBodyToAvoid_( stateFunctionOfBodyToAvoid ),
    bodyToAvoid_( bodyToAvoid ){ }

    ~BodyAvoidanceAngleCalculator( ){ }

    bool isObservationViable( const std::vector< Eigen::Vector6d >& linkEndStates,
                                          const std::vector< double >& linkEndTimes );

private:

    std::vector< std::pair< int, int > > linkEndIndices_;

    double bodyAvoidanceAngle_;

    boost::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid_;

    std::string bodyToAvoid_;
};


//! Enum defining possible checks which can be performed for observation viability,
/*!
 *  Enum defining possible checks which can be performed for observation viability, the string and double parameter used to befine it when
 *  making an ObservationViabilitySettings object is given as comments.
 */
enum ObservationViabilityType
{
    minimum_elevation_angle, //properties: no string, double = elevation angle
    body_avoidance_angle //properties: string = body to avoid, double = avoidance angle
};

struct ObservationViabilitySettings
{
public:
    ObservationViabilitySettings( const ObservationViabilityType observationViabilityType,
                                  const std::pair< std::string, std::string > associatedLinkEnd,
                                  const std::string stringParameter,
                                  const double doubleParameter ):
        observationViabilityType_( observationViabilityType ), associatedLinkEnd_( associatedLinkEnd ),
        stringParameter_( stringParameter ), doubleParameter_( doubleParameter ){ }

    ObservationViabilityType observationViabilityType_;
    std::pair< std::string, std::string > associatedLinkEnd_;
    std::string stringParameter_;
    double doubleParameter_;
};

typedef std::map< ObservableType, std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > > >
PerObservableObservationViabilityCalculatorList;

typedef std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > >
PerLinkEndsObservationViabilityCalculatorList;

}

}

#endif // OBSERVATIONVIABILITYCALCULATOR_H
