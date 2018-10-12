/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NWAYRANGEPARTIAL_H
#define TUDAT_NWAYRANGEPARTIAL_H

#include <functional>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{


//! Derived class for scaling three-dimensional position partial to n-way range observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to n-way range observable partial. Implementation is taken
 *  from Moyer(2000) for constituent one-way ranges, whicha re scaled (at O(v/c)) to n-way range implementation
 */
class NWayRangeScaling: public PositionPartialScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param constituentRangeScalings Map of consitutent one-way range scaling objects, with link end index as map key
     * \param numberOfLinkEnds Number of link ends in observable
     */
    NWayRangeScaling( const std::map< int, std::shared_ptr< OneWayRangeScaling > >& constituentRangeScalings,
                      const int numberOfLinkEnds ):
        constituentRangeScalings_( constituentRangeScalings ),
        numberOfLinkEnds_( numberOfLinkEnds ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    //! Function to get value by which to scale a constituent one-way ranges partial for it to be put into n-way range partial.
    /*!
     * Function to get value by which to scale a constituent one-way ranges partial for it to be put into n-way range partial,
     * for a single constituent one-way range (as compiuuted by last call to update function).
     * \param linkIndex Index of link for which scaling is to be returned.
     * \return Value by which to scale a constituent one-way ranges partial for it to be put into n-way range partial.
     */
    double getProjectedRelativeVelocityRatio( const int linkIndex )
    {
        return projectedRelativeVelocityRatios_.at( linkIndex );
    }

private:

    //! Map of consitutent one-way range scaling objects, with link index as map key
    std::map< int, std::shared_ptr< OneWayRangeScaling > > constituentRangeScalings_;

    //! Predeclared iterator (for efficiencty)
    std::map< int, std::shared_ptr< OneWayRangeScaling > >::iterator constituentRangeScalingIterator_;

    //! List of values by which to scale constituent one-way ranges partials for it to be put into n-way range partial.
    std::map< int, double > projectedRelativeVelocityRatios_;

    //! Number of link ends in observation
    int numberOfLinkEnds_;

};

//! Class to compute the partial derivatives of an n-way range observation partial.
class NWayRangePartial: public ObservationPartial< 1 >
{

public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > NWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param nWayRangeScaler Scaling object used for mapping partials of one-way ranges to partials of observable
     * \param rangePartialList List of one-way range partials per link index.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param numberOfLinkEnds Number of link ends in n-way observable
     */
    NWayRangePartial( const std::shared_ptr< NWayRangeScaling > nWayRangeScaler,
                      const std::map< int, std::shared_ptr< ObservationPartial< 1 > > >& rangePartialList,
                      const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                      const int numberOfLinkEnds ):
        ObservationPartial< 1 >( parameterIdentifier ), nWayRangeScaler_( nWayRangeScaler ), rangePartialList_( rangePartialList ),
        numberOfLinkEnds_( numberOfLinkEnds ){ }

    //! Destructor
    ~NWayRangePartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    NWayRangePartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

protected:

    //! Scaling object used for mapping partials of one-way ranges to partials of observable
    std::shared_ptr< NWayRangeScaling > nWayRangeScaler_;

    //! List of one-way range partials per link index.
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > > rangePartialList_;

    //! Predeclared iterator
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > >::iterator rangePartialIterator_;

    //! Number of link ends in n-way observable
    int numberOfLinkEnds_;
};

}

}


#endif // TUDAT_NWAYRANGEPARTIAL_H
