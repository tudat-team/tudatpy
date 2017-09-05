/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TWOWAYDOPPLERPARTIAL_H
#define TUDAT_TWOWAYDOPPLERPARTIAL_H

#include <iostream>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayDopplerPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{


class TwoWayDopplerScaling: public PositionPartialScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param constituentRangeScalings Map of consitutent one-way Doppler scaling objects, with link end index as map key
     * \param numberOfLinkEnds Number of link ends in observable
     */
    TwoWayDopplerScaling( const std::map< int, boost::shared_ptr< OneWayDopplerScaling > > dopplerScalings )
    {
        if( dopplerScalings.count( 0 ) == 0 )
        {
            throw std::runtime_error( "Error when making two-way Doppler scaling, no uplink scaling found" );
        }
        else
        {
            uplinkDopplerScaling_ = dopplerScalings.at( 0 );
        }

        if( dopplerScalings.count( 1 ) == 0 )
        {
            throw std::runtime_error( "Error when making two-way Doppler scaling, no downlink scaling found" );
        }
        else
        {
            downlinkDopplerScaling_ = dopplerScalings.at( 1 );
        }
    }

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

    //! Function to get value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial.
    /*!
     * Function to get value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial,
     * for a single constituent one-way Doppler (as computed by last call to update function).
     * \param linkIndex Index of link for which scaling is to be returned.
     * \return Value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial.
     */
    double getProjectedRelativeVelocityRatio( const int linkIndex )
    {
        return projectedRelativeVelocityRatios_.at( linkIndex );
    }

private:

    boost::shared_ptr< OneWayDopplerScaling > uplinkDopplerScaling_;

    boost::shared_ptr< OneWayDopplerScaling > downlinkDopplerScaling_;

    //! List of values by which to scale constituent one-way ranges partials for it to be put into two-way range partial.
    std::map< int, double > projectedRelativeVelocityRatios_;
};

//! Class to compute the partial derivatives of an two-way range observation partial.
class TwoWayDopplerPartial: public ObservationPartial< 1 >
{

public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > TwoWayDopplerPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param twoWayDopplerScaler Scaling object used for mapping partials of one-way ranges to partials of observable
     * \param dopplerPartialList List of one-way range partials per link index.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param numberOfLinkEnds Number of link ends in two-way observable
     */
    TwoWayDopplerPartial( const boost::shared_ptr< TwoWayDopplerScaling > twoWayDopplerScaler,
                          const std::map< int, boost::shared_ptr< ObservationPartial< 1 > > >& dopplerPartialList,
                          const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                          const int numberOfLinkEnds ):
        ObservationPartial< 1 >( parameterIdentifier ), twoWayDopplerScaler_( twoWayDopplerScaler ), dopplerPartialList_( dopplerPartialList ),
        numberOfLinkEnds_( numberOfLinkEnds ){ }

    //! Destructor
    ~TwoWayDopplerPartial( ) { }

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
    TwoWayDopplerPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

protected:

    //! Scaling object used for mapping partials of one-way ranges to partials of observable
    boost::shared_ptr< TwoWayDopplerScaling > twoWayDopplerScaler_;

    //! List of one-way range partials per link index.
    std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > dopplerPartialList_;

    //! Predeclared iterator
    std::map< int, boost::shared_ptr< ObservationPartial< 1 > > >::iterator dopplerPartialIterator_;

    //! Number of link ends in two-way observable
    int numberOfLinkEnds_;
};

}

}


#endif // TUDAT_TWOWAYDOPPLERPARTIAL_H
