#ifndef DIFFERENCEDONEWAYRANGERATEPARTIAL_H
#define DIFFERENCEDONEWAYRANGERATEPARTIAL_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to one-way range-rate (differenced) observable partial
class OneWayRangeRateScaling: public PositionPartialScaling
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param oneWayRangeScalerArcStart Partial scaling for arc start range observation
     * \param oneWayRangeScalerArcEnd Partial scaling for arc end range observation
     */
    OneWayRangeRateScaling(
            const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcStart,
            const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcEnd ):
        oneWayRangeScalerArcStart_( oneWayRangeScalerArcStart ), oneWayRangeScalerArcEnd_( oneWayRangeScalerArcEnd ){ }

    //! Destructor
    ~OneWayRangeRateScaling( ){ }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation Index of vector maps to link end for a
     *  given ObsevableType through getLinkEndIndex function.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd )
    {
        oneWayRangeScalerArcStart_->update(
                    std::vector< Eigen::Vector6d >(
                        linkEndStates.begin( ), linkEndStates.begin( ) + 2 ), times, fixedLinkEnd );
        oneWayRangeScalerArcEnd_->update(
                    std::vector< Eigen::Vector6d >(
                        linkEndStates.begin( ) + 2, linkEndStates.begin( ) + 4 ), times, fixedLinkEnd );
    }

    //! Partial scaling for arc start range observation
    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcStart_;

    //! Partial scaling for arc end range observation
    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcEnd_;
};

//! Class to compute the partial derivatives of a one-way range-rate (differenced) observation
class DifferencedOneWayRangeRatePartial: public ObservationPartial< 1 >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param arcStartRangePartial Partial object for arc start range observation
     * \param arcEndRangePartial Partial object for arc end range observation
     */
    DifferencedOneWayRangeRatePartial(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const boost::shared_ptr< ObservationPartial< 1 > > arcStartRangePartial,
            const boost::shared_ptr< ObservationPartial< 1 > > arcEndRangePartial):
        ObservationPartial< 1 >( parameterIdentifier ),
        arcStartRangePartial_( arcStartRangePartial ),
        arcEndRangePartial_( arcEndRangePartial ){ }

    //! Destructor
    ~DifferencedOneWayRangeRatePartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime );

protected:

    //! Partial object for arc start range observation
    boost::shared_ptr< ObservationPartial< 1 > > arcStartRangePartial_;

    //! Partial object for arc end range observation
    boost::shared_ptr< ObservationPartial< 1 > > arcEndRangePartial_;
};

}

}


#endif // DIFFERENCEDONEWAYRANGERATEPARTIAL_H
