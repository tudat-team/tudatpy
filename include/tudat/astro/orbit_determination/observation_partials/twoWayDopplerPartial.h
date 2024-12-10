/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

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
     * \param dopplerScalings List of consitutent one-way Doppler scaling objects, with entry i indicating link i
     * (0=transmitter to reflector, etc.)
     * \param rangeScalings List of consitutent one-way range scaling objects, with entry i indicating link i
     * (0=transmitter to reflector, etc.)
     * \param oneWayDopplerModels List of double precision observation functions for one-way Doppler, with entry i indicating
     * link i (0=transmitter to reflector, etc.)
     */
    TwoWayDopplerScaling( const std::vector< std::shared_ptr< OneWayDopplerScaling > > dopplerScalings,
                          const std::vector< std::shared_ptr< OneWayRangeScaling > > rangeScalings,
                          const std::vector< std::function< double( const double, const observation_models::LinkEndType ) > >
                          oneWayDopplerModels,
                          const double multiplicationTerm ): multiplicationTerm_( multiplicationTerm )
    {
        if( dopplerScalings.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making two-way Doppler scaling, one-way Doppler scaling size is inconsistent" );
        }
        uplinkDopplerScaling_ = dopplerScalings.at( 0 );
        downlinkDopplerScaling_ = dopplerScalings.at( 1 );

        if( rangeScalings.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making two-way Doppler scaling, one-way range scaling size is inconsistent" );
        }
        uplinkRangeScaling_ = rangeScalings.at( 0 );
        downlinkRangeScaling_ = rangeScalings.at( 1 );

        if( oneWayDopplerModels.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making two-way Doppler scaling, one-way modelg size is inconsistent" );
        }
        uplinkDopplerModel_ = oneWayDopplerModels.at( 0 );
        downlinkDopplerModel_ = oneWayDopplerModels.at( 1 );
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

    //! Funtion to get the relevant derivative of one-way Doppler observation, when setting up tow-way Doppler mode
    /*!
     * Funtion to get the relevant derivative of one-way Doppler observation, when setting up tow-way Doppler mode. For a fixed
     * transmitter, this is the downlink, for a fixed receiver it is the uplink. For a fixed reflector, no derivatives are
     * relevant and zero is returned
     * \param fixedLinkEnd Link end that is fixed when computing the two-way observable
     * \return Relevant one-way link Doppler obervation partial (multiplied by -1 if relevant partial derivative).
     */
    double getRelevantOneWayDopplerTimePartial( const observation_models::LinkEndType fixedLinkEnd )
    {
        if( fixedLinkEnd == observation_models::transmitter )
        {
            return downlinkOneWayDopplerTimeDerivative_;
        }
        else if( fixedLinkEnd == observation_models::receiver )
        {
            return -uplinkOneWayDopplerTimeDerivative_;
        }
        else
        {
            return 0.0;
        }
    }

    //! Function to get value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial.
    /*!
     * Function to get value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial,
     * for a single constituent one-way Doppler (as computed by last call to update function).
     * \param linkIndex Index of link for which scaling is to be returned.
     * \return Value by which to scale a constituent one-way Doppler partial for it to be put into two-way Doppler partial.
     */
    double getProjectedRelativeVelocityRatio( const int linkIndex )
    {
        return projectedRelativeVelocityRatios_.at( linkIndex ) * multiplicationTerm_;
    }

private:

    double multiplicationTerm_;

    //! One-way Doppler partial scaling object for uplink
    std::shared_ptr< OneWayDopplerScaling > uplinkDopplerScaling_;

    //! One-way Doppler partial scaling object for downlink
    std::shared_ptr< OneWayDopplerScaling > downlinkDopplerScaling_;


    //! One-way range partial scaling object for uplink
    std::shared_ptr< OneWayRangeScaling > uplinkRangeScaling_;

    //! One-way range partial scaling object for downlink
    std::shared_ptr< OneWayRangeScaling > downlinkRangeScaling_;


    //! Function to return the uplink one-way Doppler observation
    std::function< double( const double, const observation_models::LinkEndType ) > uplinkDopplerModel_;

    //! Function to return the downlink one-way Doppler observation
    std::function< double( const double, const observation_models::LinkEndType ) > downlinkDopplerModel_;


    //! Current partial time derivative of one-way Doppler for uplink.
    double uplinkOneWayDopplerTimeDerivative_;

    //! Current partial time derivative of one-way Doppler for downlink.
    double downlinkOneWayDopplerTimeDerivative_;


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
     * \param rangePartialList List of one-way range per link index.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param numberOfLinkEnds Number of link ends in two-way observable
     */
    TwoWayDopplerPartial( const std::shared_ptr< TwoWayDopplerScaling > twoWayDopplerScaler,
                          const std::map< int, std::shared_ptr< ObservationPartial< 1 > > >& dopplerPartialList,
                          const std::map< int, std::shared_ptr< ObservationPartial< 1 > > >& rangePartialList,
                          const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                          const int numberOfLinkEnds,
                          const std::function< double(
                              const observation_models::LinkEndType, const std::vector< Eigen::Vector6d >&,
                              const std::vector< double >&, const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) > scalingFactorFunction = nullptr ):
        ObservationPartial< 1 >( parameterIdentifier ), twoWayDopplerScaler_( twoWayDopplerScaler ), dopplerPartialList_( dopplerPartialList ),
        rangePartialList_( rangePartialList ), numberOfLinkEnds_( numberOfLinkEnds ), scalingFactorFunction_( scalingFactorFunction ){ }

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
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

protected:

    //! Scaling object used for mapping partials of one-way ranges to partials of observable
    std::shared_ptr< TwoWayDopplerScaling > twoWayDopplerScaler_;

    //! List of one-way range partials per link index.
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > > dopplerPartialList_;

    //! Predeclared iterator
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > >::iterator dopplerPartialIterator_;


    //! List of one-way range partials per link index.
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > > rangePartialList_;

    //! Predeclared iterator
    std::map< int, std::shared_ptr< ObservationPartial< 1 > > >::iterator rangePartialIterator_;


    //! Number of link ends in two-way observable
    int numberOfLinkEnds_;

    std::function< double(
        const observation_models::LinkEndType, const std::vector< Eigen::Vector6d >&,
        const std::vector< double >&, const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ) > scalingFactorFunction_;

};

}

}


#endif // TUDAT_TWOWAYDOPPLERPARTIAL_H
