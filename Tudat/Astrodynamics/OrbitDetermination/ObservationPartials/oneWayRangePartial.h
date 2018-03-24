/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYRANGEPARTIAL_H
#define TUDAT_ONEWAYRANGEPARTIAL_H

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to one-way range observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to one-way range observable partial. Implementation is taken
 *  from Moyer(2000) and is separately implemented for fixed receiver and transmitter.
 */
class OneWayRangeScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~OneWayRangeScaling( ){ }

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
                 const Eigen::VectorXd currentObservation = Eigen::VectorXd::Constant( 1, TUDAT_NAN ) );

    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
    }

    //! Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
    /*!
     * Function to retrieve the factor by which the light-time partials should be scaled in one-way observation partial.
     * \return Factor by which the light-time partials should be scaled in one-way observation partial.
     */
    double getLightTimePartialScalingFactor( )
    {
       return referenceLightTimeCorrectionScaling_;
    }

    //! Function to get the fixed link end for last computation of update() function.
    /*!
     * Fixed link end for last computation of update() function.
     * \return Function to get the fixed link end for last computation of update() function.
     */
    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:

    //! Computed scaling factor (at receiver)
    Eigen::Matrix< double, 1, 3 > referenceScalingFactor_;

    //! Computed light time correction scaling factor
    double referenceLightTimeCorrectionScaling_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

};

//! Class to compute the partial derivatives of a one-way range observation partial.
class OneWayRangePartial: public ObservationPartial< 1 >
{

public:

    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayRangePartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param oneWayRangeScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    OneWayRangePartial(
            const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
            const std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > >& positionPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
            lighTimeCorrectionPartials =
            std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ), oneWayRangeScaler_( oneWayRangeScaler ),
        positionPartialList_( positionPartialList )
    {
        std::pair< boost::function< SingleOneWayRangePartialReturnType(
                    const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >,
                bool > lightTimeCorrectionPartial;

        // Create light time correction partial functions
        for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
        {
            lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
                        parameterIdentifier, lighTimeCorrectionPartials.at( i ) );
            if( lightTimeCorrectionPartial.second != 0 )
            {
                lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
            }
        }
    }
\
    //! Destructor.
    ~OneWayRangePartial( ) { }

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
    virtual OneWayRangePartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

    //! Function to get scaling object used for mapping partials of positions to partials of observable
    /*!
     * Function to get scaling object used for mapping partials of positions to partials of observable
     * \return
     */
    boost::shared_ptr< OneWayRangeScaling > getOneWayRangeScaler( )
    {
        return oneWayRangeScaler_;
    }

    //! Function to get the number of light-time correction partial functions.
    /*!
     * Number of light-time correction partial functions.
     * \return Number of light-time correction partial functions.
     */
    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
        return lighTimeCorrectionPartialsFunctions_.size( );
    }

protected:

    //! Scaling object used for mapping partials of positions to partials of observable
    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > > positionPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > >::iterator positionPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< boost::function< SingleOneWayRangePartialReturnType(
            const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctions_;

    //! List of light-time correction partial objects.
    std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartials_;

    //! Pre-declared state variable to be used in calculatePartial function.
    Eigen::Vector6d currentState_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;

};

}

}



#endif // TUDAT_ONEWAYRANGEPARTIAL_H
