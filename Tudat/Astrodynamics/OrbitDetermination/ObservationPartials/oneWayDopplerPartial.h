#ifndef ONEWAYDOPPLERPARTIAL_H
#define ONEWAYDOPPLERPARTIAL_H


#include <iostream>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayDopplerObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Derived class for scaling three-dimensional position partial to one-way doppler observable partial
/*!
 *  Derived class for scaling three-dimensional position partial to one-way doppler observable partial. Implementation is taken
 *  from Moyer(2000) and is separately implemented for fixed receiver and transmitter.
 */
class OneWayDopplerScaling: public PositionPartialScaling
{
public:

    //! Destructor
    ~OneWayDopplerScaling( ){ }

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
                 const observation_models::LinkEndType fixedLinkEnd );


    //! Function to retrieve the scaling factor for specific link end
    /*!
     * Function to retrieve the scaling factor for specific link end
     * \param linkEndType Link end for which scaling factor is to be returned
     * \return Position partial scaling factor at current link end
     */
    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return positionScalingFactor_ * ( ( linkEndType == observation_models::receiver ) ? ( -1.0 ) : ( 1.0 ) );
    }

    Eigen::Matrix< double, 1, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return ( ( linkEndType == observation_models::receiver ) ? ( transmitterVelocityScalingFactor_ ) : ( receiverVelocityScalingFactor_ ) );
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
    Eigen::Matrix< double, 1, 3 > positionScalingFactor_;

    Eigen::Matrix< double, 1, 3 > receiverVelocityScalingFactor_;

    Eigen::Matrix< double, 1, 3 > transmitterVelocityScalingFactor_;

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;

};

//! Class to compute the partial derivatives of a one-way doppler observation partial.
class OneWayDopplerPartial: public ObservationPartial< 1 >
{

public:

    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayDopplerPartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayDopplerPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param oneWayDopplerScaler Scaling object used for mapping partials of positions to partials of observable
     * \param positionPartialList List of position partials per link end.
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     * \param lighTimeCorrectionPartials List if light-time correction partial objects.
     */
    OneWayDopplerPartial(
            const boost::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler,
            const std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >& positionPartialList,
            const std::map< observation_models::LinkEndType, boost::shared_ptr< VelocityPartial > >& velocityPartialList,
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
            lighTimeCorrectionPartials =
            std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ), oneWayDopplerScaler_( oneWayDopplerScaler ),
        positionPartialList_( positionPartialList ), velocityPartialList_( velocityPartialList )
    {
        std::pair< boost::function< SingleOneWayDopplerPartialReturnType(
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
    ~OneWayDopplerPartial( ) { }

    //! Function to calculate the observation partial(s) at required time and state
    /*!
     *  Function to calculate the observation partial(s) at required time and state. State and time
     *  are typically obtained from evaluation of observation model.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end time.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual OneWayDopplerPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime );

    //! Function to get scaling object used for mapping partials of positions to partials of observable
    /*!
     * Function to get scaling object used for mapping partials of positions to partials of observable
     * \return
     */
    boost::shared_ptr< OneWayDopplerScaling > getOneWayDopplerScaler( )
    {
        return oneWayDopplerScaler_;
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
    boost::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler_;

    //! List of position partials per link end.
    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartialList_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< VelocityPartial > > velocityPartialList_;

    //! Iterator over list of position partials per link end.
    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >::iterator positionPartialIterator_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< VelocityPartial > >::iterator velocityPartialIterator_;

    //! List of light-time correction partial functions.
    std::vector< boost::function< SingleOneWayDopplerPartialReturnType(
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

#endif // ONEWAYDOPPLERPARTIAL_H
