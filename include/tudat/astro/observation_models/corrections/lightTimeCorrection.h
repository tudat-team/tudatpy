/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LIGHTTIMECORRECTION_H
#define TUDAT_LIGHTTIMECORRECTION_H

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observationModel.h"

namespace tudat
{

namespace observation_models
{

//! Enum defining different types of light time corrections.
enum LightTimeCorrectionType {
    first_order_relativistic,
    function_wrapper_light_time_correction,
    tabulated_tropospheric,
    saastamoinen_tropospheric,
    vmf_tropospheric,
    tabulated_ionospheric,
    jakowski_vtec_ionospheric,
    inverse_power_series_solar_corona
};

bool requiresMultiLegIterations( const LightTimeCorrectionType& lightTimeCorrectionType );

std::string getLightTimeCorrectionName( const LightTimeCorrectionType& lightTimeCorrectionType );

//! Base class for computing deviations from the light-time w.r.t. straight-line propagation at constant velocity (c).
/*!
 *  Base class for computing deviations from the light-time w.r.t. straight-line propagation at constant velocity (c).
 *  This base class is non-functional, and each time of light-time correction must be defined in a dedicated derived class
 */
class LightTimeCorrection
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrectionType Type of light-time correction represented by instance of class.
     */
    LightTimeCorrection( const LightTimeCorrectionType lightTimeCorrectionType ):
        lightTimeCorrectionType_( lightTimeCorrectionType ), timePerturbation_( 15.0 ), positionRelativePerturbation_( 1.0E-4 )
    { }

    //! Destructor
    virtual ~LightTimeCorrection( ) { }

    //! Pure virtual function to compute the light-time correction
    /*!
     * Pure virtual function to compute the light-time correction, function is to be implemented in derived class
     * for specific correction model. The input is the states and times of the two link ends (at the current iteration of
     * the solution of the implicit light-time equation), which is the information of the link already computed by the
     * light time calculator, ensuring that no double computations are performed.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \return Light-time correction
     */
    double calculateLightTimeCorrection(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr )
    {
        std::vector< Eigen::Vector6d > linkEndsStates = { transmitterState, receiverState };
        std::vector< double > linkEndsTimes = { transmissionTime, receptionTime };

        return calculateLightTimeCorrectionWithMultiLegLinkEndStates( linkEndsStates, linkEndsTimes, 0, ancillarySettings );
    }

    //! Pure virtual function to compute the light-time correction
    /*!
     * Pure virtual function to compute the light-time correction, function is to be implemented in derived class
     * for specific correction model.
     * @param linkEndsStates List of states at each link end during observation.
     * @param linkEndsTimes List of times at each link end during observation.
     * @param currentMultiLegTransmitterIndex Index in the linkEndsStates and linkEndsTimes of the transmitter in the current link.
     * @param ancillarySettings Observation ancillary simulation settings.
     * @return
     */
    virtual double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex = 0,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr ) = 0;

    //! Pure virtual function to compute the partial derivative of the light-time correction w.r.t. observation time
    /*!
     * Pure virtual function to compute the partial derivative of the light-time correction w.r.t. observation time.
     * Function is to be implemented in derived class for specific correction model.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Light-time correction w.r.t. observation time
     */
    virtual double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings );

    //! Pure virtual function to compute the partial derivative of the light-time correction w.r.t. link end position
    /*!
     * Pure virtual function to compute the partial derivative of the light-time correction w.r.t. link end position
     * Function is to be implemented in derived class for specific correction model.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Light-time correction w.r.t. link end position
     */
    virtual Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings );

    //! Function to retrieve the type of light-time correction represented by instance of class.
    /*!
     *  Function to retrieve the type of light-time correction represented by instance of class.
     *  \return Type of light-time correction represented by instance of class.
     */
    LightTimeCorrectionType getLightTimeCorrectionType( )
    {
        return lightTimeCorrectionType_;
    }

protected:
    void getTransmissionReceptionTimesAndStates( const std::vector< Eigen::Vector6d >& linkEndsStatesInput,
                                                 const std::vector< double >& linkEndsTimesInput,
                                                 const unsigned int currentMultiLegTransmitterIndex,
                                                 Eigen::Vector6d& transmitterStateOutput,
                                                 Eigen::Vector6d& receiverStateOutput,
                                                 double& transmissionTimeOutput,
                                                 double& receptionTimeOutput )
    {
        const unsigned int currentMultiLegReceiverIndex = currentMultiLegTransmitterIndex + 1;

        if( currentMultiLegReceiverIndex >= linkEndsStatesInput.size( ) || currentMultiLegReceiverIndex >= linkEndsTimesInput.size( ) )
        {
            throw std::runtime_error(
                    "Error when getting receiver and transmitter states and times in LightTimeCorrection: "
                    "specified transmitter index is invalid" );
        }

        transmitterStateOutput = linkEndsStatesInput.at( currentMultiLegTransmitterIndex );
        receiverStateOutput = linkEndsStatesInput.at( currentMultiLegReceiverIndex );
        transmissionTimeOutput = linkEndsTimesInput.at( currentMultiLegTransmitterIndex );
        receptionTimeOutput = linkEndsTimesInput.at( currentMultiLegReceiverIndex );
    }

    //! Type of light-time correction represented by instance of class.
    LightTimeCorrectionType lightTimeCorrectionType_;

    double timePerturbation_;

    double positionRelativePerturbation_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_LIGHTTIMECORRECTION_H
