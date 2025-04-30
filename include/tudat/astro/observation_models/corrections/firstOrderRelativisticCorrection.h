/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H
#define TUDAT_FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H

#include <cmath>
#include <vector>

#include <functional>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Class to calculate first order relativistic light time correction (Shapiro time delay) due to a set of point masses.
/*!
 *  Class to calculate first order relativistic light time correction (Shapiro time delay) due to a set of point masses.
 *  This class has the properties (mass and position) of the gravitating bodies as member variables. It receives the states
 *  and times of  receiver and transmitter as input each time a light time is to be calculated.
 *  NOTE: If the perturbing body is equal to the transmitting or receiving body, the transmission or reception time,
 *  respectively, are used to evaluate this perturbing body's state when computing the correction. For all other bodies,
 *  the mid-way point of the light time is used.
 */
class FirstOrderLightTimeCorrectionCalculator : public LightTimeCorrection
{
public:
    //! Constructor, takes and sets gravitating body properties.
    /*!
     *  Constructor, takes and sets gravitating body properties.
     *  \param perturbingBodyStateFunctions Set of function returning the state of the gravitating bodies as a function
     *  of time.
     *  \param perturbingBodyGravitationalParameterFunctions Set of functions returning the gravitational parameters of
     *  the gravitating bodies.
     *  \param perturbingBodyNames Names of bodies causing light-time correction.
     *  \param transmittingBody Name of transmitting body
     *  \param receivingBody Name of receiving body
     *  \param ppnParameterGammaFunction Function returning the parametric post-Newtonian parameter gamma, a measure
     *  for the space-time curvature due to a unit rest mass (default 1.0; value from GR)
     *  \param bending Boolean flag to determine if light bending should be included.
     */
    FirstOrderLightTimeCorrectionCalculator(
            const std::vector< std::function< Eigen::Vector6d( const double ) > >& perturbingBodyStateFunctions,
            const std::vector< std::function< double( ) > >& perturbingBodyGravitationalParameterFunctions,
            const std::vector< std::string > perturbingBodyNames,
            const std::string transmittingBody,
            const std::string receivingBody,
            const std::function< double( ) >& ppnParameterGammaFunction = []( ) { return 1.0; },
            const bool bending = false ):
        LightTimeCorrection( first_order_relativistic ), perturbingBodyStateFunctions_( perturbingBodyStateFunctions ),
        perturbingBodyGravitationalParameterFunctions_( perturbingBodyGravitationalParameterFunctions ),
        perturbingBodyNames_( perturbingBodyNames ), ppnParameterGammaFunction_( ppnParameterGammaFunction ), bending_( bending )
    {
        currentTotalLightTimeCorrection_ = 0.0;
        currentLighTimeCorrectionComponents_.resize( perturbingBodyNames_.size( ) );

        // Check if perturbing body is transmitting/receiving body, and set evaluation time settings accordingly
        for( unsigned int i = 0; i < perturbingBodyNames.size( ); i++ )
        {
            if( perturbingBodyNames.at( i ) == transmittingBody )
            {
                lightTimeEvaluationContribution_.push_back( 0.0 );
            }
            else if( perturbingBodyNames.at( i ) == receivingBody )
            {
                lightTimeEvaluationContribution_.push_back( 1.0 );
            }
            else
            {
                lightTimeEvaluationContribution_.push_back( 0.5 );
            }
        }
    }

    //! Destructor
    ~FirstOrderLightTimeCorrectionCalculator( ) { }

    //! Function to calculate first-order relativistic light time correction due to set of gravitating point masses.
    /*!
     *  Function to calculate first order relativistic light time correction due to set of gravitating point masses,
     *  according to Eq. (11.17) of 2010 IERS conventions. Calculation are performed by calling ca
     *  calculateFirstOrderLightTimeCorrectionFromCentralBody function for each gravitating body.
     *  \param transmitterState State of transmitter at transmission time.
     *  \param receiverState State of receiver at reception time
     *  \param transmissionTime Time of signal transmission
     *  \param receptionTime Time of signal reception
     *  \return Total light time correction due to gravitating masses defined by perturbingBodyStateFunctions_ and
     *  perturbingBodyGravitationalParameterFunctions_ member variables.
     */
    double calculateLightTimeCorrectionWithMultiLegLinkEndStates(
            const std::vector< Eigen::Vector6d >& linkEndsStates,
            const std::vector< double >& linkEndsTimes,
            const unsigned int currentMultiLegTransmitterIndex,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override;

    //! Function to compute the partial derivative of the light-time correction w.r.t. observation time
    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. observation time, equal to zero in this
     * model.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param fixedLinkEnd Reference link end for observation
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the time partial is to be taken
     * \return Light-time correction w.r.t. observation time
     */
    double calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings ) override
    {
        return 0.0;
    }

    //! Function to compute the partial derivative of the light-time correction w.r.t. link end position
    /*!
     * Function to compute the partial derivative of the light-time correction w.r.t. link end position
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
     * \return Light-time correction w.r.t. link end position
     */
    Eigen::Matrix< double, 3, 1 > calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
            const Eigen::Vector6d& transmitterState,
            const Eigen::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime,
            const LinkEndType linkEndAtWhichPartialIsEvaluated,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr ) override;

    //! Function to get the names of bodies causing light-time correction.
    /*!
     * Function to get the names of bodies causing light-time correction.
     * \return Names of bodies causing light-time correction.
     */
    std::vector< std::string > getPerturbingBodyNames( )
    {
        return perturbingBodyNames_;
    }

    //! Function to get the set of functions returning the gravitational parameters of the gravitating bodies.
    /*!
     * Function to get the set of functions returning the gravitational parameters of the gravitating bodies.
     * \return Set of functions returning the gravitational parameters of the gravitating bodies.
     */
    std::vector< std::function< double( ) > > getPerturbingBodyGravitationalParameterFunctions( )
    {
        return perturbingBodyGravitationalParameterFunctions_;
    }

    //! Function to get the total light-time correction, as computed by last call to calculateLightTimeCorrection.
    /*!
     * Function to get the total light-time correction, as computed by last call to calculateLightTimeCorrection.
     * \return Total light-time correction, as computed by last call to calculateLightTimeCorrection.
     */
    double getCurrentTotalLightTimeCorrection( )
    {
        return currentTotalLightTimeCorrection_;
    }

    //! Function to get the light-time correction of given perturbing body, as computed by last call to
    //! calculateLightTimeCorrection.
    /*!
     * Function to get the light-time correction of given perturbing body, as computed by last call to
     * calculateLightTimeCorrection.
     * \param bodyIndex Index in list of bodies for which the correction is to be returbed
     * \return Light-time correction of given perturbing body, as computed by last call to
     * calculateLightTimeCorrection.
     */
    double getCurrentLightTimeCorrectionComponent( const int bodyIndex )
    {
        return currentLighTimeCorrectionComponents_.at( bodyIndex );
    }

    //! Function to get the function returning the parametric post-Newtonian parameter gamma
    /*!
     * Function to get the function returning the parametric post-Newtonian parameter gamma
     * \return Function returning the parametric post-Newtonian parameter gamma
     */
    std::function< double( ) > getPpnParameterGammaFunction_( )
    {
        return ppnParameterGammaFunction_;
    }

    //! Function to get the bending flag
    /*!
     * Function to get the bending flag
     * \return Bending flag
     */
    bool getBendingFlag( )
    {
        return bending_;
    }

private:
    //! Set of function returning the state of the gravitating bodies as a function of time.
    std::vector< std::function< Eigen::Vector6d( const double ) > > perturbingBodyStateFunctions_;

    //! Set of functions returning the gravitational parameters of the gravitating bodies.
    std::vector< std::function< double( ) > > perturbingBodyGravitationalParameterFunctions_;

    //! Names of bodies causing light-time correction.
    std::vector< std::string > perturbingBodyNames_;

    //! Function returning the parametric post-Newtonian parameter gamma
    /*!
     *  Function returning the parametric post-Newtonian parameter gamma, a measure for the space-time curvature due to a
     *  unit rest mass (1.0 in GR)
     */
    std::function< double( ) > ppnParameterGammaFunction_;

    //! List of values (between 0 and 1) of how far into the light-time the state of each perturbing body is to be evaluated.
    std::vector< double > lightTimeEvaluationContribution_;

    //! List of light-time correction due to each separate perturbing body, as computed by last call to
    //! calculateLightTimeCorrection.
    std::vector< double > currentLighTimeCorrectionComponents_;

    //! Total light-time correction, as computed by last call to calculateLightTimeCorrection.
    double currentTotalLightTimeCorrection_;

    //! Boolean flag to determine if light bending should be included.
    bool bending_;
};

}  // namespace observation_models

}  // namespace tudat
#endif  // TUDAT_FIRSTORDERRELATIVISTICLIGHTTIMECORRECTION_H
