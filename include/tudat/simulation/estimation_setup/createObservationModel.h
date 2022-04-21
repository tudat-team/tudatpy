/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONMODEL_H
#define TUDAT_CREATEOBSERVATIONMODEL_H

#include <map>

#include <functional>
#include <boost/make_shared.hpp>


#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/observation_models/twoWayDopplerObservationModel.h"
#include "tudat/astro/observation_models/oneWayDifferencedRangeRateObservationModel.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/astro/observation_models/positionObservationModel.h"
#include "tudat/astro/observation_models/eulerAngleObservationModel.h"
#include "tudat/astro/observation_models/velocityObservationModel.h"
#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createLightTimeCalculator.h"
#include "tudat/simulation/estimation_setup/createObservationViability.h"


namespace tudat
{

namespace observation_models
{

//! Base class to define settings for creation of an observation bias model.
/*!
 *  Base class to define settings for creation of an observation bias model. For each specific bias type, a derived class
 *  is to be implemented, in which the specific properties of the bias model are given
 */
class ObservationBiasSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observationBiasType Type of bias model that is to be created.
     */
    ObservationBiasSettings(
            const observation_models::ObservationBiasTypes observationBiasType ):
        observationBiasType_( observationBiasType ){ }

    //! Destructor
    virtual ~ObservationBiasSettings( ){ }

    //! Type of bias model that is to be created.
    observation_models::ObservationBiasTypes observationBiasType_;
};

//! Class for defining settings for the creation of a multiple biases for a single observable
class MultipleObservationBiasSettings: public ObservationBiasSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param biasSettingsList List of settings for bias objects that are to be created.
     */
    MultipleObservationBiasSettings(
            const std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList ):
        ObservationBiasSettings( multiple_observation_biases ),
        biasSettingsList_( biasSettingsList ){ }

    //! Destructor
    ~MultipleObservationBiasSettings( ){ }

    //! List of settings for bias objects that are to be created.
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList_;
};

//! Class for defining settings for the creation of a constant absolute or relative observation bias model
class ConstantObservationBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param observationBias Constant bias that is to be added to the observable. The size of this vector must be equal to the
     * size of the observable to which it is assigned.
     * \param useAbsoluteBias Boolean to denote whether an absolute or relative bias is to be created.
     */
    ConstantObservationBiasSettings(
            const Eigen::VectorXd& observationBias,
            const bool useAbsoluteBias ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ? ( constant_absolute_bias ) : ( constant_relative_bias ) ),
        observationBias_( observationBias ), useAbsoluteBias_( useAbsoluteBias ){ }

    //! Destructor
    ~ConstantObservationBiasSettings( ){ }

    //! Constant bias that is to be added to the observable.
    /*!
     *  Constant bias that is to be added to the observable. The size of this vector must be equal to the
     *  size of the observable to which it is assigned.
     */
    Eigen::VectorXd observationBias_;

    //! Boolean to denote whether an absolute or relative bias is to be created.
    bool useAbsoluteBias_;
};

//! Class for defining settings for the creation of an arc-wise constant absolute or relative observation bias model
class ArcWiseConstantObservationBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases List of observation biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     * \param useAbsoluteBias Boolean to denote whether an absolute or relative bias is to be created.
     */
    ArcWiseConstantObservationBiasSettings(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::VectorXd >& observationBiases,
            const LinkEndType linkEndForTime,
            const bool useAbsoluteBias ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ?
                                     ( arc_wise_constant_absolute_bias ) : ( arc_wise_constant_relative_bias ) ),
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndForTime_( linkEndForTime ),
        useAbsoluteBias_( useAbsoluteBias ){ }

    //! Constuctor
    /*!
     * Constuctor
     * \param observationBiases Map of observation biases per arc, with bias as map value, and arc start time as map key
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     * \param useAbsoluteBias Boolean to denote whether an absolute or relative bias is to be created.
     */
    ArcWiseConstantObservationBiasSettings(
            const std::map< double, Eigen::VectorXd >& observationBiases,
            const LinkEndType linkEndForTime,
            const bool useAbsoluteBias  ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ?
                                     ( arc_wise_constant_absolute_bias ) : ( arc_wise_constant_relative_bias ) ),
        arcStartTimes_( utilities::createVectorFromMapKeys( observationBiases ) ),
        observationBiases_( utilities::createVectorFromMapValues( observationBiases ) ), linkEndForTime_( linkEndForTime ),
        useAbsoluteBias_( useAbsoluteBias ){ }

    //! Destructor
    ~ArcWiseConstantObservationBiasSettings( ){ }

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! List of observation biases per arc
    std::vector< Eigen::VectorXd > observationBiases_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

    //! Boolean to denote whether an absolute or relative bias is to be created.
    bool useAbsoluteBias_;
};

inline std::shared_ptr< ObservationBiasSettings > constantAbsoluteBias(
        const Eigen::VectorXd& observationBias )
{
    return std::make_shared< ConstantObservationBiasSettings >(
                observationBias, true );
}

inline std::shared_ptr< ObservationBiasSettings > constantRelativeBias(
        const Eigen::VectorXd& observationBias )
{
    return std::make_shared< ConstantObservationBiasSettings >(
                observationBias, false );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseAbsoluteBias(
        const std::vector< double >& arcStartTimes,
        const std::vector< Eigen::VectorXd >& observationBiases,
        const LinkEndType linkEndForTime)
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >(
                arcStartTimes, observationBiases, linkEndForTime, true );
}


inline std::shared_ptr< ObservationBiasSettings > arcWiseAbsoluteBias(
        const std::map< double, Eigen::VectorXd >& observationBiases,
        const LinkEndType linkEndForTime)
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >(
                observationBiases, linkEndForTime, true );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseRelativeBias(
        const std::vector< double >& arcStartTimes,
        const std::vector< Eigen::VectorXd >& observationBiases,
        const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >(
                arcStartTimes, observationBiases, linkEndForTime, false );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseRelativeBias(
        const std::map< double, Eigen::VectorXd >& observationBiases,
        const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >(
                observationBiases, linkEndForTime, false );
}

inline std::shared_ptr< ObservationBiasSettings > multipleObservationBiasSettings(
        const std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList  )
{
    return std::make_shared< MultipleObservationBiasSettings >(
                biasSettingsList );
}

//! Class used for defining the settings for an observation model that is to be created.
/*!
 * Class used for defining the settings for an observation model that is to be created. This class allows the type, light-time
 * corrections and bias for the observation to be set. For observation models that require additional information (e.g.
 * integration time, retransmission time, etc.), a specific derived class must be implemented.
 */
class ObservationModelSettings
{
public:


    //! Constructor
    /*!
     * Constructor (single light-time correction)
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
     * (nullptr if none)
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    ObservationModelSettings(
            const observation_models::ObservableType observableType,
            const LinkEnds linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        biasSettings_( biasSettings )
    {
        if( lightTimeCorrections != nullptr )
        {
            lightTimeCorrectionsList_.push_back( lightTimeCorrections );
        }
    }

    //! Constructor
    /*!
     * Constructor (multiple light-time correction)
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrectionsList List of settings for a single light-time correction that is to be used for the observation
     * model
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    ObservationModelSettings(
            const observation_models::ObservableType observableType,
            const LinkEnds linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        lightTimeCorrectionsList_( lightTimeCorrectionsList ),
        biasSettings_( biasSettings ){ }

    //! Destructor
    virtual ~ObservationModelSettings( ){ }

    //! Type of observation model that is to be created
    observation_models::ObservableType observableType_;

    LinkEnds linkEnds_;

    //! List of settings for a single light-time correction that is to be used for the observation model
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList_;

    //! Settings for the observation bias model that is to be used (default none: nullptr)
    std::shared_ptr< ObservationBiasSettings > biasSettings_;
};

std::vector< LinkEnds > getObservationModelListLinkEnds(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelList );



//! Enum defining all possible types of proper time rate computations in one-way Doppler
enum DopplerProperTimeRateType
{
    custom_doppler_proper_time_rate,
    direct_first_order_doppler_proper_time_rate
};

//! Base class to define the settings for proper time rate (at a single link end) in one-way Doppler mode.
class DopplerProperTimeRateSettings
{
public:
    DopplerProperTimeRateSettings( const DopplerProperTimeRateType dopplerProperTimeRateType ):
        dopplerProperTimeRateType_( dopplerProperTimeRateType ){ }

    virtual ~DopplerProperTimeRateSettings( ){ }

    DopplerProperTimeRateType dopplerProperTimeRateType_;
};

//! Class to define the settings for first-order, single body, proper time rate (at a single link end) in one-way Doppler mode.
class DirectFirstOrderDopplerProperTimeRateSettings: public DopplerProperTimeRateSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param centralBodyName Name of central body, fromw which the mass monopole is retrieved to compute the proper time rate,
     * and w.r.t. which the velocity of the point at which proper time rate is computed is taken
     */
    DirectFirstOrderDopplerProperTimeRateSettings(
            const std::string centralBodyName ):
        DopplerProperTimeRateSettings( direct_first_order_doppler_proper_time_rate ),
        centralBodyName_( centralBodyName ){ }

    //! Destructor.
    ~DirectFirstOrderDopplerProperTimeRateSettings( ){ }

    //! Name of central body
    /*!
     * Name of central body, fromw which the mass monopole is retrieved to compute the proper time rate,
     * and w.r.t. which the velocity of the point at which proper time rate is computed is taken
     */
    std::string centralBodyName_;
};

//! Class to define the settings for one-way Doppler observable
class OneWayDopplerObservationSettings: public ObservationModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
     * (nullptr if none)
     * \param transmitterProperTimeRateSettings Settings for proper time rate at transmitter
     * \param receiverProperTimeRateSettings Settings for proper time rate at receiver
     * \param biasSettings Settings for the observation bias model that is to be used (default none: NUL
     */
    OneWayDopplerObservationSettings(
            const LinkEnds& linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( one_way_doppler, linkEnds, lightTimeCorrections, biasSettings ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ){ }

    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrectionsList List of settings for a single light-time correction that is to be used for the observation
     * model (empty if none)
     * \param transmitterProperTimeRateSettings Settings for proper time rate at transmitter
     * \param receiverProperTimeRateSettings Settings for proper time rate at receiver
     * \param biasSettings Settings for the observation bias model that is to be used (default none: NUL
     */
    OneWayDopplerObservationSettings(
            const LinkEnds& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( one_way_doppler, linkEnds, lightTimeCorrectionsList, biasSettings ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ){ }

    //! Destructor
    ~OneWayDopplerObservationSettings( ){ }

    //! Settings for proper time rate at transmitter
    std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings_;

    //! Settings for proper time rate at receiver
    std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings_;
};



//! Class to define the settings for one-way Doppler observable
class TwoWayDopplerObservationSettings: public ObservationModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param uplinkOneWayDopplerSettings Settings for the one-way Doppler model of the uplink
     * \param downlinkOneWayDopplerSettings Settings for the one-way Doppler model of the downlink
     * \param biasSettings Settings for the observation bias model that is to be used (default none: NUL
     */
    TwoWayDopplerObservationSettings(
            const std::shared_ptr< OneWayDopplerObservationSettings > uplinkOneWayDopplerSettings,
            const std::shared_ptr< OneWayDopplerObservationSettings > downlinkOneWayDopplerSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( two_way_doppler, mergeUpDownLink(
                                      uplinkOneWayDopplerSettings->linkEnds_, downlinkOneWayDopplerSettings->linkEnds_ ),
                                  std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ),
        uplinkOneWayDopplerSettings_( uplinkOneWayDopplerSettings ),
        downlinkOneWayDopplerSettings_( downlinkOneWayDopplerSettings )
    {

    }

    TwoWayDopplerObservationSettings(
            const LinkEnds& linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( two_way_doppler, linkEnds,
                                  lightTimeCorrections, biasSettings )
    {
        uplinkOneWayDopplerSettings_ = std::make_shared< OneWayDopplerObservationSettings >(
                    getUplinkFromTwoWayLinkEnds( linkEnds ), lightTimeCorrections );
        downlinkOneWayDopplerSettings_ = std::make_shared< OneWayDopplerObservationSettings >(
                    getDownlinkFromTwoWayLinkEnds( linkEnds ), lightTimeCorrections );
    }


    //! Destructor
    ~TwoWayDopplerObservationSettings( ){ }

    //! Settings for the one-way Doppler model of the uplink
    std::shared_ptr< OneWayDopplerObservationSettings > uplinkOneWayDopplerSettings_;

    //! Settings for the one-way Doppler model of the downlink
    std::shared_ptr< OneWayDopplerObservationSettings > downlinkOneWayDopplerSettings_;
};




//! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler) observable
class OneWayDifferencedRangeRateObservationSettings: public ObservationModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param integrationTimeFunction Function that returns the integration time of observable as a function of time
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
     * (nullptr if none)
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    OneWayDifferencedRangeRateObservationSettings(
            const LinkEnds& linkEnds,
            const std::function< double( const double ) > integrationTimeFunction,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrections, biasSettings ),
        integrationTimeFunction_( integrationTimeFunction ){ }

    //! Constructor
    /*!
     * Constructor
     * \param integrationTimeFunction Function that returns the integration time of observable as a function of time
     * \param lightTimeCorrectionsList List of ettings for a single light-time correction that is to be used for the observation model
     * (empty if none)
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    OneWayDifferencedRangeRateObservationSettings(
            const LinkEnds& linkEnds,
            const std::function< double( const double ) > integrationTimeFunction,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings ),
        integrationTimeFunction_( integrationTimeFunction ){ }

    //! Destructor
    ~OneWayDifferencedRangeRateObservationSettings( ){ }

    //! Function that returns the integration time of observable as a function of time
    const std::function< double( const double ) > integrationTimeFunction_;

};


//! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler) observable
class NWayRangeObservationSettings: public ObservationModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param oneWayRangeObsevationSettings List of settings for one-way observables that make up n-way link (each must be for
     * one_way_range_
     * \param retransmissionTimesFunction Function that returns the retransmission delay time of the signal as a function of
     * observation timew.
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
            const std::function< std::vector< double >( const double ) > retransmissionTimesFunction =
            std::function< std::vector< double >( const double  ) >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( n_way_range, mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObsevationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ),
        oneWayRangeObsevationSettings_( oneWayRangeObsevationSettings ),
        retransmissionTimesFunction_( retransmissionTimesFunction ){ }

    //! Constructor
    /*!
     * Constructor for same light-time corrections per link
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
     * (nullptr if none)
     * \param numberOfLinkEnds Number of link ends in observable (equal to n+1 for 'n'-way observable)
     * \param retransmissionTimesFunction Function that returns the retransmission delay time of the signal as a function of
     * observation timew.
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationSettings(
            const LinkEnds& linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const int numberOfLinkEnds,
            const std::function< std::vector< double >( const double ) > retransmissionTimesFunction =
            std::function< std::vector< double >( const double  ) >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( n_way_range, linkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ),
        retransmissionTimesFunction_( retransmissionTimesFunction )
    {
        for( int i = 0; i < numberOfLinkEnds - 1; i++ )
        {
            oneWayRangeObsevationSettings_.push_back( std::make_shared< ObservationModelSettings >(
                                                          one_way_range, linkEnds, lightTimeCorrections ) );
        }
    }

    //! Destructor
    ~NWayRangeObservationSettings( ){ }

    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings_;

    //! Function that returns the integration time of observable as a function of time
    std::function< std::vector< double >( const double ) > retransmissionTimesFunction_;

};


inline std::shared_ptr< ObservationModelSettings > oneWayRangeSettings(
        const LinkEnds& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                one_way_range, linkEnds, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > angularPositionSettings(
        const LinkEnds& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                angular_position, linkEnds, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > positionObservableSettings(
        const LinkEnds& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                position_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > velocityObservableSettings(
        const LinkEnds& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                velocity_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > eulerAngle313ObservableSettings(
        const LinkEnds& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                euler_angle_313_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > oneWayOpenLoopDoppler(
        const LinkEnds& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr )
{
    return std::make_shared< OneWayDopplerObservationSettings >(
                linkEnds, lightTimeCorrectionsList, transmitterProperTimeRateSettings, receiverProperTimeRateSettings,
                biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > twoWayOpenLoopDoppler(
        const std::shared_ptr< OneWayDopplerObservationSettings > uplinkOneWayDopplerSettings,
        const std::shared_ptr< OneWayDopplerObservationSettings > downlinkOneWayDopplerSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< TwoWayDopplerObservationSettings >(
                uplinkOneWayDopplerSettings, downlinkOneWayDopplerSettings, biasSettings );
}


inline std::shared_ptr< ObservationModelSettings > twoWayOpenLoopDoppler(
        const LinkEnds& linkEnds,
        const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections = nullptr,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< TwoWayDopplerObservationSettings >(
                linkEnds, lightTimeCorrections, biasSettings );
}


inline std::shared_ptr< ObservationModelSettings > oneWayClosedLoopDoppler(
        const LinkEnds& linkEnds,
        const std::function< double( const double ) > integrationTimeFunction,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< OneWayDifferencedRangeRateObservationSettings >(
                linkEnds, integrationTimeFunction, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > oneWayClosedLoopDoppler(
        const LinkEnds& linkEnds,
        const double integrationTime,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< OneWayDifferencedRangeRateObservationSettings >(
                linkEnds, [=](const double){ return integrationTime; }, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > nWayRange(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::function< std::vector< double >( const double ) > retransmissionTimesFunction =
        std::function< std::vector< double >( const double  ) >( ))
{
    return std::make_shared< NWayRangeObservationSettings >(
                oneWayRangeObsevationSettings, retransmissionTimesFunction, biasSettings );
}


inline std::shared_ptr< ObservationModelSettings > nWayRangeSimple(
        const LinkEnds& linkEnds,
        const int numberOfLinkEnds,
        const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections = std::shared_ptr< LightTimeCorrectionSettings > ( ),
        const std::function< std::vector< double >( const double ) > retransmissionTimesFunction =
        std::function< std::vector< double >( const double  ) >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    // change order of input args from FF to accomodate default (empty) lightTimeCorrections
    return std::make_shared< NWayRangeObservationSettings >(
            linkEnds, lightTimeCorrections, numberOfLinkEnds, retransmissionTimesFunction, biasSettings );
}



//! Function to create the proper time rate calculator for use in one-way Doppler
/*!
 *  Function to create the proper time rate calculator for use in one-way Doppler
 *  \param properTimeRateSettings Settings for proper time rate model
 *  \param linkEnds Link ends of one-way Doppler observation  model
 *  \param bodies List of body objects that constitutes the environment
 *  \param linkEndForCalculator Link end for which the proper time rate is to be computed
 *  \return Proper time rate calculator for use in one-way Doppler
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< DopplerProperTimeRateInterface > createOneWayDopplerProperTimeCalculator(
        std::shared_ptr< DopplerProperTimeRateSettings > properTimeRateSettings,
        const LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies &bodies,
        const LinkEndType linkEndForCalculator )
{
    std::shared_ptr< DopplerProperTimeRateInterface > properTimeRateInterface;

    // Check tyope of proper time rate model
    switch( properTimeRateSettings->dopplerProperTimeRateType_ )
    {
    case direct_first_order_doppler_proper_time_rate:
    {
        // Check input consistency
        std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > directFirstOrderDopplerProperTimeRateSettings =
                std::dynamic_pointer_cast< DirectFirstOrderDopplerProperTimeRateSettings >( properTimeRateSettings );
        if( directFirstOrderDopplerProperTimeRateSettings == nullptr )
        {
            throw std::runtime_error( "Error when making DopplerProperTimeRateInterface, input type (direct_first_order_doppler_proper_time_rate) is inconsistent" );
        }
        else if( linkEnds.count( linkEndForCalculator ) == 0 )
        {
            std::string errorMessage = "Error when creating one-way Doppler proper time calculator, did not find link end " +
                    std::to_string( linkEndForCalculator );
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( bodies.at( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_ )->getGravityFieldModel( ) == nullptr )
            {
                throw std::runtime_error( "Error when making DirectFirstOrderDopplerProperTimeRateInterface, no gravity field found for " +
                                          directFirstOrderDopplerProperTimeRateSettings->centralBodyName_ );
            }
            else
            {
                // Retrieve gravitational parameter
                std::function< double( ) > gravitationalParameterFunction =
                        std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                   bodies.at( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_ )->
                                   getGravityFieldModel( ) );

                // Create calculation object.
                LinkEndId referencePointId =
                        std::make_pair( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, "" );
                if( ( linkEnds.at( receiver ) != referencePointId ) && ( linkEnds.at( transmitter ) != referencePointId ) )
                {
                    properTimeRateInterface = std::make_shared<
                            DirectFirstOrderDopplerProperTimeRateInterface >(
                                linkEndForCalculator, gravitationalParameterFunction,
                                directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, unidentified_link_end,
                                simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                                    std::make_pair( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, ""), bodies ) );
                }
                else
                {
                    throw std::runtime_error(
                                "Error, proper time reference point as link end not yet implemented for DopplerProperTimeRateInterface creation" );
                }
            }
        }
        break;
    }
    default:
        std::string errorMessage = "Error when creating one-way Doppler proper time calculator, did not recognize type " +
                std::to_string( properTimeRateSettings->dopplerProperTimeRateType_ );
        throw std::runtime_error( errorMessage );
    }
    return properTimeRateInterface;
}

////! Typedef of list of observation models per obserable type and link ends: note that ObservableType key must be consistent
////! with contents of ObservationModelSettings pointers. The ObservationSettingsMap may be used as well, which contains the same
////! type of information. This typedef, however, has some advantages in terms of book-keeping when creating observation models.
//typedef std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationModelSettings > > > SortedObservationSettingsMap;

////! Typedef of list of observation models per link ends. Multiple observation models for a single set of link ends are allowed,
////! since this typedef represents a multimap.
//typedef std::multimap< LinkEnds, std::shared_ptr< ObservationModelSettings > > ObservationSettingsMap;

//typedef std::vector< std::pair< LinkEnds, std::shared_ptr< ObservationModelSettings > > > ObservationSettingsVector;

//typedef std::map< LinkEnds, std::vector< std::shared_ptr< ObservationModelSettings > > > ObservationSettingsListPerLinkEnd;


std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortObservationModelSettingsByType(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings );

//! Function to create an object that computes an observation bias
/*!
 *  Function to create an object that computes an observation bias, which can represent any type of system-dependent influence
 *  on the observed value (e.g. absolute bias, relative bias, clock drift, etc.)
 *  \param linkEnds Observation link ends for which the bias is to be created.
 *  \param observableType Observable type for which bias is to be created
 *  \param biasSettings Settings for the observation bias that is to be created.
 *  \param bodies List of body objects that comprises the environment.
 *  \return Object that computes an observation bias according to requested settings.
 */
template< int ObservationSize = 1 >
std::shared_ptr< ObservationBias< ObservationSize > > createObservationBiasCalculator(
        const LinkEnds linkEnds,
        const ObservableType observableType,
        const std::shared_ptr< ObservationBiasSettings > biasSettings,
        const simulation_setup::SystemOfBodies &bodies )
{
    std::shared_ptr< ObservationBias< ObservationSize > > observationBias;
    switch( biasSettings->observationBiasType_ )
    {
    case constant_absolute_bias:
    {
        // Check input consistency
        std::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings = std::dynamic_pointer_cast<
                ConstantObservationBiasSettings >( biasSettings );
        if( constantBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making constant observation bias, settings are inconsistent" );
        }

        if( !constantBiasSettings->useAbsoluteBias_ )
        {
            throw std::runtime_error( "Error when making constant observation bias, class settings are inconsistent" );
        }

        // Check if size of bias is consistent with requested observable size
        if( constantBiasSettings->observationBias_.rows( ) != ObservationSize )
        {
            throw std::runtime_error( "Error when making constant observation bias, bias size is inconsistent" );
        }
        observationBias = std::make_shared< ConstantObservationBias< ObservationSize > >(
                    constantBiasSettings->observationBias_ );
        break;
    }
    case arc_wise_constant_absolute_bias:
    {
        // Check input consistency
        std::shared_ptr< ArcWiseConstantObservationBiasSettings > arcwiseBiasSettings = std::dynamic_pointer_cast<
                ArcWiseConstantObservationBiasSettings >( biasSettings );
        if( arcwiseBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making arc-wise observation bias, settings are inconsistent" );
        }
        else if( !arcwiseBiasSettings->useAbsoluteBias_ )
        {
            throw std::runtime_error( "Error when making arc-wise observation bias, class contents are inconsistent" );
        }

        std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases;
        for( unsigned int i = 0; i < arcwiseBiasSettings->observationBiases_.size( ); i++ )
        {
            // Check if size of bias is consistent with requested observable size
            if( arcwiseBiasSettings->observationBiases_.at( i ).rows( ) != ObservationSize )
            {
                throw std::runtime_error( "Error when making arc-wise observation bias, bias size is inconsistent" );
            }
            else
            {
                observationBiases.push_back( arcwiseBiasSettings->observationBiases_.at( i ) );
            }
        }
        observationBias = std::make_shared< ConstantArcWiseObservationBias< ObservationSize > >(
                    arcwiseBiasSettings->arcStartTimes_, observationBiases,
                    observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                        observableType, arcwiseBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ) );
        break;
    }
    case constant_relative_bias:
    {
        // Check input consistency
        std::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings = std::dynamic_pointer_cast<
                ConstantObservationBiasSettings >( biasSettings );
        if( constantBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making constant relative observation bias, settings are inconsistent" );
        }

        if( constantBiasSettings->useAbsoluteBias_ )
        {
            throw std::runtime_error( "Error when making constant relative observation bias, class settings are inconsistent" );
        }

        // Check if size of bias is consistent with requested observable size
        if( constantBiasSettings->observationBias_.rows( ) != ObservationSize )
        {
            throw std::runtime_error( "Error when making constant relative observation bias, bias size is inconsistent" );
        }
        observationBias = std::make_shared< ConstantRelativeObservationBias< ObservationSize > >(
                    constantBiasSettings->observationBias_ );
        break;
    }
    case arc_wise_constant_relative_bias:
    {
        // Check input consistency
        std::shared_ptr< ArcWiseConstantObservationBiasSettings > arcwiseBiasSettings = std::dynamic_pointer_cast<
                ArcWiseConstantObservationBiasSettings >( biasSettings );
        if( arcwiseBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making arc-wise relative observation bias, settings are inconsistent" );
        }
        else if( arcwiseBiasSettings->useAbsoluteBias_ )
        {
            throw std::runtime_error( "Error when making arc-wise relative observation bias, class contents are inconsistent" );
        }

        std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases;
        for( unsigned int i = 0; i < arcwiseBiasSettings->observationBiases_.size( ); i++ )
        {
            // Check if size of bias is consistent with requested observable size
            if( arcwiseBiasSettings->observationBiases_.at( i ).rows( ) != ObservationSize )
            {
                throw std::runtime_error( "Error when making arc-wise observation bias, bias size is inconsistent" );
            }
            else
            {
                observationBiases.push_back( arcwiseBiasSettings->observationBiases_.at( i ) );
            }
        }
        observationBias = std::make_shared< ConstantRelativeArcWiseObservationBias< ObservationSize > >(
                    arcwiseBiasSettings->arcStartTimes_, observationBiases,
                    observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                        observableType, arcwiseBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ) );
        break;
    }
    case multiple_observation_biases:
    {
        // Check input consistency
        std::shared_ptr< MultipleObservationBiasSettings > multiBiasSettings = std::dynamic_pointer_cast<
                MultipleObservationBiasSettings >( biasSettings );
        if( multiBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making multiple observation biases, settings are inconsistent" );
        }

        // Create list of biases
        std::vector< std::shared_ptr< ObservationBias< ObservationSize > > > observationBiasList;
        for( unsigned int i = 0; i < multiBiasSettings->biasSettingsList_.size( ); i++ )
        {
            observationBiasList.push_back( createObservationBiasCalculator< ObservationSize >(
                                               linkEnds, observableType, multiBiasSettings->biasSettingsList_.at( i ) , bodies ) );
        }

        // Create combined bias object
        observationBias = std::make_shared< MultiTypeObservationBias< ObservationSize > >(
                    observationBiasList );
        break;
    }
    default:
    {
        std::string errorMessage = "Error when making observation bias, bias type " +
                std::to_string( biasSettings->observationBiasType_  ) + " not recognized";
        throw std::runtime_error( errorMessage );
    }
    }
    return observationBias;
}

//! Interface class for creating observation models
/*!
 *  Interface class for creating observation models. This class is used instead of a single templated free function to
 *  allow ObservationModel deroved classed with different ObservationSize template arguments to be created using the same
 *  interface. This class has template specializations for each value of ObservationSize, and contains a single
 *  createObservationModel function that performs the required operation.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
class ObservationModelCreator
{
public:

    //! Function to create an observation model.
    /*!
     * Function to create an observation model.
     * \param linkEnds Link ends for observation model that is to be created
     * \param observationSettings Settings for observation model that is to be created.
     * \param bodies List of body objects that comprises the environment
     * \return Observation model of required settings.
     */
    static std::shared_ptr< observation_models::ObservationModel<
    ObservationSize, ObservationScalarType, TimeType > > createObservationModel(
            const std::shared_ptr< ObservationModelSettings > observationSettings,
            const simulation_setup::SystemOfBodies &bodies );
};

//! Interface class for creating observation models of size 1.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 1, ObservationScalarType, TimeType >
{
public:

    //! Function to create an observation model of size 1.
    /*!
     * Function to create an observation model of size 1.
     * \param linkEnds Link ends for observation model that is to be created
     * \param observationSettings Settings for observation model that is to be created (must be for observation model if size 1).
     * \param bodies List of body objects that comprises the environment
     * \return Observation model of required settings.
     */
    static std::shared_ptr< observation_models::ObservationModel<
    1, ObservationScalarType, TimeType > > createObservationModel(
            const std::shared_ptr< ObservationModelSettings > observationSettings,
            const simulation_setup::SystemOfBodies &bodies )
    {
        using namespace observation_models;

        std::shared_ptr< observation_models::ObservationModel<
                1, ObservationScalarType, TimeType > > observationModel;
        LinkEnds linkEnds = observationSettings->linkEnds_;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
        {
        case one_way_range:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making 1 way range model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no transmitter found" );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model
            observationModel = std::make_shared< OneWayRangeObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds, createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodies, observationSettings->lightTimeCorrectionsList_ ),
                        observationBias );

            break;
        }
        case one_way_doppler:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making 1 way Doppler model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way Doppler model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way Doppler model, no transmitter found" );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            if( std::dynamic_pointer_cast< OneWayDopplerObservationSettings >( observationSettings ) == nullptr )
            {
                // Create observation model
                observationModel = std::make_shared< OneWayDopplerObservationModel<
                        ObservationScalarType, TimeType > >(
                            linkEnds,
                            createLightTimeCalculator< ObservationScalarType, TimeType >(
                                linkEnds.at( transmitter ), linkEnds.at( receiver ),
                                bodies, observationSettings->lightTimeCorrectionsList_ ),
                            observationBias );
            }
            else
            {
                std::shared_ptr< OneWayDopplerObservationSettings > oneWayDopplerSettings =
                        std::dynamic_pointer_cast< OneWayDopplerObservationSettings >( observationSettings );
                std::shared_ptr< DopplerProperTimeRateInterface > transmitterProperTimeRate = nullptr;
                std::shared_ptr< DopplerProperTimeRateInterface > receiverProperTimeRate = nullptr;
                if( oneWayDopplerSettings->transmitterProperTimeRateSettings_ != nullptr )
                {
                    transmitterProperTimeRate =
                            createOneWayDopplerProperTimeCalculator< ObservationScalarType, TimeType >(
                                oneWayDopplerSettings->transmitterProperTimeRateSettings_, linkEnds, bodies, transmitter );
                }

                if( oneWayDopplerSettings->transmitterProperTimeRateSettings_ != nullptr )
                {
                    receiverProperTimeRate =
                            createOneWayDopplerProperTimeCalculator< ObservationScalarType, TimeType >(
                                oneWayDopplerSettings->receiverProperTimeRateSettings_, linkEnds, bodies, receiver );
                }

                // Create observation model
                observationModel = std::make_shared< OneWayDopplerObservationModel<
                        ObservationScalarType, TimeType > >(
                            linkEnds,
                            createLightTimeCalculator< ObservationScalarType, TimeType >(
                                linkEnds.at( transmitter ), linkEnds.at( receiver ),
                                bodies, observationSettings->lightTimeCorrectionsList_ ),
                                transmitterProperTimeRate,
                                receiverProperTimeRate,
                            observationBias );
            }

            break;
        }

        case two_way_doppler:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 3 )
            {
                std::string errorMessage =
                        "Error when making 2 way Doppler model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making 2 way Doppler model, no receiver found" );
            }

            if( linkEnds.count( reflector1 ) == 0 )
            {
                throw std::runtime_error( "Error when making 2 way Doppler model, no retransmitter found" );
            }

            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making 2 way Doppler model, no transmitter found" );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model

            LinkEnds uplinkLinkEnds;
            uplinkLinkEnds[ transmitter ] = linkEnds.at( transmitter );
            uplinkLinkEnds[ receiver ] = linkEnds.at( reflector1 );

            LinkEnds downlinkLinkEnds;
            downlinkLinkEnds[ transmitter ] = linkEnds.at( reflector1 );
            downlinkLinkEnds[ receiver ] = linkEnds.at( receiver );

            std::shared_ptr< TwoWayDopplerObservationSettings > twoWayDopplerSettings =
                    std::dynamic_pointer_cast< TwoWayDopplerObservationSettings >( observationSettings );

            if( twoWayDopplerSettings == nullptr )
            {
                observationModel = std::make_shared< TwoWayDopplerObservationModel<
                        ObservationScalarType, TimeType > >(
                            linkEnds,
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    std::make_shared< ObservationModelSettings >(
                                        one_way_doppler, uplinkLinkEnds, observationSettings->lightTimeCorrectionsList_ ), bodies ) ),
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    std::make_shared< ObservationModelSettings >(
                                        one_way_doppler, downlinkLinkEnds, observationSettings->lightTimeCorrectionsList_ ), bodies ) ),
                            observationBias );
            }
            else
            {
                observationModel = std::make_shared< TwoWayDopplerObservationModel<
                        ObservationScalarType, TimeType > >(
                            linkEnds,
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                     twoWayDopplerSettings->uplinkOneWayDopplerSettings_, bodies ) ),
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    twoWayDopplerSettings->downlinkOneWayDopplerSettings_, bodies ) ),
                            observationBias );
            }

            break;
        }

        case one_way_differenced_range:
        {
            std::shared_ptr< OneWayDifferencedRangeRateObservationSettings > rangeRateObservationSettings =
                    std::dynamic_pointer_cast< OneWayDifferencedRangeRateObservationSettings >( observationSettings );
            if( rangeRateObservationSettings == nullptr )
            {
                throw std::runtime_error( "Error when making differenced one-way range rate, input type is inconsistent" );
            }
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making 1 way range model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making 1 way range model, no transmitter found" );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model
            observationModel = std::make_shared< OneWayDifferencedRangeObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodies, observationSettings->lightTimeCorrectionsList_ ),
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodies, observationSettings->lightTimeCorrectionsList_ ),
                        rangeRateObservationSettings->integrationTimeFunction_,
                        observationBias );

            break;
        }
        case n_way_range:
        {
            // Check consistency input.
            if( linkEnds.size( ) < 2 )
            {
                std::string errorMessage =
                        "Error when making n way range model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making n way range model, no receiver found" );
            }

            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making n way range model, no transmitter found" );
            }

            // Check link end consistency.
            for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
                 linkEndIterator++ )
            {
                if( ( linkEndIterator->first != transmitter ) && ( linkEndIterator->first != receiver ) )
                {
                    int linkEndIndex = static_cast< int >( linkEndIterator->first );
                    LinkEndType previousLinkEndType = static_cast< LinkEndType >( linkEndIndex - 1 );

                    if( linkEnds.count( previousLinkEndType ) == 0 )
                    {
                        throw std::runtime_error( "Error when making n-way range model, did not find link end type " +
                                                  std::to_string( previousLinkEndType ) );
                    }
                }
            }

            // Create observation bias object
            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList;

            std::function< std::vector< double >( const double ) > retransmissionTimesFunction_;

            std::shared_ptr< NWayRangeObservationSettings > nWayRangeObservationSettings =
                    std::dynamic_pointer_cast< NWayRangeObservationSettings >( observationSettings );

            if( nWayRangeObservationSettings == nullptr )
            {
                lightTimeCorrectionsList = observationSettings->lightTimeCorrectionsList_;
            }
            else if( nWayRangeObservationSettings->oneWayRangeObsevationSettings_.size( ) != linkEnds.size( ) - 1 )
            {
                throw std::runtime_error( "Error whaen making n-way range, input data is inconsistent" );
            }
            else
            {
                retransmissionTimesFunction_ = nWayRangeObservationSettings->retransmissionTimesFunction_;
            }

            // Define light-time calculator list
            std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > lightTimeCalculators;

            // Iterate over all link ends and create light-time calculators
            LinkEnds::const_iterator transmitterIterator = linkEnds.begin( );
            LinkEnds::const_iterator receiverIterator = linkEnds.begin( );
            receiverIterator++;
            for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
            {
                if( nWayRangeObservationSettings != nullptr )
                {
                    if( nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->observableType_ != one_way_range )
                    {
                        throw std::runtime_error( "Error in n-way observable creation, consituent link is not of type 1-way" );
                    }
                    lightTimeCalculators.push_back(
                                createLightTimeCalculator< ObservationScalarType, TimeType >(
                                    transmitterIterator->second, receiverIterator->second,
                                    bodies, nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->
                                    lightTimeCorrectionsList_ ) );
                }
                else
                {
                    lightTimeCalculators.push_back(
                                createLightTimeCalculator< ObservationScalarType, TimeType >(
                                    transmitterIterator->second, receiverIterator->second,
                                    bodies, observationSettings->lightTimeCorrectionsList_ ) );
                }

                transmitterIterator++;
                receiverIterator++;
            }

            // Create observation model
            observationModel = std::make_shared< NWayRangeObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        lightTimeCalculators, retransmissionTimesFunction_,
                        observationBias );
            break;
        }

        default:
            std::string errorMessage = "Error, observable " + std::to_string(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 1 observation model.";
            throw std::runtime_error( errorMessage );
        }
        return observationModel;
    }

};

//! Interface class for creating observation models of size 2.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 2, ObservationScalarType, TimeType >
{
public:

    //! Function to create an observation model of size 2.
    /*!
     * Function to create an observation model of size 2.
     * \param linkEnds Link ends for observation model that is to be created
     * \param observationSettings Settings for observation model that is to be created (must be for observation model if size 1).
     * \param bodies List of body objects that comprises the environment
     * \return Observation model of required settings.
     */
    static std::shared_ptr< observation_models::ObservationModel<
    2, ObservationScalarType, TimeType > > createObservationModel(
            const std::shared_ptr< ObservationModelSettings > observationSettings,
            const simulation_setup::SystemOfBodies &bodies )
    {
        using namespace observation_models;
        std::shared_ptr< observation_models::ObservationModel<
                2, ObservationScalarType, TimeType > > observationModel;
        LinkEnds linkEnds = observationSettings->linkEnds_;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
        {
        case angular_position:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage =
                        "Error when making angular position model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }
            if( linkEnds.count( receiver ) == 0 )
            {
                throw std::runtime_error( "Error when making angular position model, no receiver found" );
            }
            if( linkEnds.count( transmitter ) == 0 )
            {
                throw std::runtime_error( "Error when making angular position model, no transmitter found" );
            }


            std::shared_ptr< ObservationBias< 2 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator< 2 >(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model
            observationModel = std::make_shared< AngularPositionObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds.at( transmitter ), linkEnds.at( receiver ),
                            bodies, observationSettings->lightTimeCorrectionsList_ ),
                        observationBias );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + std::to_string(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 2 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }

        return observationModel;
    }

};

//! Interface class for creating observation models of size 3.
template< typename ObservationScalarType, typename TimeType >
class ObservationModelCreator< 3, ObservationScalarType, TimeType >
{
public:

    //! Function to create an observation model of size 3.
    /*!
     * Function to create an observation model of size 3.
     * \param linkEnds Link ends for observation model that is to be created
     * \param observationSettings Settings for observation model that is to be created (must be for observation model if size 1).
     * \param bodies List of body objects that comprises the environment
     * \return Observation model of required settings.
     */
    static std::shared_ptr< observation_models::ObservationModel<
    3, ObservationScalarType, TimeType > > createObservationModel(
            const std::shared_ptr< ObservationModelSettings > observationSettings,
            const simulation_setup::SystemOfBodies &bodies )
    {
        using namespace observation_models;
        std::shared_ptr< observation_models::ObservationModel<
                3, ObservationScalarType, TimeType > > observationModel;

        LinkEnds linkEnds = observationSettings->linkEnds_;

        // Check type of observation model.
        switch( observationSettings->observableType_ )
        {
        case position_observable:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 1 )
            {
                std::string errorMessage =
                        "Error when making position observable model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                throw std::runtime_error( "Error when making position observable model, no observed_body found" );
            }

            if( observationSettings->lightTimeCorrectionsList_.size( ) > 0 )
            {
                throw std::runtime_error( "Error when making position observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).second != "" )
            {
                throw std::runtime_error( "Error, cannot yet create position function for reference point" );
            }

            std::shared_ptr< ObservationBias< 3 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator< 3 >(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }


            // Create observation model
            observationModel = std::make_shared< PositionObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris<
                                   ObservationScalarType, TimeType >,
                                   bodies.at( linkEnds.at( observed_body ).first ), std::placeholders::_1 ),
                        observationBias );

            break;
        }
        case euler_angle_313_observable:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 1 )
            {
                std::string errorMessage =
                        "Error when making euler angle observable model, " +
                        std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                throw std::runtime_error( "Error when making euler angle observable model, no observed_body found" );
            }

            if( observationSettings->lightTimeCorrectionsList_.size( ) > 0 )
            {
                throw std::runtime_error( "Error when making euler angle observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).second != "" )
            {
                throw std::runtime_error( "Error, cannot yet create euler angle function for reference point" );
            }

            std::shared_ptr< ObservationBias< 3 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator< 3 >(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            std::function< Eigen::Quaterniond( const TimeType ) > toBodyFixedFrameFunction;
            if( bodies.at( linkEnds.at( observed_body ).first )->getRotationalEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error, cannot euler angle observable; no rotation model found" );
            }
            else
            {
                toBodyFixedFrameFunction = std::bind(
                            &ephemerides::RotationalEphemeris::getRotationToTargetFrameTemplated< TimeType >,
                            bodies.at( linkEnds.at( observed_body ).first )->getRotationalEphemeris( ),
                            std::placeholders::_1 );
            }


            // Create observation model
            observationModel = std::make_shared< EulerAngle313ObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        toBodyFixedFrameFunction, observationBias );

            break;
        }
        case velocity_observable:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 1 )
            {
                std::string errorMessage =
                        "Error when making velocity observable model, " + std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                throw std::runtime_error( "Error when making velocity observable model, no observed_body found" );
            }

            if( observationSettings->lightTimeCorrectionsList_.size( ) > 0 )
            {

                throw std::runtime_error( "Error when making velocity observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).second != "" )
            {
                throw std::runtime_error( "Error, cannot yet create velocity function for reference point" );
            }

            std::shared_ptr< ObservationBias< 3 > > observationBias;
            if( observationSettings->biasSettings_ != NULL )
            {
                observationBias =
                        createObservationBiasCalculator< 3 >(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }


            // Create observation model
            observationModel = std::make_shared< VelocityObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds,
                        std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris<
                                   ObservationScalarType, TimeType >,
                                   bodies.at( linkEnds.at( observed_body ).first ), std::placeholders::_1 ),
                        observationBias );

            break;
        }
        default:
            std::string errorMessage = "Error, observable " + std::to_string(
                        observationSettings->observableType_ ) +
                    "  not recognized when making size 3 observation model.";
            throw std::runtime_error( errorMessage );
            break;
        }
        return observationModel;
    }
};

//! Function to create an object to simulate observations of a given type
/*!
 *  Function to create an object to simulate observations of a given type
 *  \param observableType Type of observable for which object is to simulate ObservationSimulator
 *  \param settingsPerLinkEnds Map of settings for the observation models that are to be created in the simulator object: one
 *  for each required set of link ends (each settings object must be consistent with observableType).
 *  \param bodies Map of Body objects that comprise the environment
 *  \return Object that simulates the observables according to the provided settings.
 */
template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > createObservationSimulator(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings  > > settingsList,
        const simulation_setup::SystemOfBodies &bodies )
{
    std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
            observationModels;

    // Iterate over all link ends
    for( unsigned int i = 0; i < settingsList.size( ); i++ )
    {
        observationModels[ settingsList.at( i )->linkEnds_ ] = ObservationModelCreator<
                ObservationSize, ObservationScalarType, TimeType >::createObservationModel(
                    settingsList.at( i ), bodies );
    }

    return std::make_shared< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationModels );
}

//! Function to create a map of object to simulate observations (one object for each type of observable).
/*!
 *  Function to create a map of object to simulate observations (one object for each type of observable).
 *  \param observationSettingsList List of settings for the observation models that are to be created in the simulator object
 *  \param bodies Map of Body objects that comprise the environment
 *  \return List of objects that simulate the observables according to the provided settings.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > createObservationSimulators(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies )
{
    std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators;
    std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationSettingsList =
            sortObservationModelSettingsByType( observationSettingsList );

    // Iterate over all observables
    for( auto it : sortedObservationSettingsList )
    {
        // Call createObservationSimulator of required observation size
        ObservableType observableType = it.first;
        int observableSize = getObservableSize( observableType );
        switch( observableSize )
        {
        case 1:
        {
            observationSimulators.push_back( createObservationSimulator< 1, ObservationScalarType, TimeType >(
                        observableType, it.second, bodies ) );
            break;
        }
        case 2:
        {
            observationSimulators.push_back( createObservationSimulator< 2, ObservationScalarType, TimeType >(
                        observableType, it.second, bodies ) );
            break;
        }
        case 3:
        {
            observationSimulators.push_back( createObservationSimulator< 3, ObservationScalarType, TimeType >(
                        observableType, it.second, bodies ) );
            break;
        }
        default:
            throw std::runtime_error( "Error, cannot create observation simulator for size other than 1,2 and 3 ");
        }
    }
    return observationSimulators;
}




} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATEOBSERVATIONMODEL_H
