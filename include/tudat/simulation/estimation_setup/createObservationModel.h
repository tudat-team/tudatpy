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



#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"
#include "tudat/astro/observation_models/nWayRangeRateObservationModel.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/observation_models/twoWayDopplerObservationModel.h"
#include "tudat/astro/observation_models/oneWayDifferencedRangeRateObservationModel.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/astro/observation_models/relativeAngularPositionObservationModel.h"
#include "tudat/astro/observation_models/positionObservationModel.h"
#include "tudat/astro/observation_models/relativePositionObservationModel.h"
#include "tudat/astro/observation_models/eulerAngleObservationModel.h"
#include "tudat/astro/observation_models/velocityObservationModel.h"
#include "tudat/astro/observation_models/observationSimulator.h"
#include "tudat/astro/observation_models/dsnNWayAveragedDopplerObservationModel.h"
#include "tudat/astro/observation_models/dopplerMeasuredFrequencyObservationModel.h"
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

//! Class for defining settings for the creation of a constant time drift model
class ConstantTimeDriftBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param observationBias Constant time drift bias that is to be considered for the observation time. The size of this vector must be equal to the
     * size of the observable to which it is assigned.
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     * \param referenceEpoch Reference epoch at which the effect of the time drift is initialised.
     */
    ConstantTimeDriftBiasSettings(
            const Eigen::VectorXd& timeDriftBias,
            const LinkEndType linkEndForTime,
            const double referenceEpoch ):
            ObservationBiasSettings( constant_time_drift_bias ),
            timeDriftBias_( timeDriftBias ), linkEndForTime_( linkEndForTime ), referenceEpoch_( referenceEpoch ){ }

    //! Destructor
    ~ConstantTimeDriftBiasSettings( ){ }

    //! Constant time drift bias that is to be considered for the observation time.
    /*!
     *  Constant time drift bias that is to be considered for the observation time. The size of this vector must be equal to the
     *  size of the observable to which it is assigned.
     */
    Eigen::VectorXd timeDriftBias_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

    //! Reference epoch at which the effect of the time drift is supposed to be zero
    double referenceEpoch_;
};

//! Class for defining settings for the creation of an arc-wise time drift bias model
class ArcWiseTimeDriftBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (timeDriftBiases) are used
     * \param timeDriftBiases List of time drift biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     * \param referenceEpochs Reference epochs (per arc) at which the time drifts are initialised.
     */
    ArcWiseTimeDriftBiasSettings(
            const std::vector< double >& arcStartTimes,
            const std::vector< Eigen::VectorXd >& timeDriftBiases,
            const LinkEndType linkEndForTime,
            const std::vector< double > referenceEpochs ):
            ObservationBiasSettings( arc_wise_time_drift_bias ),
            arcStartTimes_( arcStartTimes ), timeDriftBiases_( timeDriftBiases ), linkEndForTime_( linkEndForTime ),
            referenceEpochs_( referenceEpochs ){ }

    //! Constuctor
    /*!
     * Constuctor
     * \param timeDriftBiases Map of observation biases per arc, with bias as map value, and arc start time as map key
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     * \param referenceEpochs Reference epochs (per arc) at which the time drifts are initialised.
     */
    ArcWiseTimeDriftBiasSettings(
            const std::map< double, Eigen::VectorXd >& timeDriftBiases,
            const LinkEndType linkEndForTime,
            const std::vector< double > referenceEpochs ):
            ObservationBiasSettings( arc_wise_time_drift_bias ),
            arcStartTimes_( utilities::createVectorFromMapKeys( timeDriftBiases ) ),
            timeDriftBiases_( utilities::createVectorFromMapValues( timeDriftBiases ) ), linkEndForTime_( linkEndForTime ),
            referenceEpochs_( referenceEpochs ){ }

    //! Destructor
    ~ArcWiseTimeDriftBiasSettings( ){ }

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! List of time drift biases per arc
    std::vector< Eigen::VectorXd > timeDriftBiases_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

    //! Reference epochs at which the time drifts are initialised.
    std::vector< double > referenceEpochs_;
};

//! Class for defining settings for the creation of a constant time bias
class ConstantTimeBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param timeBias Constant time bias that is to be considered for the observation time. The size of this vector must be equal to the
     * size of the observable to which it is assigned.
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     */
    ConstantTimeBiasSettings(
            const double timeBias,
            const LinkEndType linkEndForTime ):
            ObservationBiasSettings( constant_time_bias ),
            timeBias_( timeBias ), linkEndForTime_( linkEndForTime ){ }

    //! Destructor
    ~ConstantTimeBiasSettings( ){ }

    //! Constant time bias that is to be considered for the observation time.
    /*!
     *  Constant time bias that is to be considered for the observation time.
     */
    double timeBias_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

};

//! Class for defining settings for the creation of an arc-wise time bias model
class ArcWiseTimeBiasSettings: public ObservationBiasSettings
{
public:

    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (timeBiases) are used
     * \param timeBiases List of time biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     */
    ArcWiseTimeBiasSettings(
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& timeBiases,
            const LinkEndType linkEndForTime ):
            ObservationBiasSettings( arc_wise_time_bias ),
            arcStartTimes_( arcStartTimes ), timeBiases_( timeBiases ), linkEndForTime_( linkEndForTime ){ }

    //! Constuctor
    /*!
     * Constuctor
     * \param timeBiases Map of observation biases per arc, with bias as map value, and arc start time as map key
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time (and current arc)
     */
    ArcWiseTimeBiasSettings(
            const std::map< double, double >& timeBiases,
            const LinkEndType linkEndForTime ):
            ObservationBiasSettings( arc_wise_time_bias ),
            arcStartTimes_( utilities::createVectorFromMapKeys( timeBiases ) ),
            timeBiases_( utilities::createVectorFromMapValues( timeBiases ) ), linkEndForTime_( linkEndForTime ){ }

    //! Destructor
    ~ArcWiseTimeBiasSettings( ){ }

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! List of time biases per arc
    std::vector< double > timeBiases_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

};


//! Class for defining settings for the creation of a multiple biases for a single observable
class TiminigSystemBiasSettings: public ObservationBiasSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param biasSettingsList List of settings for bias objects that are to be created.
     */
    TiminigSystemBiasSettings(
            const std::string& bodyName,
            const std::string& stationName ):
            ObservationBiasSettings( clock_induced_bias ),
            bodyName_( bodyName ), stationName_( stationName ){ }

    //! Destructor
    ~TiminigSystemBiasSettings( ){ }

    const std::string bodyName_;

    const std::string stationName_;
};

inline std::shared_ptr< ObservationBiasSettings > clockInducedBias(
        const std::string& bodyName, const std::string& stationName  )
{
    return std::make_shared< TiminigSystemBiasSettings >(
            bodyName, stationName );
}

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

inline std::shared_ptr< ObservationBiasSettings > constantTimeDriftBias(
        const Eigen::VectorXd& timeDriftBias,
        const LinkEndType linkEndForTime,
        const double referenceEpoch )
{
    return std::make_shared< ConstantTimeDriftBiasSettings >( timeDriftBias, linkEndForTime, referenceEpoch );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeDriftBias(
        const std::vector< Eigen::VectorXd >& timeDriftBiases,
        const std::vector< double >& arcStartTimes,
        const LinkEndType linkEndForTime,
        const std::vector< double >& referenceEpochs )
{
    return std::make_shared< ArcWiseTimeDriftBiasSettings >( arcStartTimes, timeDriftBiases, linkEndForTime, referenceEpochs );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeDriftBias(
        const std::map< double, Eigen::VectorXd >& timeDriftBiases,
        const LinkEndType linkEndForTime,
        const std::vector< double > referenceEpochs )
{
    return std::make_shared< ArcWiseTimeDriftBiasSettings >( timeDriftBiases, linkEndForTime, referenceEpochs );
}

inline std::shared_ptr< ObservationBiasSettings > constantTimeBias(
        const double timeBias,
        const LinkEndType linkEndForTime )
{
    return std::make_shared< ConstantTimeBiasSettings >( timeBias, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias(
    const std::vector< double >& arcStartTimes,
    const LinkEndType linkEndForTime )
{
    std::vector< double > timeBiases = std::vector< double >( arcStartTimes.size( ), 0.0 );
    return std::make_shared< ArcWiseTimeBiasSettings >( arcStartTimes, timeBiases, linkEndForTime );
}


inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias(
        const std::vector< double >& timeBiases,
        const std::vector< double >& arcStartTimes,
        const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseTimeBiasSettings >( arcStartTimes, timeBiases, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias(
        const std::map< double, double >& timeBiases,
        const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseTimeBiasSettings >( timeBiases, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > twoWayTimeScaleRangeBias(  )
{
    return std::make_shared< ObservationBiasSettings >( two_way_range_time_scale_bias );
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
            const LinkDefinition linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        biasSettings_( biasSettings ),
        lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
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
            const LinkDefinition linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        observableType_( observableType ),
        linkEnds_( linkEnds ),
        lightTimeCorrectionsList_( lightTimeCorrectionsList ),
        biasSettings_( biasSettings ),
    lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ){ }

    //! Destructor
    virtual ~ObservationModelSettings( ){ }

    //! Type of observation model that is to be created
    observation_models::ObservableType observableType_;

    LinkDefinition linkEnds_;

    //! List of settings for a single light-time correction that is to be used for the observation model
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList_;

    //! Settings for the observation bias model that is to be used (default none: nullptr)
    std::shared_ptr< ObservationBiasSettings > biasSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria_;
};

std::vector< LinkDefinition > getObservationModelListLinkEnds(
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
            const LinkDefinition& linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings(
            one_way_doppler, linkEnds, lightTimeCorrections, biasSettings, lightTimeConvergenceCriteria ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ),
    normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight ){ }

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
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings( one_way_doppler, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ),
        normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight ){ }

    //! Destructor
    ~OneWayDopplerObservationSettings( ){ }

    //! Settings for proper time rate at transmitter
    std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings_;

    //! Settings for proper time rate at receiver
    std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings_;

    bool normalizeWithSpeedOfLight_;
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
                                  std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings, nullptr ),
        uplinkOneWayDopplerSettings_( uplinkOneWayDopplerSettings ),
        downlinkOneWayDopplerSettings_( downlinkOneWayDopplerSettings )
    {
        if( uplinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ != downlinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ )
        {
            throw std::runtime_error( "Error when making two-way Doppler observable; up and downlink normalization not consistent." );
        }
        else
        {
            normalizeWithSpeedOfLight_ = uplinkOneWayDopplerSettings->normalizeWithSpeedOfLight_;
        }
    }

    TwoWayDopplerObservationSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings( two_way_doppler, linkEnds,
                                  lightTimeCorrectionsList, biasSettings ), normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight )
    {
        uplinkOneWayDopplerSettings_ = std::make_shared< OneWayDopplerObservationSettings >(
                    getUplinkFromTwoWayLinkEnds( linkEnds ), lightTimeCorrectionsList,
                    nullptr, nullptr, nullptr, lightTimeConvergenceCriteria, normalizeWithSpeedOfLight );
        downlinkOneWayDopplerSettings_ = std::make_shared< OneWayDopplerObservationSettings >(
                    getDownlinkFromTwoWayLinkEnds( linkEnds ), lightTimeCorrectionsList,
                    nullptr, nullptr, nullptr, lightTimeConvergenceCriteria, normalizeWithSpeedOfLight );
    }


    //! Destructor
    ~TwoWayDopplerObservationSettings( ){ }

    //! Settings for the one-way Doppler model of the uplink
    std::shared_ptr< OneWayDopplerObservationSettings > uplinkOneWayDopplerSettings_;

    //! Settings for the one-way Doppler model of the downlink
    std::shared_ptr< OneWayDopplerObservationSettings > downlinkOneWayDopplerSettings_;

    bool normalizeWithSpeedOfLight_;
};




////! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler) observable
//class OneWayDifferencedRangeRateObservationSettings: public ObservationModelSettings
//{
//public:

//    //! Constructor
//    /*!
//     * Constructor
//     * \param integrationTimeFunction Function that returns the integration time of observable as a function of time
//     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
//     * (nullptr if none)
//     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
//     */
//    OneWayDifferencedRangeRateObservationSettings(
//            const LinkDefinition& linkEnds,
//            const std::function< double( const double ) > integrationTimeFunction,
//            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
//            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
//        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrections, biasSettings ),
//        integrationTimeFunction_( integrationTimeFunction ){ }

//    //! Constructor
//    /*!
//     * Constructor
//     * \param integrationTimeFunction Function that returns the integration time of observable as a function of time
//     * \param lightTimeCorrectionsList List of ettings for a single light-time correction that is to be used for the observation model
//     * (empty if none)
//     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
//     */
//    OneWayDifferencedRangeRateObservationSettings(
//            const LinkDefinition& linkEnds,
//            const std::function< double( const double ) > integrationTimeFunction,
//            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
//            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
//            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
//        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings ),
//        integrationTimeFunction_( integrationTimeFunction ){ }

//    //! Destructor
//    ~OneWayDifferencedRangeRateObservationSettings( ){ }

//    //! Function that returns the integration time of observable as a function of time
//    const std::function< double( const double ) > integrationTimeFunction_;

//};


//! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler) observable
class NWayRangeObservationSettings: public ObservationModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param oneWayRangeObsevationSettings List of settings for one-way observables that make up n-way link (each must be for
     * one_way_range_
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_range, mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObsevationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings, nullptr ),
        oneWayRangeObsevationSettings_( oneWayRangeObsevationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ){ }

    //! Constructor
    /*!
     * Constructor for same light-time corrections per link
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used for the observation model
     * (nullptr if none)
     * \param biasSettings Settings for the observation bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_range, linkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObsevationSettings_.push_back(
                        std::make_shared< ObservationModelSettings >(
                            one_way_range, getSingleLegLinkEnds( linkEnds.linkEnds_, i ), lightTimeCorrectionsList,
                            nullptr, lightTimeConvergenceCriteria ) );
        }
    }

    //! Destructor
    ~NWayRangeObservationSettings( ){ }

    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;

};

class NWayDifferencedRangeObservationSettings: public ObservationModelSettings
{
public:

    NWayDifferencedRangeObservationSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObsevationSettings_.push_back(
                        std::make_shared< ObservationModelSettings >(
                            one_way_range, getSingleLegLinkEnds( linkEnds.linkEnds_, i ), lightTimeCorrectionsList, nullptr,
                            lightTimeConvergenceCriteria ) );
        }
    }

    NWayDifferencedRangeObservationSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_differenced_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObsevationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ),
        oneWayRangeObsevationSettings_( oneWayRangeObsevationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ) { }

    std::shared_ptr< ObservationModelSettings > getUndifferencedObservationSettings( )
    {
        return std::make_shared< NWayRangeObservationSettings >(
                oneWayRangeObsevationSettings_, nullptr, multiLegLightTimeConvergenceCriteria_ );
    }

private:

    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;

};

class DsnNWayAveragedDopplerObservationSettings: public ObservationModelSettings
{
public:

    DsnNWayAveragedDopplerObservationSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool subtractDopplerSignature = true ):
        ObservationModelSettings( dsn_n_way_averaged_doppler, linkEnds, lightTimeCorrectionsList, biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
        subtractDopplerSignature_( subtractDopplerSignature )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObsevationSettings_.push_back(
                        std::make_shared< ObservationModelSettings >(
                            one_way_range, getSingleLegLinkEnds( linkEnds.linkEnds_, i ), lightTimeCorrectionsList, nullptr,
                            lightTimeConvergenceCriteria ) );
        }
    }

    DsnNWayAveragedDopplerObservationSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool subtractDopplerSignature = true ):
        ObservationModelSettings( n_way_differenced_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObsevationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ),
        oneWayRangeObsevationSettings_( oneWayRangeObsevationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
        subtractDopplerSignature_( subtractDopplerSignature ){ }

    std::shared_ptr< ObservationModelSettings > getNWayRangeObservationSettings( )
    {
        return std::make_shared< NWayRangeObservationSettings >(
                oneWayRangeObsevationSettings_, nullptr, multiLegLightTimeConvergenceCriteria_ );
    }

    bool getSubtractDopplerSignature( )
    {
        return subtractDopplerSignature_;
    }


private:
    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;

    bool subtractDopplerSignature_;

};

inline std::shared_ptr< ObservationModelSettings > oneWayRangeSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
                one_way_range, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > angularPositionSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
                angular_position, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > relativeAngularPositionSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
                relative_angular_position, linkEnds, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > positionObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                position_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > relativePositionObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
            relative_position_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > velocityObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                velocity_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > eulerAngle313ObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr)
{
    return std::make_shared< ObservationModelSettings >(
                euler_angle_313_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > oneWayOpenLoopDoppler(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool normalizeWithSpeedOfLight = false )
{
    return std::make_shared< OneWayDopplerObservationSettings >(
                linkEnds, lightTimeCorrectionsList, transmitterProperTimeRateSettings, receiverProperTimeRateSettings,
                biasSettings, lightTimeConvergenceCriteria, normalizeWithSpeedOfLight );
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
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList = std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool normalizeWithSpeedOfLight = false )
{
    return std::make_shared< TwoWayDopplerObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, normalizeWithSpeedOfLight );
}

inline std::shared_ptr< ObservationModelSettings > dopplerMeasuredFrequencyObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList = std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ))
{
    return std::make_shared< ObservationModelSettings >(
                doppler_measured_frequency,
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria);
}

inline std::shared_ptr< ObservationModelSettings > oneWayClosedLoopDoppler(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
                one_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > twoWayDifferencedRangeObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    if( linkEnds.linkEnds_.size( ) != 3 )
    {
        throw std::runtime_error( "Error when creating two-way averaged Doppler settings, number of input link ends (" +
                                  std::to_string( linkEnds.linkEnds_.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayDifferencedRangeObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > twoWayDifferencedRangeObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    if( oneWayRangeObsevationSettings.size( ) != 2 )
    {
        throw std::runtime_error( "Error when creating two-way averaged Doppler settings, number of input one-way links (" +
                                  std::to_string( oneWayRangeObsevationSettings.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayDifferencedRangeObservationSettings >(
                oneWayRangeObsevationSettings, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > nWayDifferencedRangeObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< NWayDifferencedRangeObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > nWayDifferencedRangeObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< NWayDifferencedRangeObservationSettings >(
                oneWayRangeObsevationSettings, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > dsnNWayAveragedDopplerObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool subtractDopplerSignature = true )
{
    return std::make_shared< DsnNWayAveragedDopplerObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, subtractDopplerSignature );
}

inline std::shared_ptr< ObservationModelSettings > dsnNWayAveragedDopplerObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool subtractDopplerSignature = true )
{
    return std::make_shared< DsnNWayAveragedDopplerObservationSettings >(
                oneWayRangeObsevationSettings, biasSettings, lightTimeConvergenceCriteria, subtractDopplerSignature );
}


inline std::shared_ptr< ObservationModelSettings > twoWayRange(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    if( oneWayRangeObsevationSettings.size( ) != 2 )
    {
        throw std::runtime_error( "Error when creating two-way range settings, number of input one-way links (" +
                                  std::to_string( oneWayRangeObsevationSettings.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayRangeObservationSettings >(
                oneWayRangeObsevationSettings, biasSettings );
}


inline std::shared_ptr< ObservationModelSettings > twoWayRangeSimple(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    if( linkEnds.linkEnds_.size( ) != 3 )
    {
        throw std::runtime_error( "Error when creating two-way range settings, number of input link ends (" +
                                  std::to_string( linkEnds.linkEnds_.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayRangeObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}


inline std::shared_ptr< ObservationModelSettings > nWayRange(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObsevationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< NWayRangeObservationSettings >(
                oneWayRangeObsevationSettings, biasSettings );
}


inline std::shared_ptr< ObservationModelSettings > nWayRangeSimple(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
        = std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< NWayRangeObservationSettings >(
                linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria(
        const bool iterateCorrections = false,
        const int maximumNumberOfIterations = 50,
        const double fractionOfLightTimeTolerance = TUDAT_NAN,
        const LightTimeFailureHandling failureHandling = accept_without_warning )
{
    return std::make_shared< LightTimeConvergenceCriteria >(
            iterateCorrections, maximumNumberOfIterations, fractionOfLightTimeTolerance, failureHandling );
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
        const LinkDefinition& linkEndsDefinitions,
        const simulation_setup::SystemOfBodies &bodies,
        const LinkEndType linkEndForCalculator )
{
    LinkEnds linkEnds = linkEndsDefinitions.linkEnds_;
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
                        LinkEndId( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, "" );
                if( ( linkEnds.at( receiver ) != referencePointId ) && ( linkEnds.at( transmitter ) != referencePointId ) )
                {
                    properTimeRateInterface = std::make_shared<
                            DirectFirstOrderDopplerProperTimeRateInterface >(
                                linkEndForCalculator, gravitationalParameterFunction,
                                directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, unidentified_link_end,
                                simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                                    LinkEndId( directFirstOrderDopplerProperTimeRateSettings->centralBodyName_, ""), bodies ) );
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
        const LinkDefinition linkEnds,
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
    case constant_time_drift_bias:
    {
        // Check input consistency
        std::shared_ptr< ConstantTimeDriftBiasSettings > constantTimeDriftBiasSettings = std::dynamic_pointer_cast<
                ConstantTimeDriftBiasSettings >( biasSettings );
        if( constantTimeDriftBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making constant time drift observation bias, settings are inconsistent" );
        }

        // Check if size of bias is consistent with requested observable size
        if( constantTimeDriftBiasSettings->timeDriftBias_.rows( ) != ObservationSize )
        {
            throw std::runtime_error( "Error when making constant observation time drift bias, bias size is inconsistent" );
        }
        observationBias = std::make_shared< ConstantTimeDriftBias< ObservationSize > >(
                    constantTimeDriftBiasSettings->timeDriftBias_, observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                    observableType, constantTimeDriftBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ), constantTimeDriftBiasSettings->referenceEpoch_ );
        break;
    }
    case arc_wise_time_drift_bias:
    {
        // Check input consistency
        std::shared_ptr< ArcWiseTimeDriftBiasSettings > arcwiseBiasSettings = std::dynamic_pointer_cast<
                ArcWiseTimeDriftBiasSettings >( biasSettings );
        if( arcwiseBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making arc-wise time drift bias, settings are inconsistent" );
        }

        std::vector< Eigen::Matrix< double, ObservationSize, 1 > > observationBiases;
        for( unsigned int i = 0; i < arcwiseBiasSettings->timeDriftBiases_.size( ); i++ )
        {
            // Check if size of bias is consistent with requested observable size
            if( arcwiseBiasSettings->timeDriftBiases_.at( i ).rows( ) != ObservationSize )
            {
                throw std::runtime_error( "Error when making arc-wise time drift bias, bias size is inconsistent" );
            }
            else
            {
                observationBiases.push_back( arcwiseBiasSettings->timeDriftBiases_.at( i ) );
            }
        }
        observationBias = std::make_shared< ArcWiseTimeDriftBias< ObservationSize > >(
                arcwiseBiasSettings->arcStartTimes_, observationBiases,
                observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                        observableType, arcwiseBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ), arcwiseBiasSettings->referenceEpochs_ );
        break;
    }
    case constant_time_bias:
    {
        // Check input consistency
        std::shared_ptr< ConstantTimeBiasSettings > constantTimeBiasSettings = std::dynamic_pointer_cast< ConstantTimeBiasSettings >( biasSettings );
        if( constantTimeBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making constant time observation bias, settings are inconsistent" );
        }

        observationBias = std::make_shared< ConstantTimeBias< ObservationSize > >(
                constantTimeBiasSettings->timeBias_, observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                        observableType, constantTimeBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ) );
        break;
    }
    case arc_wise_time_bias:
    {
        // Check input consistency
        std::shared_ptr< ArcWiseTimeBiasSettings > arcwiseBiasSettings = std::dynamic_pointer_cast<
                ArcWiseTimeBiasSettings >( biasSettings );
        if( arcwiseBiasSettings == nullptr )
        {
            throw std::runtime_error( "Error when making arc-wise time bias, settings are inconsistent" );
        }

        observationBias = std::make_shared< ArcWiseTimeBias< ObservationSize > >(
                arcwiseBiasSettings->arcStartTimes_, arcwiseBiasSettings->timeBiases_,
                observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                        observableType, arcwiseBiasSettings->linkEndForTime_, linkEnds.size( ) ).at( 0 ) );
        break;
    }
    case two_way_range_time_scale_bias:
    {
        if( observableType != n_way_range )
        {
            throw std::runtime_error( "Error when making two-way range time scale bias, bias is to be applied to wrong observable: " +
                getObservableName( observableType, linkEnds.size( ) ) );
        }
        else if( linkEnds.size( ) != 3 )
        {
            throw std::runtime_error( "Error when making two-way range time scale bias, bias is to be applied to observable with wrong number of link ends: " +
                                      std::to_string( linkEnds.size( ) ) );
        }
        else
        {
            if( linkEnds.at( transmitter ).bodyName_ != "Earth" )
            {
                std::cerr<<" Warning when making two-way range time scale bias, bias is to be applied observable with transmitter not on Earth: "<<linkEnds.at( transmitter ).bodyName_ <<std::endl;
            }

            if( linkEnds.at( receiver ).bodyName_ != "Earth" )
            {
                std::cerr<<" Warning when making two-way range time scale bias, bias is to be applied observable with receiver not on Earth: "<<linkEnds.at( receiver ).bodyName_ <<std::endl;
            }
        }

        observationBias = std::make_shared< TwoWayTimeScaleRangeBias< ObservationSize > >(
            earth_orientation::createDefaultTimeConverter( ),
            bodies.at( linkEnds.at( transmitter ).bodyName_ )->getGroundStation( linkEnds.at( transmitter ).stationName_ )->getNominalStationState( ),
            bodies.at( linkEnds.at( transmitter ).bodyName_ )->getGroundStation( linkEnds.at( receiver ).stationName_ )->getNominalStationState( ) );

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
    case clock_induced_bias:
    {
        if( observableType != one_way_range && observableType != n_way_range )
        {
            throw std::runtime_error( "Error, clock-induced observation bias currently only supported for one- and n-way range" );
        }
        else
        {
            // Check input consistency
            std::shared_ptr< TiminigSystemBiasSettings > clockInducedBiasSettings = std::dynamic_pointer_cast<
                    TiminigSystemBiasSettings >( biasSettings );

            std::vector< int > relevantLinkEndIndices = getLinkEndIndicesForLinkEndIdAtObservable(
                    observableType, linkEnds.linkEnds_, LinkEndId( clockInducedBiasSettings->bodyName_, clockInducedBiasSettings->stationName_ ) );

            if( clockInducedBiasSettings == nullptr )
            {
                throw std::runtime_error( "Error when making clock-induced bias, settings are inconsistent" );
            }

            std::shared_ptr< system_models::TimingSystem > timingSystem;
            if( bodies.count( clockInducedBiasSettings->bodyName_ ) == 0 )
            {
                throw std::runtime_error( "Error when getting timing system of body " + clockInducedBiasSettings->bodyName_ +
                    " for observation bias, no such body exists." );
            }
            else
            {
                if( clockInducedBiasSettings->stationName_ == "" )
                {
                    if( bodies.at( clockInducedBiasSettings->bodyName_ )->getVehicleSystems( ) == nullptr )
                    {
                        throw std::runtime_error( "Error when getting timing system of body " + clockInducedBiasSettings->bodyName_ +
                                                  " for observation bias, body has no vehicle systems." );
                    }
                    else
                    {
                        timingSystem = bodies.at( clockInducedBiasSettings->bodyName_ )->getVehicleSystems( )->getTimingSystem( );
                        if( timingSystem == nullptr )
                        {
                            throw std::runtime_error( "Error when getting timing system of body " + clockInducedBiasSettings->bodyName_ +
                                                      " for observation bias, body has no timing system." );
                        }
                    }
                }
                else
                {
                    if( bodies.at( clockInducedBiasSettings->bodyName_ )->getGroundStationMap( ).count( clockInducedBiasSettings->stationName_ ) == 0 )
                    {
                        throw std::runtime_error( "Error when getting timing system of body " + clockInducedBiasSettings->bodyName_ +
                                                  " and station " + clockInducedBiasSettings->stationName_ +
                                                  " no such station exists" );
                    }
                    else
                    {
                        timingSystem = bodies.at( clockInducedBiasSettings->bodyName_ )->getGroundStation( clockInducedBiasSettings->stationName_ )->getTimingSystem( );
                        if( timingSystem == nullptr )
                        {
                            throw std::runtime_error( "Error when getting timing system of body " + clockInducedBiasSettings->bodyName_ +
                                                      " and station " + clockInducedBiasSettings->stationName_ +
                                                      " station has no timing system" );
                        }
                    }
                }

                // Create combined bias object
                observationBias = std::make_shared< ClockInducedRangeBias< ObservationSize > >(
                        timingSystem, relevantLinkEndIndices,
                        LinkEndId( clockInducedBiasSettings->bodyName_, clockInducedBiasSettings->stationName_ ) );
            }
        }
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

std::function< double ( observation_models::FrequencyBands, observation_models::FrequencyBands ) > getTurnaroundFunction(
    const simulation_setup::SystemOfBodies &bodies,
    const LinkEnds& linkEnds );

std::shared_ptr< ground_stations::StationFrequencyInterpolator > getTransmittingFrequencyInterpolator(
    const simulation_setup::SystemOfBodies &bodies,
    const LinkEnds& linkEnds );

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
            const simulation_setup::SystemOfBodies &bodies,
            ObservableType topLevelObservableType = undefined_observation_model );
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
            const simulation_setup::SystemOfBodies &bodies,
            ObservableType topLevelObservableType = undefined_observation_model )
    {
        using namespace observation_models;


        std::shared_ptr< observation_models::ObservationModel<
                1, ObservationScalarType, TimeType > > observationModel;
        LinkEnds linkEnds = observationSettings->linkEnds_.linkEnds_;


        if ( topLevelObservableType == undefined_observation_model )
        {
            topLevelObservableType = observationSettings->observableType_;
        }


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
                            linkEnds, transmitter, receiver,
                            bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ),
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
                                linkEnds, transmitter, receiver,
                                bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                                observationSettings->lightTimeConvergenceCriteria_ ),
                            observationBias,
                            std::function< ObservationScalarType( const TimeType ) >( ),
                            std::function< ObservationScalarType( const TimeType ) >( ),
                            false );
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

                if( oneWayDopplerSettings->receiverProperTimeRateSettings_ != nullptr )
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
                                linkEnds, transmitter, receiver,
                                bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                                observationSettings->lightTimeConvergenceCriteria_ ),
                            transmitterProperTimeRate,
                            receiverProperTimeRate,
                            observationBias,
                            oneWayDopplerSettings->normalizeWithSpeedOfLight_ );
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
            LinkDefinition uplinkLinkEnds;
            uplinkLinkEnds[ transmitter ] = linkEnds.at( transmitter );
            uplinkLinkEnds[ receiver ] = linkEnds.at( reflector1 );

            LinkDefinition downlinkLinkEnds;
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
                                        one_way_doppler, uplinkLinkEnds, observationSettings->lightTimeCorrectionsList_ ), bodies, topLevelObservableType ) ),
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    std::make_shared< ObservationModelSettings >(
                                        one_way_doppler, downlinkLinkEnds, observationSettings->lightTimeCorrectionsList_ ), bodies, topLevelObservableType ) ),
                            observationBias );
            }
            else
            {
                observationModel = std::make_shared< TwoWayDopplerObservationModel<
                        ObservationScalarType, TimeType > >(
                            linkEnds,
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    twoWayDopplerSettings->uplinkOneWayDopplerSettings_, bodies, topLevelObservableType ) ),
                            std::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                    twoWayDopplerSettings->downlinkOneWayDopplerSettings_, bodies, topLevelObservableType ) ),
                            observationBias, twoWayDopplerSettings->normalizeWithSpeedOfLight_ );
            }

            break;
        }

        case one_way_differenced_range:
        {
//            std::shared_ptr< OneWayDifferencedRangeRateObservationSettings > rangeRateObservationSettings =
//                    std::dynamic_pointer_cast< OneWayDifferencedRangeRateObservationSettings >( observationSettings );
//            if( rangeRateObservationSettings == nullptr )
//            {
//                throw std::runtime_error( "Error when making differenced one-way range rate, input type is inconsistent" );
//            }
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
                            linkEnds, transmitter, receiver,
                            bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ),
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds, transmitter, receiver,
                            bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ),
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

            // Create observation bias object
            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            std::shared_ptr< NWayRangeObservationSettings > nWayRangeObservationSettings =
                    std::dynamic_pointer_cast< NWayRangeObservationSettings >( observationSettings );
            if( nWayRangeObservationSettings != nullptr && nWayRangeObservationSettings->oneWayRangeObsevationSettings_.size( )
                    != linkEnds.size( ) - 1 )
            {
                throw std::runtime_error( "Error when making n-way range, input data is inconsistent" );
            }

            // Create vector of convergence criteria and light time corrections
            std::vector< std::shared_ptr< LightTimeConvergenceCriteria > > singleLegsLightTimeConvergenceCriteriaList;
            std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria;
            std::vector< std::vector< std::shared_ptr< LightTimeCorrectionSettings > > > lightTimeCorrectionsList;
            if ( nWayRangeObservationSettings != nullptr )
            {
                multiLegLightTimeConvergenceCriteria = nWayRangeObservationSettings->multiLegLightTimeConvergenceCriteria_;
                for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
                {
                    if( nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->observableType_ != one_way_range )
                    {
                        throw std::runtime_error( "Error in n-way observable creation, constituent link is not of type 1-way." );
                    }
                    lightTimeCorrectionsList.push_back(
                            nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->lightTimeCorrectionsList_ );
                    singleLegsLightTimeConvergenceCriteriaList.push_back(
                            nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->lightTimeConvergenceCriteria_ );
                    if( nWayRangeObservationSettings->oneWayRangeObsevationSettings_.at( i )->biasSettings_ != nullptr )
                    {
                        throw std::runtime_error( "Error when creating n-way range observable, cannot process single-leg bias, but setting is not null" );
                    }
                }
            }
            else
            {
                multiLegLightTimeConvergenceCriteria = observationSettings->lightTimeConvergenceCriteria_;
                for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
                {
                    lightTimeCorrectionsList.push_back(
                            observationSettings->lightTimeCorrectionsList_ );
                    singleLegsLightTimeConvergenceCriteriaList.push_back(
                            observationSettings->lightTimeConvergenceCriteria_ );
                }
            }

            // Create multi-leg light time calculator
            std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >
                    multiLegLightTimeCalculator = createMultiLegLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds, bodies, topLevelObservableType, lightTimeCorrectionsList,
                            singleLegsLightTimeConvergenceCriteriaList,
                            multiLegLightTimeConvergenceCriteria );

            // Create observation model
            observationModel = std::make_shared< NWayRangeObservationModel< ObservationScalarType, TimeType > >(
                        linkEnds, multiLegLightTimeCalculator, observationBias );
            break;
        }
        case n_way_differenced_range:
        {
            std::shared_ptr< NWayDifferencedRangeObservationSettings > nWayDifferencedRangeObservationSettings =
                    std::dynamic_pointer_cast< NWayDifferencedRangeObservationSettings >( observationSettings );
            if( nWayDifferencedRangeObservationSettings == nullptr )
            {
                throw std::runtime_error( "Error when making n-way differenced range observation model, input type inconsistent" );
            }
            std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel;
            std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel;
            try
            {
                std::shared_ptr< ObservationModelSettings > undifferencedObservationSettings =
                        nWayDifferencedRangeObservationSettings->getUndifferencedObservationSettings( );

                arcStartObservationModel =
                        std::dynamic_pointer_cast< NWayRangeObservationModel< ObservationScalarType, TimeType > >(
                            ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                undifferencedObservationSettings, bodies, topLevelObservableType ) );
                arcEndObservationModel =
                        std::dynamic_pointer_cast< NWayRangeObservationModel< ObservationScalarType, TimeType > >(
                            ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                undifferencedObservationSettings, bodies, topLevelObservableType ) );
            }
            catch( const std::exception& caughtException )
            {
                std::string exceptionText = std::string( caughtException.what( ) );
                throw std::runtime_error( "Error when creating n-way differenced range observation model, error: " + exceptionText );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            observationModel = std::make_shared< NWayDifferencedRangeObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds, arcStartObservationModel, arcEndObservationModel, observationBias );
            break;
        }
        case dsn_n_way_averaged_doppler:
        {
            std::shared_ptr< DsnNWayAveragedDopplerObservationSettings > dsnNWayAveragedDopplerObservationSettings =
                    std::dynamic_pointer_cast< DsnNWayAveragedDopplerObservationSettings >( observationSettings );
            if( dsnNWayAveragedDopplerObservationSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating DSN N-way averaged Doppler observation model, input type "
                                          "inconsistent." );
            }

            std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel;
            std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel;
            try
            {
                std::shared_ptr< ObservationModelSettings > nWayRangeObservationSettings =
                        dsnNWayAveragedDopplerObservationSettings->getNWayRangeObservationSettings( );

                arcStartObservationModel =
                        std::dynamic_pointer_cast< NWayRangeObservationModel< ObservationScalarType, TimeType > >(
                            ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                nWayRangeObservationSettings, bodies, topLevelObservableType ) );
                arcEndObservationModel =
                        std::dynamic_pointer_cast< NWayRangeObservationModel< ObservationScalarType, TimeType > >(
                            ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                                nWayRangeObservationSettings, bodies, topLevelObservableType ) );
            }
            catch( const std::exception& caughtException )
            {
                std::string exceptionText = std::string( caughtException.what( ) );
                throw std::runtime_error( "Error when creating DSN N-way averaged Doppler observation model, error: " +
                exceptionText );
            }

            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias = createObservationBiasCalculator(
                        linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            std::function< double ( observation_models::FrequencyBands, observation_models::FrequencyBands ) > turnaroundRatioFunction =
                getTurnaroundFunction( bodies, linkEnds );

            std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyInterpolator =
                getTransmittingFrequencyInterpolator( bodies, linkEnds );


            std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates;
            for( auto it : linkEnds )
            {
                if( bodies.at( linkEnds.at( it.first ).bodyName_ )->getGroundStationMap( ).count( linkEnds.at( it.first ).stationName_ ) > 0 )
                {
                    stationStates[ it.first ] =
                        bodies.at( linkEnds.at( it.first ).bodyName_ )->getGroundStation( linkEnds.at( it.first ).stationName_ )->getNominalStationState( );
                }
            }

            std::shared_ptr< ground_stations::GroundStationState > receivingStationState =
                    bodies.at( linkEnds.at( receiver ).bodyName_ )->getGroundStation( linkEnds.at( receiver ).stationName_ )->getNominalStationState( );
            observationModel = std::make_shared<
                    DsnNWayAveragedDopplerObservationModel< ObservationScalarType, TimeType > >(
                        linkEnds, arcStartObservationModel, arcEndObservationModel,
                        bodies.getBody( linkEnds.at( observation_models::transmitter ).bodyName_ )->getGroundStation(
                                linkEnds.at( observation_models::transmitter ).stationName_ )->getTransmittingFrequencyCalculator( ),
                        turnaroundRatioFunction, observationBias, stationStates, dsnNWayAveragedDopplerObservationSettings->getSubtractDopplerSignature( ) );
            break;
        }
        case doppler_measured_frequency:
        {

            std::shared_ptr< TwoWayDopplerObservationModel< ObservationScalarType, TimeType > > twoWayDopplerModel;
            try
            {
                auto uplinkOneWaySettings = std::make_shared< OneWayDopplerObservationSettings >(
                    getUplinkFromTwoWayLinkEnds( linkEnds ), observationSettings->lightTimeCorrectionsList_,
                    std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ), nullptr );
                uplinkOneWaySettings->normalizeWithSpeedOfLight_ = false;
                auto downlinkOneWaySettings = std::make_shared< OneWayDopplerObservationSettings >(
                    getDownlinkFromTwoWayLinkEnds( linkEnds ), observationSettings->lightTimeCorrectionsList_,
                    nullptr, std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );
                downlinkOneWaySettings->normalizeWithSpeedOfLight_ = false;

                auto twoWaySettings = std::make_shared< TwoWayDopplerObservationSettings >( uplinkOneWaySettings, downlinkOneWaySettings );

                twoWayDopplerModel =
                    std::dynamic_pointer_cast< TwoWayDopplerObservationModel< ObservationScalarType, TimeType > >(
                    ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                        twoWaySettings, bodies ) );
                if( twoWayDopplerModel == nullptr )
                {
                    throw std::runtime_error( "Error when getting two-way Doppler model for measured frequency observable, model in null." );
                }
            }
            catch( const std::exception& caughtException )
            {
                std::string exceptionText = std::string( caughtException.what( ) );
                throw std::runtime_error( "Error when creating Doppler Measured Frequency observation model, error: " +
                exceptionText );
            }

            // Determine the bias settings
            std::shared_ptr< ObservationBias< 1 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias = createObservationBiasCalculator(
                        linkEnds, observationSettings->observableType_, observationSettings->biasSettings_, bodies );
            }

            // Determine the turnaround ratio function
            std::function< double ( observation_models::FrequencyBands, observation_models::FrequencyBands ) > turnaroundRatioFunction =
                getTurnaroundFunction( bodies, linkEnds );

            // Check if transmitter has frequency calculator
            std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyInterpolator =
                getTransmittingFrequencyInterpolator( bodies, linkEnds );

            std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates;
            for( auto it : linkEnds )
            {
                if( bodies.at( linkEnds.at( it.first ).bodyName_ )->getGroundStationMap( ).count( linkEnds.at( it.first ).stationName_ ) > 0 )
                {
                    stationStates[ it.first ] =
                        bodies.at( linkEnds.at( it.first ).bodyName_ )->getGroundStation( linkEnds.at( it.first ).stationName_ )->getNominalStationState( );
                }
            }

            observationModel = std::make_shared<
                    DopplerMeasuredFrequencyObservationModel< ObservationScalarType, TimeType > >(
                        linkEnds, twoWayDopplerModel,
                        transmittingFrequencyInterpolator,
                        turnaroundRatioFunction,
                        observationBias,
                        stationStates );

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
            const simulation_setup::SystemOfBodies &bodies,
            ObservableType topLevelObservableType = undefined_observation_model )
    {
        using namespace observation_models;
        std::shared_ptr< observation_models::ObservationModel<
                2, ObservationScalarType, TimeType > > observationModel;
        LinkEnds linkEnds = observationSettings->linkEnds_.linkEnds_;

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
                            linkEnds, transmitter, receiver,
                            bodies, topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ),
                        observationBias );

            break;
        }
        case relative_angular_position:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 3 )
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
            if( linkEnds.count( transmitter2 ) == 0 )
            {
                throw std::runtime_error( "Error when making angular position model, no second transmitter found" );
            }


            std::shared_ptr< ObservationBias< 2 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias =
                        createObservationBiasCalculator< 2 >(
                            linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model
            observationModel = std::make_shared< RelativeAngularPositionObservationModel<
                    ObservationScalarType, TimeType > >(
                        linkEnds, createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds, transmitter, receiver, bodies,
                            topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ),
                        createLightTimeCalculator< ObservationScalarType, TimeType >(
                            linkEnds, transmitter2, receiver, bodies,
                            topLevelObservableType, observationSettings->lightTimeCorrectionsList_,
                            observationSettings->lightTimeConvergenceCriteria_ ), observationBias );

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
    static std::shared_ptr< observation_models::ObservationModel< 3, ObservationScalarType, TimeType > > createObservationModel(
            const std::shared_ptr< ObservationModelSettings > observationSettings,
            const simulation_setup::SystemOfBodies &bodies,
            ObservableType topLevelObservableType = undefined_observation_model )
    {
        using namespace observation_models;
        std::shared_ptr< observation_models::ObservationModel<
                3, ObservationScalarType, TimeType > > observationModel;

        LinkEnds linkEnds = observationSettings->linkEnds_.linkEnds_;

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
            if( linkEnds.at( observed_body ).stationName_ != "" )
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
                                   bodies.at( linkEnds.at( observed_body ).bodyName_ ), std::placeholders::_1 ),
                        observationBias );

            break;
        }
        case relative_position_observable:
        {
            // Check consistency input.
            if( linkEnds.size( ) != 2 )
            {
                std::string errorMessage = "Error when making relative position observable model, " + std::to_string( linkEnds.size( ) ) + " link ends found";
                throw std::runtime_error( errorMessage );
            }

            if( linkEnds.count( observed_body ) == 0 )
            {
                throw std::runtime_error( "Error when making relative position observable model, no observed_body found" );
            }
            if( linkEnds.count( observer ) == 0 )
            {
                throw std::runtime_error( "Error when making relative position observable model, no observer found" );
            }

            if( observationSettings->lightTimeCorrectionsList_.size( ) > 0 )
            {
                throw std::runtime_error( "Error when making relative position observable model, found light time corrections" );
            }
            if( linkEnds.at( observed_body ).stationName_ != "" )
            {
                throw std::runtime_error( "Error, cannot yet create position function for reference point" );
            }

            std::shared_ptr< ObservationBias< 3 > > observationBias;
            if( observationSettings->biasSettings_ != nullptr )
            {
                observationBias = createObservationBiasCalculator< 3 >(
                        linkEnds, observationSettings->observableType_, observationSettings->biasSettings_,bodies );
            }

            // Create observation model
            observationModel = std::make_shared< RelativePositionObservationModel<
                    ObservationScalarType, TimeType > >(
                    linkEnds,
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< ObservationScalarType, TimeType >,
                               bodies.at( linkEnds.at( observed_body ).bodyName_ ), std::placeholders::_1 ),
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< ObservationScalarType, TimeType >,
                               bodies.at( linkEnds.at( observer ).bodyName_ ), std::placeholders::_1 ),
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
            if( linkEnds.at( observed_body ).stationName_ != "" )
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
            if( bodies.at( linkEnds.at( observed_body ).bodyName_ )->getRotationalEphemeris( ) == nullptr )
            {
                throw std::runtime_error( "Error, cannot euler angle observable; no rotation model found" );
            }
            else
            {
                toBodyFixedFrameFunction = std::bind(
                            &ephemerides::RotationalEphemeris::getRotationToTargetFrameTemplated< TimeType >,
                            bodies.at( linkEnds.at( observed_body ).bodyName_ )->getRotationalEphemeris( ),
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
            if( linkEnds.at( observed_body ).stationName_ != "" )
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
                                   bodies.at( linkEnds.at( observed_body ).bodyName_ ), std::placeholders::_1 ),
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
        observationModels[ settingsList.at( i )->linkEnds_.linkEnds_ ] = ObservationModelCreator<
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
std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > >
    createObservationSimulators(
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

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > createSortedObservationSimulators(
    const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
    const simulation_setup::SystemOfBodies& bodies )
{
    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > sortedObservationSimulators;

    std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators =
        createObservationSimulators< ObservationScalarType, TimeType >( observationSettingsList, bodies );
    for( unsigned int i = 0; i < observationSimulators.size( ); i++ )
    {
        sortedObservationSimulators[ observationSimulators.at( i )->getObservableType( ) ] = observationSimulators.at( i );
    }
    return sortedObservationSimulators;

}


//! Typedef for list of light time corrections for a list of link ends
typedef std::map< observation_models::LinkEnds,
std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >
PerLinkEndPerLightTimeSolutionCorrections;

template< typename ObservationScalarType, typename TimeType, int ObservationSize  >
std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > getLightTimeCorrections(
        const std::shared_ptr< observation_models::ObservationModel<
        ObservationSize, ObservationScalarType, TimeType> > observationModel )
{
    // Clear list, for current link ends.
    std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;

    std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > singleObservableCorrectionList;

    // Check type of observable
    switch( observationModel->getObservableType( ) )
    {
    case observation_models::one_way_range:
    {
        std::shared_ptr< observation_models::OneWayRangeObservationModel
                < ObservationScalarType, TimeType> > oneWayRangeModel =
                std::dynamic_pointer_cast< observation_models::OneWayRangeObservationModel
                < ObservationScalarType, TimeType> >
                ( observationModel );
        singleObservableCorrectionList = (
                    oneWayRangeModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        break;
    }
    case observation_models::one_way_doppler:
    {
        std::shared_ptr< observation_models::OneWayDopplerObservationModel
                < ObservationScalarType, TimeType> > oneWayRangeModel =
                std::dynamic_pointer_cast< observation_models::OneWayDopplerObservationModel
                < ObservationScalarType, TimeType> >
                ( observationModel );
        singleObservableCorrectionList = (
                    oneWayRangeModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        break;
    }
    case observation_models::two_way_doppler:
    {
        std::shared_ptr< observation_models::TwoWayDopplerObservationModel
                < ObservationScalarType, TimeType> > twoWaDopplerModel =
                std::dynamic_pointer_cast< observation_models::TwoWayDopplerObservationModel
                < ObservationScalarType, TimeType> >
                ( observationModel );
        currentLightTimeCorrections.push_back(
                    twoWaDopplerModel->getUplinkDopplerCalculator( )->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        currentLightTimeCorrections.push_back(
                    twoWaDopplerModel->getDownlinkDopplerCalculator( )->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        break;
    }
    case observation_models::angular_position:
    {
        std::shared_ptr< observation_models::AngularPositionObservationModel
                < ObservationScalarType, TimeType> > angularPositionModel =
                std::dynamic_pointer_cast< observation_models::AngularPositionObservationModel
                < ObservationScalarType, TimeType> >
                ( observationModel );
        singleObservableCorrectionList = (
                    angularPositionModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        break;
    }
    case observation_models::one_way_differenced_range:
    {
        std::shared_ptr< observation_models::OneWayDifferencedRangeObservationModel
                < ObservationScalarType, TimeType> > oneWayDifferencedRangeObservationModel =
                std::dynamic_pointer_cast< observation_models::OneWayDifferencedRangeObservationModel
                < ObservationScalarType, TimeType> >
                ( observationModel );
        currentLightTimeCorrections.push_back(
                    oneWayDifferencedRangeObservationModel->getArcStartLightTimeCalculator( )->
                    getLightTimeCorrection( ) );
        currentLightTimeCorrections.push_back(
                    oneWayDifferencedRangeObservationModel->getArcEndLightTimeCalculator( )->
                    getLightTimeCorrection( ) );

        break;
    }
    case observation_models::n_way_range:
    {
        std::shared_ptr< observation_models::NWayRangeObservationModel< ObservationScalarType, TimeType > > nWayRangeObservationModel =
                std::dynamic_pointer_cast< observation_models::NWayRangeObservationModel< ObservationScalarType, TimeType > >
                ( observationModel );
        std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > > lightTimeCalculatorList =
                nWayRangeObservationModel->getLightTimeCalculators( );
        for( unsigned int i = 0; i < lightTimeCalculatorList.size( ); i++ )
        {
            currentLightTimeCorrections.push_back( lightTimeCalculatorList.at( i )->getLightTimeCorrection( ) );
        }
        break;
    }
    case observation_models::position_observable:
    {
        break;
    }
    case observation_models::euler_angle_313_observable:
    {
        break;
    }
    case observation_models::velocity_observable:
    {
        break;
    }
    case observation_models::relative_angular_position:
    {
        std::shared_ptr< observation_models::RelativeAngularPositionObservationModel< ObservationScalarType, TimeType> > relativeAngularPositionModel =
                std::dynamic_pointer_cast< observation_models::RelativeAngularPositionObservationModel
                < ObservationScalarType, TimeType> >( observationModel );
        currentLightTimeCorrections.push_back(
                    relativeAngularPositionModel->getLightTimeCalculatorFirstTransmitter( )->getLightTimeCorrection( ) );
        currentLightTimeCorrections.push_back(
                    relativeAngularPositionModel->getLightTimeCalculatorSecondTransmitter( )->getLightTimeCorrection( ) );
        break;
    }
    case observation_models::relative_position_observable:
    {
        break;
    }
    default:
        std::string errorMessage =
                "Error in light time correction list creation, observable type " +
                std::to_string( observationModel->getObservableType( ) ) + " not recognized.";
        throw std::runtime_error( errorMessage );
    }

    if( singleObservableCorrectionList.size( ) > 0 )
    {
        currentLightTimeCorrections.push_back( singleObservableCorrectionList );
    }

    return currentLightTimeCorrections;

}

//! Function to retrieve a list of light-time corrections per link end from a list of observation models.
/*!
 *  Function to retrieve a list of light-time corrections per link end from a list of observation models.
 *  \param observationModels Map of observation models (may not be of mixed type) with LinkEnds of observable as map key
 *  \return Map of light-time corrections, with associated link ends as key.
 */
template< typename ObservationScalarType, typename TimeType, int ObservationSize  >
PerLinkEndPerLightTimeSolutionCorrections getLightTimeCorrectionsList(
        const std::map< observation_models::LinkEnds, std::shared_ptr< observation_models::ObservationModel<
        ObservationSize, ObservationScalarType, TimeType> > > observationModels )
{
    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionsList;

    // Retrieve type of observable
    observation_models::ObservableType observableType = observationModels.begin( )->second->getObservableType( );

    // Iterate over link ends
    for( typename  std::map< observation_models::LinkEnds, std::shared_ptr< observation_models::ObservationModel<
         ObservationSize, ObservationScalarType, TimeType> > >::const_iterator
         observationModelIterator = observationModels.begin( );
         observationModelIterator != observationModels.end( ); observationModelIterator++ )
    {
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections =
                getLightTimeCorrections( observationModelIterator->second );

        // Add light-time corrections for current link ends.
        if( currentLightTimeCorrections.size( ) > 0 )
        {
            lightTimeCorrectionsList[ observationModelIterator->first ] = currentLightTimeCorrections;
        }

    }
    return lightTimeCorrectionsList;
}

template< int ObservationSize >
class UndifferencedObservationModelExtractor
{
public:
    template< typename ObservationScalarType, typename TimeType >
    static std::pair< std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > >,
    std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > >
    extract( const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > differencedObservationModel );
};

template< >
class UndifferencedObservationModelExtractor< 1 >
{
public:
    template< typename ObservationScalarType, typename TimeType >
    static std::pair< std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > >,
    std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > >
    extract( const std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > differencedObservationModel )
    {
        std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > firstObservationModel;
        std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > secondObservationModel;

        // Check type of observable
        switch( differencedObservationModel->getObservableType( ) )
        {
        case observation_models::one_way_differenced_range:
        {
            std::shared_ptr< observation_models::OneWayDifferencedRangeObservationModel
                    < ObservationScalarType, TimeType> > oneWayDifferencedRangeObservationModel =
                    std::dynamic_pointer_cast< observation_models::OneWayDifferencedRangeObservationModel
                    < ObservationScalarType, TimeType > >( differencedObservationModel );
            firstObservationModel = std::make_shared< observation_models::OneWayRangeObservationModel< ObservationScalarType, TimeType > >(
                        oneWayDifferencedRangeObservationModel->getLinkEnds( ),
                        oneWayDifferencedRangeObservationModel->getArcStartLightTimeCalculator( ) );
            secondObservationModel = std::make_shared< observation_models::OneWayRangeObservationModel< ObservationScalarType, TimeType > >(
                        oneWayDifferencedRangeObservationModel->getLinkEnds( ),
                        oneWayDifferencedRangeObservationModel->getArcEndLightTimeCalculator( ) );
            break;
        }
        case observation_models::n_way_differenced_range:
        {
            std::shared_ptr< observation_models::NWayDifferencedRangeObservationModel
                    < ObservationScalarType, TimeType> > nWayDifferencedRangeObservationModel =
                    std::dynamic_pointer_cast< observation_models::NWayDifferencedRangeObservationModel
                    < ObservationScalarType, TimeType > >( differencedObservationModel );
            firstObservationModel = nWayDifferencedRangeObservationModel->getArcStartObservationModel( );
            secondObservationModel = nWayDifferencedRangeObservationModel->getArcEndObservationModel( );
            break;
        }
        case observation_models::dsn_n_way_averaged_doppler:
        {
            std::shared_ptr< observation_models::DsnNWayAveragedDopplerObservationModel< ObservationScalarType, TimeType > >
                    dsnNWayAveragedDopplerObservationModel =
                    std::dynamic_pointer_cast< observation_models::DsnNWayAveragedDopplerObservationModel< ObservationScalarType, TimeType > >(
                            differencedObservationModel );
            firstObservationModel = dsnNWayAveragedDopplerObservationModel->getArcStartObservationModel( );
            secondObservationModel = dsnNWayAveragedDopplerObservationModel->getArcEndObservationModel( );
            break;
        }
        case observation_models::doppler_measured_frequency:
        {
            throw std::runtime_error( "Error when extract undifferenced observation model. Doppler Measured Frequency model not implemented." );
        }
        default:
            std::string errorMessage =
                    "Error when getting size 1 undifferenced observation models: observable type " +
                    std::to_string( differencedObservationModel->getObservableType( ) ) + " not recognized.";
            throw std::runtime_error( errorMessage );
        }
        return std::make_pair( firstObservationModel, secondObservationModel );
    }
};

template< >
class UndifferencedObservationModelExtractor< 2 >
{
public:
    template< typename ObservationScalarType, typename TimeType >
    static std::pair< std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > >,
    std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > >
    extract( const std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > differencedObservationModel )
    {
        std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > firstObservationModel;
        std::shared_ptr< observation_models::ObservationModel< 2, ObservationScalarType, TimeType > > secondObservationModel;

        // Check type of observable
        switch( differencedObservationModel->getObservableType( ) )
        {
        case observation_models::relative_angular_position:
        {
            std::shared_ptr< observation_models::RelativeAngularPositionObservationModel
                    < ObservationScalarType, TimeType> > relativeAngularPositionModel =
                    std::dynamic_pointer_cast< observation_models::RelativeAngularPositionObservationModel
                    < ObservationScalarType, TimeType > >( differencedObservationModel );
            firstObservationModel = std::make_shared< observation_models::AngularPositionObservationModel< ObservationScalarType, TimeType > >(
                        relativeAngularPositionModel->getFirstLinkEnds( ),
                        relativeAngularPositionModel->getLightTimeCalculatorFirstTransmitter( ) );
            secondObservationModel = std::make_shared< observation_models::AngularPositionObservationModel< ObservationScalarType, TimeType > >(
                        relativeAngularPositionModel->getSecondLinkEnds( ),
                        relativeAngularPositionModel->getLightTimeCalculatorSecondTransmitter() );
            break;
        }
        default:
            std::string errorMessage =
                    "Error when getting size 2 undifferenced observation models: observable type " +
                    std::to_string( differencedObservationModel->getObservableType( ) ) + " not recognized.";
            throw std::runtime_error( errorMessage );
        }
        return std::make_pair( firstObservationModel, secondObservationModel );
    }
};




} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATEOBSERVATIONMODEL_H
