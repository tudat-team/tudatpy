/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOBSERVATIONMODELSETTINGS_H
#define TUDAT_CREATEOBSERVATIONMODELSETTINGS_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationBias.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_models
{

#ifndef TUDAT_OBSERVATION_MODEL_SETTINGS_DEFINED
#define TUDAT_OBSERVATION_MODEL_SETTINGS_DEFINED

class LightTimeCorrectionSettings;
//! Base class to define settings for creation of an observation bias model.
/*!
 *  Base class to define settings for creation of an observation bias model. For each specific bias
 * type, a derived class is to be implemented, in which the specific properties of the bias model
 * are given
 */
class ObservationBiasSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param observationBiasType Type of bias model that is to be created.
     */
    ObservationBiasSettings( const observation_models::ObservationBiasTypes observationBiasType ):
        observationBiasType_( observationBiasType )
    {}

    //! Destructor
    virtual ~ObservationBiasSettings( ) {}

    //! Type of bias model that is to be created.
    observation_models::ObservationBiasTypes observationBiasType_;
};

//! Class for defining settings for the creation of a multiple biases for a single observable
class MultipleObservationBiasSettings : public ObservationBiasSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param biasSettingsList List of settings for bias objects that are to be created.
     */
    MultipleObservationBiasSettings( const std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList ):
        ObservationBiasSettings( multiple_observation_biases ), biasSettingsList_( biasSettingsList )
    {}

    //! Destructor
    ~MultipleObservationBiasSettings( ) {}

    //! List of settings for bias objects that are to be created.
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList_;
};

//! Class for defining settings for the creation of a constant absolute or relative observation bias model
class ConstantObservationBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param observationBias Constant bias that is to be added to the observable. The size of this
     * vector must be equal to the size of the observable to which it is assigned. \param
     * useAbsoluteBias Boolean to denote whether an absolute or relative bias is to be created.
     */
    ConstantObservationBiasSettings( const Eigen::VectorXd& observationBias, const bool useAbsoluteBias ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ? ( constant_absolute_bias ) : ( constant_relative_bias ) ),
        observationBias_( observationBias ), useAbsoluteBias_( useAbsoluteBias )
    {}

    //! Destructor
    ~ConstantObservationBiasSettings( ) {}

    //! Constant bias that is to be added to the observable.
    /*!
     *  Constant bias that is to be added to the observable. The size of this vector must be equal
     * to the size of the observable to which it is assigned.
     */
    Eigen::VectorXd observationBias_;

    //! Boolean to denote whether an absolute or relative bias is to be created.
    bool useAbsoluteBias_;
};

//! Class for defining settings for the creation of an arc-wise constant absolute or relative observation bias model
class ArcWiseConstantObservationBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (observationBiases) are used
     * \param observationBiases List of observation biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     * (and current arc) \param useAbsoluteBias Boolean to denote whether an absolute or relative
     * bias is to be created.
     */
    ArcWiseConstantObservationBiasSettings( const std::vector< double >& arcStartTimes,
                                            const std::vector< Eigen::VectorXd >& observationBiases,
                                            const LinkEndType linkEndForTime,
                                            const bool useAbsoluteBias ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ? ( arc_wise_constant_absolute_bias ) : ( arc_wise_constant_relative_bias ) ),
        arcStartTimes_( arcStartTimes ), observationBiases_( observationBiases ), linkEndForTime_( linkEndForTime ),
        useAbsoluteBias_( useAbsoluteBias )
    {}

    //! Constuctor
    /*!
     * Constuctor
     * \param observationBiases Map of observation biases per arc, with bias as map value, and arc
     * start time as map key \param linkEndForTime Link end at which time is to be evaluated to
     * determine current time (and current arc) \param useAbsoluteBias Boolean to denote whether an
     * absolute or relative bias is to be created.
     */
    ArcWiseConstantObservationBiasSettings( const std::map< double, Eigen::VectorXd >& observationBiases,
                                            const LinkEndType linkEndForTime,
                                            const bool useAbsoluteBias ):
        ObservationBiasSettings( ( useAbsoluteBias == true ) ? ( arc_wise_constant_absolute_bias ) : ( arc_wise_constant_relative_bias ) ),
        arcStartTimes_( utilities::createVectorFromMapKeys( observationBiases ) ),
        observationBiases_( utilities::createVectorFromMapValues( observationBiases ) ), linkEndForTime_( linkEndForTime ),
        useAbsoluteBias_( useAbsoluteBias )
    {}

    //! Destructor
    ~ArcWiseConstantObservationBiasSettings( ) {}

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
class ConstantTimeDriftBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param observationBias Constant time drift bias that is to be considered for the observation
     * time. The size of this vector must be equal to the size of the observable to which it is
     * assigned. \param linkEndForTime Link end at which time is to be evaluated to determine
     * current time \param referenceEpoch Reference epoch at which the effect of the time drift is
     * initialised.
     */
    ConstantTimeDriftBiasSettings( const Eigen::VectorXd& timeDriftBias, const LinkEndType linkEndForTime, const double referenceEpoch ):
        ObservationBiasSettings( constant_time_drift_bias ), timeDriftBias_( timeDriftBias ), linkEndForTime_( linkEndForTime ),
        referenceEpoch_( referenceEpoch )
    {}

    //! Destructor
    ~ConstantTimeDriftBiasSettings( ) {}

    //! Constant time drift bias that is to be considered for the observation time.
    /*!
     *  Constant time drift bias that is to be considered for the observation time. The size of this
     * vector must be equal to the size of the observable to which it is assigned.
     */
    Eigen::VectorXd timeDriftBias_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;

    //! Reference epoch at which the effect of the time drift is supposed to be zero
    double referenceEpoch_;
};

//! Class for defining settings for the creation of an arc-wise time drift bias model
class ArcWiseTimeDriftBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (timeDriftBiases) are used
     * \param timeDriftBiases List of time drift biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     * (and current arc) \param referenceEpochs Reference epochs (per arc) at which the time drifts
     * are initialised.
     */
    ArcWiseTimeDriftBiasSettings( const std::vector< double >& arcStartTimes,
                                  const std::vector< Eigen::VectorXd >& timeDriftBiases,
                                  const LinkEndType linkEndForTime,
                                  const std::vector< double > referenceEpochs ):
        ObservationBiasSettings( arc_wise_time_drift_bias ), arcStartTimes_( arcStartTimes ), timeDriftBiases_( timeDriftBiases ),
        linkEndForTime_( linkEndForTime ), referenceEpochs_( referenceEpochs )
    {}

    //! Constuctor
    /*!
     * Constuctor
     * \param timeDriftBiases Map of observation biases per arc, with bias as map value, and arc
     * start time as map key \param linkEndForTime Link end at which time is to be evaluated to
     * determine current time (and current arc) \param referenceEpochs Reference epochs (per arc) at
     * which the time drifts are initialised.
     */
    ArcWiseTimeDriftBiasSettings( const std::map< double, Eigen::VectorXd >& timeDriftBiases,
                                  const LinkEndType linkEndForTime,
                                  const std::vector< double > referenceEpochs ):
        ObservationBiasSettings( arc_wise_time_drift_bias ), arcStartTimes_( utilities::createVectorFromMapKeys( timeDriftBiases ) ),
        timeDriftBiases_( utilities::createVectorFromMapValues( timeDriftBiases ) ), linkEndForTime_( linkEndForTime ),
        referenceEpochs_( referenceEpochs )
    {}

    //! Destructor
    ~ArcWiseTimeDriftBiasSettings( ) {}

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
class ConstantTimeBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param timeBias Constant time bias that is to be considered for the observation time. The
     * size of this vector must be equal to the size of the observable to which it is assigned.
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     */
    ConstantTimeBiasSettings( const double timeBias, const LinkEndType linkEndForTime ):
        ObservationBiasSettings( constant_time_bias ), timeBias_( timeBias ), linkEndForTime_( linkEndForTime )
    {}

    //! Destructor
    ~ConstantTimeBiasSettings( ) {}

    //! Constant time bias that is to be considered for the observation time.
    /*!
     *  Constant time bias that is to be considered for the observation time.
     */
    double timeBias_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;
};

//! Class for defining settings for the creation of an arc-wise time bias model
class ArcWiseTimeBiasSettings : public ObservationBiasSettings
{
public:
    //! Constuctor
    /*!
     * Constuctor
     * \param arcStartTimes Start times for arcs in which biases (timeBiases) are used
     * \param timeBiases List of time biases per arc
     * \param linkEndForTime Link end at which time is to be evaluated to determine current time
     * (and current arc)
     */
    ArcWiseTimeBiasSettings( const std::vector< double >& arcStartTimes,
                             const std::vector< double >& timeBiases,
                             const LinkEndType linkEndForTime ):
        ObservationBiasSettings( arc_wise_time_bias ), arcStartTimes_( arcStartTimes ), timeBiases_( timeBiases ),
        linkEndForTime_( linkEndForTime )
    {}

    //! Constuctor
    /*!
     * Constuctor
     * \param timeBiases Map of observation biases per arc, with bias as map value, and arc start
     * time as map key \param linkEndForTime Link end at which time is to be evaluated to determine
     * current time (and current arc)
     */
    ArcWiseTimeBiasSettings( const std::map< double, double >& timeBiases, const LinkEndType linkEndForTime ):
        ObservationBiasSettings( arc_wise_time_bias ), arcStartTimes_( utilities::createVectorFromMapKeys( timeBiases ) ),
        timeBiases_( utilities::createVectorFromMapValues( timeBiases ) ), linkEndForTime_( linkEndForTime )
    {}

    //! Destructor
    ~ArcWiseTimeBiasSettings( ) {}

    //! Start times for arcs in which biases (observationBiases) are used
    std::vector< double > arcStartTimes_;

    //! List of time biases per arc
    std::vector< double > timeBiases_;

    //! Link end at which time is to be evaluated to determine current time (and current arc)
    LinkEndType linkEndForTime_;
};

//! Class for defining settings for the creation of a multiple biases for a single observable
class TiminigSystemBiasSettings : public ObservationBiasSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param biasSettingsList List of settings for bias objects that are to be created.
     */
    TiminigSystemBiasSettings( const std::string& bodyName, const std::string& stationName ):
        ObservationBiasSettings( clock_induced_bias ), bodyName_( bodyName ), stationName_( stationName )
    {}

    //! Destructor
    ~TiminigSystemBiasSettings( ) {}

    const std::string bodyName_;

    const std::string stationName_;
};

inline std::shared_ptr< ObservationBiasSettings > clockInducedBias( const std::string& bodyName, const std::string& stationName )
{
    return std::make_shared< TiminigSystemBiasSettings >( bodyName, stationName );
}

inline std::shared_ptr< ObservationBiasSettings > constantAbsoluteBias( const Eigen::VectorXd& observationBias )
{
    return std::make_shared< ConstantObservationBiasSettings >( observationBias, true );
}

inline std::shared_ptr< ObservationBiasSettings > constantRelativeBias( const Eigen::VectorXd& observationBias )
{
    return std::make_shared< ConstantObservationBiasSettings >( observationBias, false );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseAbsoluteBias( const std::vector< double >& arcStartTimes,
                                                                       const std::vector< Eigen::VectorXd >& observationBiases,
                                                                       const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >( arcStartTimes, observationBiases, linkEndForTime, true );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseAbsoluteBias( const std::map< double, Eigen::VectorXd >& observationBiases,
                                                                       const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >( observationBiases, linkEndForTime, true );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseRelativeBias( const std::vector< double >& arcStartTimes,
                                                                       const std::vector< Eigen::VectorXd >& observationBiases,
                                                                       const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >( arcStartTimes, observationBiases, linkEndForTime, false );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseRelativeBias( const std::map< double, Eigen::VectorXd >& observationBiases,
                                                                       const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseConstantObservationBiasSettings >( observationBiases, linkEndForTime, false );
}

inline std::shared_ptr< ObservationBiasSettings > multipleObservationBiasSettings(
        const std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList )
{
    return std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
}

inline std::shared_ptr< ObservationBiasSettings > constantTimeDriftBias( const Eigen::VectorXd& timeDriftBias,
                                                                         const LinkEndType linkEndForTime,
                                                                         const double referenceEpoch )
{
    return std::make_shared< ConstantTimeDriftBiasSettings >( timeDriftBias, linkEndForTime, referenceEpoch );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeDriftBias( const std::vector< Eigen::VectorXd >& timeDriftBiases,
                                                                        const std::vector< double >& arcStartTimes,
                                                                        const LinkEndType linkEndForTime,
                                                                        const std::vector< double >& referenceEpochs )
{
    return std::make_shared< ArcWiseTimeDriftBiasSettings >( arcStartTimes, timeDriftBiases, linkEndForTime, referenceEpochs );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeDriftBias( const std::map< double, Eigen::VectorXd >& timeDriftBiases,
                                                                        const LinkEndType linkEndForTime,
                                                                        const std::vector< double > referenceEpochs )
{
    return std::make_shared< ArcWiseTimeDriftBiasSettings >( timeDriftBiases, linkEndForTime, referenceEpochs );
}

inline std::shared_ptr< ObservationBiasSettings > constantTimeBias( const double timeBias, const LinkEndType linkEndForTime )
{
    return std::make_shared< ConstantTimeBiasSettings >( timeBias, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias( const std::vector< double >& arcStartTimes,
                                                                   const LinkEndType linkEndForTime )
{
    std::vector< double > timeBiases = std::vector< double >( arcStartTimes.size( ), 0.0 );
    return std::make_shared< ArcWiseTimeBiasSettings >( arcStartTimes, timeBiases, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias( const std::vector< double >& timeBiases,
                                                                   const std::vector< double >& arcStartTimes,
                                                                   const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseTimeBiasSettings >( arcStartTimes, timeBiases, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > arcWiseTimeBias( const std::map< double, double >& timeBiases,
                                                                   const LinkEndType linkEndForTime )
{
    return std::make_shared< ArcWiseTimeBiasSettings >( timeBiases, linkEndForTime );
}

inline std::shared_ptr< ObservationBiasSettings > twoWayTimeScaleRangeBias( )
{
    return std::make_shared< ObservationBiasSettings >( two_way_range_time_scale_bias );
}

//! Class used for defining the settings for an observation model that is to be created.
/*!
 * Class used for defining the settings for an observation model that is to be created. This class
 * allows the type, light-time corrections and bias for the observation to be set. For observation
 * models that require additional information (e.g. integration time, retransmission time, etc.), a
 * specific derived class must be implemented.
 */
class ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor (single light-time correction)
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    ObservationModelSettings( const observation_models::ObservableType observableType,
                              const LinkDefinition linkEnds,
                              const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
                              const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                              const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                      std::make_shared< LightTimeConvergenceCriteria >( ),
                              const bool useNormalization = false ):
        observableType_( observableType ), linkEnds_( linkEnds ), biasSettings_( biasSettings ),
        lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ), observableTimeScale_( basic_astrodynamics::tdb_scale ),
        useNormalization_( useNormalization )
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
     * \param lightTimeCorrectionsList List of settings for a single light-time correction that is
     * to be used for the observation model \param biasSettings Settings for the observation bias
     * model that is to be used (default none: nullptr)
     */
    ObservationModelSettings( const observation_models::ObservableType observableType,
                              const LinkDefinition linkEnds,
                              const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                                      std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                              const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                              const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                      std::make_shared< LightTimeConvergenceCriteria >( ),
                              const basic_astrodynamics::TimeScales observableTimeScale = basic_astrodynamics::tdb_scale,
                              const bool useNormalization = false ):
        observableType_( observableType ), linkEnds_( linkEnds ), lightTimeCorrectionsList_( lightTimeCorrectionsList ),
        biasSettings_( biasSettings ), lightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ),
        observableTimeScale_( observableTimeScale ), useNormalization_( useNormalization )
    {}

    //! Destructor
    virtual ~ObservationModelSettings( ) {}

    //! Type of observation model that is to be created
    observation_models::ObservableType observableType_;

    LinkDefinition linkEnds_;

    //! List of settings for a single light-time correction that is to be used for the observation model
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList_;

    //! Settings for the observation bias model that is to be used (default none: nullptr)
    std::shared_ptr< ObservationBiasSettings > biasSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria_;

    basic_astrodynamics::TimeScales observableTimeScale_;

    bool useNormalization_;
};

std::vector< LinkDefinition > getObservationModelListLinkEnds(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelList );

//! Enum defining all possible types of proper time rate computations in one-way Doppler
enum DopplerProperTimeRateType { custom_doppler_proper_time_rate, direct_first_order_doppler_proper_time_rate };

//! Base class to define the settings for proper time rate (at a single link end) in one-way Doppler mode.
class DopplerProperTimeRateSettings
{
public:
    DopplerProperTimeRateSettings( const DopplerProperTimeRateType dopplerProperTimeRateType ):
        dopplerProperTimeRateType_( dopplerProperTimeRateType )
    {}

    virtual ~DopplerProperTimeRateSettings( ) {}

    DopplerProperTimeRateType dopplerProperTimeRateType_;
};

//! Class to define the settings for first-order, single body, proper time rate (at a single link end) in one-way Doppler mode.
class DirectFirstOrderDopplerProperTimeRateSettings : public DopplerProperTimeRateSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param centralBodyName Name of central body, fromw which the mass monopole is retrieved to
     * compute the proper time rate, and w.r.t. which the velocity of the point at which proper time
     * rate is computed is taken
     */
    DirectFirstOrderDopplerProperTimeRateSettings( const std::string& centralBodyName ):
        DopplerProperTimeRateSettings( direct_first_order_doppler_proper_time_rate ), centralBodyNames_( { centralBodyName } )
    {}

    DirectFirstOrderDopplerProperTimeRateSettings( const std::vector< std::string >& centralBodyNames ):
        DopplerProperTimeRateSettings( direct_first_order_doppler_proper_time_rate ), centralBodyNames_( centralBodyNames )
    {}

    //! Destructor.
    ~DirectFirstOrderDopplerProperTimeRateSettings( ) {}

    //! Name of central body
    /*!
     * Name of central body, fromw which the mass monopole is retrieved to compute the proper time
     * rate, and w.r.t. which the velocity of the point at which proper time rate is computed is
     * taken
     */
    std::vector< std::string > centralBodyNames_;
};

//! Class to define the settings for one-way Doppler observable
class OneWayDopplerObservationModelSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param transmitterProperTimeRateSettings Settings
     * for proper time rate at transmitter \param receiverProperTimeRateSettings Settings for proper
     * time rate at receiver \param biasSettings Settings for the observation bias model that is to
     * be used (default none: NUL
     */
    OneWayDopplerObservationModelSettings(
            const LinkDefinition& linkEnds,
            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings( one_way_doppler, linkEnds, lightTimeCorrections, biasSettings, lightTimeConvergenceCriteria ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ), normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight )
    {}

    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrectionsList List of settings for a single light-time correction that is
     * to be used for the observation model (empty if none) \param transmitterProperTimeRateSettings
     * Settings for proper time rate at transmitter \param receiverProperTimeRateSettings Settings
     * for proper time rate at receiver \param biasSettings Settings for the observation bias model
     * that is to be used (default none: NUL
     */
    OneWayDopplerObservationModelSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
            const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings( one_way_doppler, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria ),
        transmitterProperTimeRateSettings_( transmitterProperTimeRateSettings ),
        receiverProperTimeRateSettings_( receiverProperTimeRateSettings ), normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight )
    {}

    //! Destructor
    ~OneWayDopplerObservationModelSettings( ) {}

    //! Settings for proper time rate at transmitter
    std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings_;

    //! Settings for proper time rate at receiver
    std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings_;

    bool normalizeWithSpeedOfLight_;
};

//! Class to define the settings for one-way Doppler observable
class TwoWayDopplerObservationModelSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param uplinkOneWayDopplerSettings Settings for the one-way Doppler model of the uplink
     * \param downlinkOneWayDopplerSettings Settings for the one-way Doppler model of the downlink
     * \param biasSettings Settings for the observation bias model that is to be used (default none:
     * NUL
     */
    TwoWayDopplerObservationModelSettings( const std::shared_ptr< OneWayDopplerObservationModelSettings > uplinkOneWayDopplerSettings,
                                           const std::shared_ptr< OneWayDopplerObservationModelSettings > downlinkOneWayDopplerSettings,
                                           const ObservableType observableType = two_way_doppler,
                                           const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
        ObservationModelSettings( observableType,
                                  mergeUpDownLink( uplinkOneWayDopplerSettings->linkEnds_, downlinkOneWayDopplerSettings->linkEnds_ ),
                                  std::shared_ptr< LightTimeCorrectionSettings >( ),
                                  biasSettings,
                                  std::make_shared< LightTimeConvergenceCriteria >( ) ),
        uplinkOneWayDopplerSettings_( uplinkOneWayDopplerSettings ), downlinkOneWayDopplerSettings_( downlinkOneWayDopplerSettings )
    {
        if( observableType != two_way_doppler && observableType != doppler_measured_frequency )
        {
            throw std::runtime_error( "Error when  creating 2-way Doppler settings, input type is inconsistent" );
        }

        if( observableType == doppler_measured_frequency )
        {
            uplinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ = false;
            downlinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ = false;
        }
        if( uplinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ != downlinkOneWayDopplerSettings->normalizeWithSpeedOfLight_ )
        {
            throw std::runtime_error(
                    "Error when making two-way Doppler observable; up and downlink normalization "
                    "not consistent." );
        }
        else
        {
            normalizeWithSpeedOfLight_ = uplinkOneWayDopplerSettings->normalizeWithSpeedOfLight_;
        }
    }

    TwoWayDopplerObservationModelSettings( const LinkDefinition& linkEnds,
                                           const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                                                   std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                           const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                                           const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                                   std::make_shared< LightTimeConvergenceCriteria >( ),
                                           const bool normalizeWithSpeedOfLight = false ):
        ObservationModelSettings( two_way_doppler, linkEnds, lightTimeCorrectionsList, biasSettings ),
        normalizeWithSpeedOfLight_( normalizeWithSpeedOfLight )
    {
        if( linkEnds.getLinkEnds( ).count( transmitter ) == 0 )
        {
            throw std::runtime_error( "Error when creating 2-way Doppler model, no transmitter found" );
        }
        if( linkEnds.getLinkEnds( ).count( receiver ) == 0 )
        {
            throw std::runtime_error( "Error when creating 2-way Doppler model, no receiver found" );
        }
        if( linkEnds.getLinkEnds( ).count( retransmitter ) == 0 )
        {
            throw std::runtime_error( "Error when creating 2-way Doppler model, no retransmitter found" );
        }
        uplinkOneWayDopplerSettings_ = std::make_shared< OneWayDopplerObservationModelSettings >( getUplinkFromTwoWayLinkEnds( linkEnds ),
                                                                                                  lightTimeCorrectionsList,
                                                                                                  nullptr,
                                                                                                  nullptr,
                                                                                                  nullptr,
                                                                                                  lightTimeConvergenceCriteria,
                                                                                                  normalizeWithSpeedOfLight );
        downlinkOneWayDopplerSettings_ =
                std::make_shared< OneWayDopplerObservationModelSettings >( getDownlinkFromTwoWayLinkEnds( linkEnds ),
                                                                           lightTimeCorrectionsList,
                                                                           nullptr,
                                                                           nullptr,
                                                                           nullptr,
                                                                           lightTimeConvergenceCriteria,
                                                                           normalizeWithSpeedOfLight );
    }

    //! Destructor
    ~TwoWayDopplerObservationModelSettings( ) {}

    //! Settings for the one-way Doppler model of the uplink
    std::shared_ptr< OneWayDopplerObservationModelSettings > uplinkOneWayDopplerSettings_;

    //! Settings for the one-way Doppler model of the downlink
    std::shared_ptr< OneWayDopplerObservationModelSettings > downlinkOneWayDopplerSettings_;

    bool normalizeWithSpeedOfLight_;
};

////! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler)
/// observable
// class OneWayDifferencedRangeRateObservationSettings: public ObservationModelSettings
//{
// public:

//    //! Constructor
//    /*!
//     * Constructor
//     * \param integrationTimeFunction Function that returns the integration time of observable as
//     a function of time
//     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
//     for the observation model
//     * (nullptr if none)
//     * \param biasSettings Settings for the observation bias model that is to be used (default
//     none: nullptr)
//     */
//    OneWayDifferencedRangeRateObservationSettings(
//            const LinkDefinition& linkEnds,
//            const std::function< double( const double ) > integrationTimeFunction,
//            const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections,
//            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
//        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrections,
//        biasSettings ), integrationTimeFunction_( integrationTimeFunction ){ }

//    //! Constructor
//    /*!
//     * Constructor
//     * \param integrationTimeFunction Function that returns the integration time of observable as
//     a function of time
//     * \param lightTimeCorrectionsList List of ettings for a single light-time correction that is
//     to be used for the observation model
//     * (empty if none)
//     * \param biasSettings Settings for the observation bias model that is to be used (default
//     none: nullptr)
//     */
//    OneWayDifferencedRangeRateObservationSettings(
//            const LinkDefinition& linkEnds,
//            const std::function< double( const double ) > integrationTimeFunction,
//            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >
//            lightTimeCorrectionsList = std::vector< std::shared_ptr< LightTimeCorrectionSettings >
//            >( ), const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr ):
//        ObservationModelSettings( one_way_differenced_range, linkEnds, lightTimeCorrectionsList,
//        biasSettings ), integrationTimeFunction_( integrationTimeFunction ){ }

//    //! Destructor
//    ~OneWayDifferencedRangeRateObservationSettings( ){ }

//    //! Function that returns the integration time of observable as a function of time
//    const std::function< double( const double ) > integrationTimeFunction_;

//};

//! Class to define the settings for one-way differenced range-rate (e.g. closed-loop Doppler) observable
class NWayRangeObservationModelSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param oneWayRangeObservationSettings List of settings for one-way observables that make up
     * n-way link (each must be for one_way_range_ \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationModelSettings( const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
                                       const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                                       const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                               std::make_shared< LightTimeConvergenceCriteria >( ),
                                       const basic_astrodynamics::TimeScales observableTimeScale = basic_astrodynamics::tdb_scale ):
        ObservationModelSettings( n_way_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObservationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                  biasSettings,
                                  nullptr,
                                  observableTimeScale ),
        oneWayRangeObservationSettings_( oneWayRangeObservationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {}

    //! Constructor
    /*!
     * Constructor for same light-time corrections per link
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    NWayRangeObservationModelSettings( const LinkDefinition& linkEnds,
                                       const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList,
                                       const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                                       const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                               std::make_shared< LightTimeConvergenceCriteria >( ),
                                       const basic_astrodynamics::TimeScales observableTimeScale = basic_astrodynamics::tdb_scale ):
        ObservationModelSettings( n_way_range,
                                  linkEnds,
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                  biasSettings,
                                  std::make_shared< LightTimeConvergenceCriteria >( ),
                                  observableTimeScale ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObservationSettings_.push_back(
                    std::make_shared< ObservationModelSettings >( one_way_range,
                                                                  getSingleLegLinkEnds( linkEnds.linkEnds_, i ),
                                                                  lightTimeCorrectionsList,
                                                                  nullptr,
                                                                  lightTimeConvergenceCriteria ) );
        }
    }

    //! Destructor
    ~NWayRangeObservationModelSettings( ) {}

    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;
};

class NWayDifferencedRangeObservationModelSettings : public ObservationModelSettings
{
public:
    NWayDifferencedRangeObservationModelSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObservationSettings_.push_back(
                    std::make_shared< ObservationModelSettings >( one_way_range,
                                                                  getSingleLegLinkEnds( linkEnds.linkEnds_, i ),
                                                                  lightTimeCorrectionsList,
                                                                  nullptr,
                                                                  lightTimeConvergenceCriteria ) );
        }
    }

    NWayDifferencedRangeObservationModelSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( n_way_differenced_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObservationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                  biasSettings ),
        oneWayRangeObservationSettings_( oneWayRangeObservationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {}

    std::shared_ptr< ObservationModelSettings > getUndifferencedObservationSettings( )
    {
        return std::make_shared< NWayRangeObservationModelSettings >(
                oneWayRangeObservationSettings_, nullptr, multiLegLightTimeConvergenceCriteria_ );
    }

private:
    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;
};

class DsnNWayAveragedDopplerObservationModelSettings : public ObservationModelSettings
{
public:
    DsnNWayAveragedDopplerObservationModelSettings(
            const LinkDefinition& linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool subtractDopplerSignature = true ):
        ObservationModelSettings( dsn_n_way_averaged_doppler, linkEnds, lightTimeCorrectionsList, biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ), subtractDopplerSignature_( subtractDopplerSignature )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObservationSettings_.push_back(
                    std::make_shared< ObservationModelSettings >( one_way_range,
                                                                  getSingleLegLinkEnds( linkEnds.linkEnds_, i ),
                                                                  lightTimeCorrectionsList,
                                                                  nullptr,
                                                                  lightTimeConvergenceCriteria ) );
        }
    }

    DsnNWayAveragedDopplerObservationModelSettings(
            const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ),
            const bool subtractDopplerSignature = true ):
        ObservationModelSettings( n_way_differenced_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObservationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                  biasSettings ),
        oneWayRangeObservationSettings_( oneWayRangeObservationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria ), subtractDopplerSignature_( subtractDopplerSignature )
    {}

    std::shared_ptr< ObservationModelSettings > getNWayRangeObservationSettings( )
    {
        return std::make_shared< NWayRangeObservationModelSettings >(
                oneWayRangeObservationSettings_, nullptr, multiLegLightTimeConvergenceCriteria_ );
    }

    bool getSubtractDopplerSignature( )
    {
        return subtractDopplerSignature_;
    }

private:
    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings_;

    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;

    bool subtractDopplerSignature_;
};

class DsnNWayRangeObservationModelSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor for same light-time corrections per link
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none)
     * \param biasSettings Settings for the observation bias model that is to be used (default none:
     * nullptr)
     */
    DsnNWayRangeObservationModelSettings( const LinkDefinition& linkEnds,
                                          const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                          const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                                          const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                                  std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( dsn_n_way_range, linkEnds, lightTimeCorrectionsList, biasSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {
        for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
        {
            oneWayRangeObservationSettings_.push_back(
                    std::make_shared< ObservationModelSettings >( one_way_range,
                                                                  getSingleLegLinkEnds( linkEnds.linkEnds_, i ),
                                                                  lightTimeCorrectionsList,
                                                                  nullptr,
                                                                  lightTimeConvergenceCriteria ) );
        }
    }

    //! Constructor
    /*!
     * Constructor for different light-time corrections per link
     * \param oneWayRangeObservationSettings List of settings for the one-way range observation
     * models
     * \param biasSettings Settings for the observation bias model that is to be used (default none:
     * nullptr)
     * \param lightTimeConvergenceCriteria Settings for the light-time convergence criteria
     */
    DsnNWayRangeObservationModelSettings( const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
                                          const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
                                          const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                                                  std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( dsn_n_way_range,
                                  mergeOneWayLinkEnds( getObservationModelListLinkEnds( oneWayRangeObservationSettings ) ),
                                  std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                  biasSettings ),
        oneWayRangeObservationSettings_( oneWayRangeObservationSettings ),
        multiLegLightTimeConvergenceCriteria_( lightTimeConvergenceCriteria )
    {}

    //! Destructor
    ~DsnNWayRangeObservationModelSettings( ) {}

    std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings_;
    std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria_;
};

class DifferencedTimeOfArrivalObservationSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor (single light-time correction)
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    DifferencedTimeOfArrivalObservationSettings(
            const LinkDefinition linkEnds,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections,
            const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( differenced_time_of_arrival, linkEnds, lightTimeCorrections, biasSettings, lightTimeConvergenceCriteria ),
        differencedTimeScale_( differencedTimeScale )
    {}

    basic_astrodynamics::TimeScales differencedTimeScale_;
};

class OneWayDopplerMeasuredFrequencyObservationSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    OneWayDopplerMeasuredFrequencyObservationSettings(
            const LinkDefinition linkEnds,
            const std::shared_ptr< OneWayDopplerObservationModelSettings > dopplerModelSettings,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections,
            const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( one_way_doppler_measured_frequency,
                                  linkEnds,
                                  lightTimeCorrections,
                                  biasSettings,
                                  lightTimeConvergenceCriteria,
                                  differencedTimeScale ),
        differencedTimeScale_( differencedTimeScale ), dopplerModelSettings_( dopplerModelSettings )
    {
        if( differencedTimeScale_ != basic_astrodynamics::tdb_scale &&
            differencedTimeScale_ != basic_astrodynamics::utc_scale )
        {
            throw std::runtime_error(
                    "Error when creating one-way Doppler measured frequency observation settings, "
                    "only TDB and UTC time scales are currently supported" );
        }
    }

    [[nodiscard]] std::shared_ptr< OneWayDopplerObservationModelSettings > getDopplerModelSettings( )
    {
        return dopplerModelSettings_;
    }

private:
    basic_astrodynamics::TimeScales differencedTimeScale_;

    std::shared_ptr< OneWayDopplerObservationModelSettings > dopplerModelSettings_;
};

class DifferencedFrequencyOfArrivalObservationSettings : public ObservationModelSettings
{
public:
    //! Constructor
    /*!
     * Constructor (single light-time correction)
     * \param observableType Type of observation model that is to be created
     * \param lightTimeCorrections Settings for a single light-time correction that is to be used
     * for the observation model (nullptr if none) \param biasSettings Settings for the observation
     * bias model that is to be used (default none: nullptr)
     */
    DifferencedFrequencyOfArrivalObservationSettings(
            const LinkDefinition linkEnds,
            const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > firstDopplerModelSettings,
            const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > secondDopplerModelSettings,
            const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections,
            const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
            const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                    std::make_shared< LightTimeConvergenceCriteria >( ) ):
        ObservationModelSettings( differenced_frequency_of_arrival,
                                  linkEnds,
                                  lightTimeCorrections,
                                  biasSettings,
                                  lightTimeConvergenceCriteria,
                                  differencedTimeScale ),
        differencedTimeScale_( differencedTimeScale ), firstDopplerModelSettings_( firstDopplerModelSettings ),
        secondDopplerModelSettings_( secondDopplerModelSettings )
    {
        if( differencedTimeScale_ != basic_astrodynamics::tdb_scale &&
            differencedTimeScale_ != basic_astrodynamics::utc_scale )
        {
            throw std::runtime_error(
                    "Error when creating differenced frequency of arrival observation settings, "
                    "only TDB and UTC time scales are currently supported" );
        }
    }

    [[nodiscard]] std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > getFirstDopplerModelSettings( ) const
    {
        return firstDopplerModelSettings_;
    }

    [[nodiscard]] std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > getSecondDopplerModelSettings( ) const
    {
        return secondDopplerModelSettings_;
    }

private:
    basic_astrodynamics::TimeScales differencedTimeScale_;

    std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > firstDopplerModelSettings_;
    std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationSettings > secondDopplerModelSettings_;
};

inline std::shared_ptr< ObservationModelSettings > oneWayRangeSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale )
{
    return std::make_shared< ObservationModelSettings >(
            one_way_range, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, timeScale );
}

inline std::shared_ptr< ObservationModelSettings > angularPositionSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool normalizeRightAscension = false )
{
    return std::make_shared< ObservationModelSettings >( angular_position,
                                                         linkEnds,
                                                         lightTimeCorrectionsList,
                                                         biasSettings,
                                                         lightTimeConvergenceCriteria,
                                                         basic_astrodynamics::tdb_scale,
                                                         normalizeRightAscension );
}

inline std::shared_ptr< ObservationModelSettings > relativeAngularPositionSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >( relative_angular_position, linkEnds, lightTimeCorrectionsList, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > positionObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< ObservationModelSettings >( position_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > relativePositionObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< ObservationModelSettings >( relative_position_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > velocityObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< ObservationModelSettings >( velocity_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > eulerAngle313ObservableSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< ObservationModelSettings >( euler_angle_313_observable, linkEnds, nullptr, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > oneWayOpenLoopDoppler(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings = nullptr,
        const std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool normalizeWithSpeedOfLight = false )
{
    return std::make_shared< OneWayDopplerObservationModelSettings >( linkEnds,
                                                                      lightTimeCorrectionsList,
                                                                      transmitterProperTimeRateSettings,
                                                                      receiverProperTimeRateSettings,
                                                                      biasSettings,
                                                                      lightTimeConvergenceCriteria,
                                                                      normalizeWithSpeedOfLight );
}

inline std::shared_ptr< ObservationModelSettings > twoWayOpenLoopDoppler(
        const std::shared_ptr< OneWayDopplerObservationModelSettings > uplinkOneWayDopplerSettings,
        const std::shared_ptr< OneWayDopplerObservationModelSettings > downlinkOneWayDopplerSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    return std::make_shared< TwoWayDopplerObservationModelSettings >(
            uplinkOneWayDopplerSettings, downlinkOneWayDopplerSettings, two_way_doppler, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > twoWayOpenLoopDoppler(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool normalizeWithSpeedOfLight = false )
{
    return std::make_shared< TwoWayDopplerObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, normalizeWithSpeedOfLight );
}

inline std::shared_ptr< ObservationModelSettings > dopplerMeasuredFrequencyObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
            doppler_measured_frequency, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > dopplerMeasuredFrequencyObservationSettings(
        const LinkDefinition& linkEnds,
        const std::shared_ptr< DopplerProperTimeRateSettings > properTimeRateSettings,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    if( properTimeRateSettings != nullptr )
    {
        std::shared_ptr< OneWayDopplerObservationModelSettings > uplinkOneWaySettings =
                std::make_shared< OneWayDopplerObservationModelSettings >( getUplinkFromTwoWayLinkEnds( linkEnds ),
                                                                           lightTimeCorrectionsList,
                                                                           properTimeRateSettings,
                                                                           properTimeRateSettings,
                                                                           nullptr,
                                                                           lightTimeConvergenceCriteria );
        uplinkOneWaySettings->normalizeWithSpeedOfLight_ = false;

        std::shared_ptr< OneWayDopplerObservationModelSettings > downlinkOneWaySettings =
                std::make_shared< OneWayDopplerObservationModelSettings >( getDownlinkFromTwoWayLinkEnds( linkEnds ),
                                                                           lightTimeCorrectionsList,
                                                                           properTimeRateSettings,
                                                                           properTimeRateSettings,
                                                                           nullptr,
                                                                           lightTimeConvergenceCriteria );
        downlinkOneWaySettings->normalizeWithSpeedOfLight_ = false;

        return std::make_shared< TwoWayDopplerObservationModelSettings >(
                uplinkOneWaySettings, downlinkOneWaySettings, doppler_measured_frequency, biasSettings );
    }
    else
    {
        return std::make_shared< ObservationModelSettings >(
                doppler_measured_frequency, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
    }
}

inline std::shared_ptr< ObservationModelSettings > oneWayClosedLoopDoppler(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< ObservationModelSettings >(
            one_way_differenced_range, linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > twoWayDifferencedRangeObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    if( linkEnds.linkEnds_.size( ) != 3 )
    {
        throw std::runtime_error(
                "Error when creating two-way averaged Doppler settings, number of input link ends "
                "(" +
                std::to_string( linkEnds.linkEnds_.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayDifferencedRangeObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > twoWayDifferencedRangeObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr )
{
    if( oneWayRangeObservationSettings.size( ) != 2 )
    {
        throw std::runtime_error(
                "Error when creating two-way averaged Doppler settings, number of input one-way "
                "links (" +
                std::to_string( oneWayRangeObservationSettings.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayDifferencedRangeObservationModelSettings >( oneWayRangeObservationSettings, biasSettings );
}

inline std::shared_ptr< ObservationModelSettings > nWayDifferencedRangeObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< NWayDifferencedRangeObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > nWayDifferencedRangeObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< NWayDifferencedRangeObservationModelSettings >(
            oneWayRangeObservationSettings, biasSettings, lightTimeConvergenceCriteria );
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
    return std::make_shared< DsnNWayAveragedDopplerObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, subtractDopplerSignature );
}

inline std::shared_ptr< ObservationModelSettings > dsnNWayAveragedDopplerObservationSettings(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const bool subtractDopplerSignature = true )
{
    return std::make_shared< DsnNWayAveragedDopplerObservationModelSettings >(
            oneWayRangeObservationSettings, biasSettings, lightTimeConvergenceCriteria, subtractDopplerSignature );
}

inline std::shared_ptr< ObservationModelSettings > dsnNWayRangeObservationSettings(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< DsnNWayRangeObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > twoWayRange(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale )
{
    if( oneWayRangeObservationSettings.size( ) != 2 )
    {
        throw std::runtime_error( "Error when creating two-way range settings, number of input one-way links (" +
                                  std::to_string( oneWayRangeObservationSettings.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayRangeObservationModelSettings >(
            oneWayRangeObservationSettings, biasSettings, std::make_shared< LightTimeConvergenceCriteria >( ), timeScale );
}

inline std::shared_ptr< ObservationModelSettings > twoWayRangeSimple(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale )
{
    if( linkEnds.linkEnds_.size( ) != 3 )
    {
        throw std::runtime_error( "Error when creating two-way range settings, number of input link ends (" +
                                  std::to_string( linkEnds.linkEnds_.size( ) ) + ") is incompatible." );
    }
    return std::make_shared< NWayRangeObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, timeScale );
}

inline std::shared_ptr< ObservationModelSettings > nWayRange(
        const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayRangeObservationSettings,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale )
{
    return std::make_shared< NWayRangeObservationModelSettings >(
            oneWayRangeObservationSettings, biasSettings, std::make_shared< LightTimeConvergenceCriteria >( ), timeScale );
}

inline std::shared_ptr< ObservationModelSettings > nWayRangeSimple(
        const LinkDefinition& linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrectionsList =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ),
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale )
{
    return std::make_shared< NWayRangeObservationModelSettings >(
            linkEnds, lightTimeCorrectionsList, biasSettings, lightTimeConvergenceCriteria, timeScale );
}

inline std::shared_ptr< ObservationModelSettings > differencedTimeOfArrivalObservationSettings(
        const LinkDefinition linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections,
        const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    return std::make_shared< DifferencedTimeOfArrivalObservationSettings >(
            linkEnds, lightTimeCorrections, differencedTimeScale, biasSettings, lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > differencedFrequencyOfArrivalObservationSettings(
        const LinkDefinition linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
    const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
        std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    // Separate link ends for both doppler models
    LinkEnds firstDopplerLinkEnds;
    firstDopplerLinkEnds[ transmitter ] = linkEnds.at( transmitter );
    firstDopplerLinkEnds[ receiver ] = linkEnds.at( receiver );
    LinkEnds secondDopplerLinkEnds;
    secondDopplerLinkEnds[ transmitter ] = linkEnds.at( transmitter );
    secondDopplerLinkEnds[ receiver ] = linkEnds.at( receiver2 );

    // Create one-way Doppler model settings (normalized by c)
    auto firstOneWayDopplerSettings = std::make_shared< OneWayDopplerObservationModelSettings >(
        firstDopplerLinkEnds, lightTimeCorrections, nullptr, nullptr, nullptr, lightTimeConvergenceCriteria, true );
    auto secondOneWayDopplerSettings = std::make_shared< OneWayDopplerObservationModelSettings >(
        secondDopplerLinkEnds, lightTimeCorrections, nullptr, nullptr, nullptr, lightTimeConvergenceCriteria, true );

    // Create OneWayDopplerMeasuredFrequencyObservationSettings for each link
    auto firstOneWayDopplerModelSettings = std::make_shared< OneWayDopplerMeasuredFrequencyObservationSettings >(
        firstDopplerLinkEnds, firstOneWayDopplerSettings, lightTimeCorrections, differencedTimeScale, nullptr, lightTimeConvergenceCriteria );
    auto secondOneWayDopplerModelSettings = std::make_shared< OneWayDopplerMeasuredFrequencyObservationSettings >(
        secondDopplerLinkEnds, secondOneWayDopplerSettings, lightTimeCorrections, differencedTimeScale, nullptr, lightTimeConvergenceCriteria );

    return std::make_shared< DifferencedFrequencyOfArrivalObservationSettings >( linkEnds,
                                         firstOneWayDopplerModelSettings,
                                         secondOneWayDopplerModelSettings,
                                         lightTimeCorrections,
                                         differencedTimeScale,
                                         biasSettings,
                                         lightTimeConvergenceCriteria );
}

inline std::shared_ptr< ObservationModelSettings > oneWayDopplerMeasuredFrequencySettings(
        const LinkDefinition linkEnds,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const basic_astrodynamics::TimeScales differencedTimeScale = basic_astrodynamics::tdb_scale,
        const std::shared_ptr< ObservationBiasSettings > biasSettings = nullptr,
    const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
        std::make_shared< LightTimeConvergenceCriteria >( ) )
{
    // Create one-way Doppler model settings (normalized by c, i.e. returns fractional Doppler)
    auto oneWayDopplerModelSettings =
        std::make_shared< OneWayDopplerObservationModelSettings >( linkEnds,
                                       lightTimeCorrections,
                                       nullptr,  // transmitterProperTimeRateSettings
                                       nullptr,  // receiverProperTimeRateSettings
                                       nullptr,  // biasSettings
                                       lightTimeConvergenceCriteria,
                                       true );  // normalizeWithSpeedOfLight

    return std::make_shared< OneWayDopplerMeasuredFrequencyObservationSettings >(
        linkEnds, oneWayDopplerModelSettings, lightTimeCorrections, differencedTimeScale, biasSettings, lightTimeConvergenceCriteria );
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

#endif  // TUDAT_OBSERVATION_MODEL_SETTINGS_DEFINED

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CREATEOBSERVATIONMODELSETTINGS_H
