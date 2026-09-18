/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONPROCESSINGSETTINGS_H
#define TUDAT_PROPAGATIONPROCESSINGSETTINGS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <unordered_map>

#include <Eigen/Core>

#include "tudat/simulation/propagation_setup/propagationPrintSettings.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/io/serialization/core.h"
#include "tudat/io/serialization/file_io_declarations.h"

namespace tudat
{

namespace propagators
{

//! Base class for defining output and processing settings for propagation.
//! This class is inherited for the separate cases of single, multi and hybrid
//! arc. Each derived class defines whether the propagation results are to be
//! used to reset the environment, and whether the numerical solution is to be
//! deleted after the propagation.
class PropagatorProcessingSettings
{
public:
    PropagatorProcessingSettings( const bool clearNumericalSolutions = false,
                                  const bool setIntegratedResult = false,
                                  const bool updateDependentVariableInterpolator = false,
                                  const bool setIntegratedVariationalResult = true ):
        clearNumericalSolutions_( clearNumericalSolutions ), setIntegratedResult_( setIntegratedResult ),
        createStateProcessors_( setIntegratedResult ), updateDependentVariableInterpolator_( updateDependentVariableInterpolator ),
        setIntegratedVariationalResult_( setIntegratedVariationalResult ),
        interpolatorSettings_( std::make_shared< interpolators::LagrangeInterpolatorSettings >(
                6,
                false,
                interpolators::huntingAlgorithm,
                interpolators::lagrange_cubic_spline_boundary_interpolation,
                interpolators::throw_exception_at_boundary ) )
    {}

    virtual ~PropagatorProcessingSettings( ) {}

    bool getClearNumericalSolutions( )
    {
        return clearNumericalSolutions_;
    }

    bool getSetIntegratedResult( )
    {
        return setIntegratedResult_;
    }

    bool getSetIntegratedVariationalResult( )
    {
        return setIntegratedVariationalResult_;
    }

    bool getUpdateDependentVariableInterpolator( )
    {
        return updateDependentVariableInterpolator_;
    }

    virtual void setClearNumericalSolutions( const bool clearNumericalSolutions )
    {
        clearNumericalSolutions_ = clearNumericalSolutions;
    }

    virtual void setIntegratedResult( const bool setIntegratedResult )
    {
        setIntegratedResult_ = setIntegratedResult;
        createStateProcessors_ = setIntegratedResult;
    }

    virtual void setIntegratedVariationalResult( const bool setIntegratedVariationalResult )
    {
        setIntegratedVariationalResult_ = setIntegratedVariationalResult;
    }
    bool getCreateStateProcessors( )
    {
        return createStateProcessors_;
    }

    virtual void setCreateStateProcessors( const bool createStateProcessors )
    {
        createStateProcessors_ = createStateProcessors;
    }

    virtual void setUpdateDependentVariableInterpolator( const bool updateDependentVariableInterpolator )
    {
        updateDependentVariableInterpolator_ = updateDependentVariableInterpolator;
    }

    virtual void setInterpolatorSettings( const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
    {
        if( interpolatorSettings == nullptr )
        {
            throw std::runtime_error( "Error when setting propagation processing interpolator settings, settings are null." );
        }
        interpolatorSettings_ = interpolatorSettings;
    }

    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( ) const
    {
        return interpolatorSettings_;
    }

    virtual bool printAnyOutput( ) = 0;

    virtual std::string getPropagationStartHeader( ) = 0;

    virtual std::string getPropagationEndHeader( ) = 0;

    // Used for serialization testing
    bool operator==( const PropagatorProcessingSettings& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const PropagatorProcessingSettings& rhs ) const
    {
        return !equals( rhs );
    }

    //! Save processing settings to a file
    TUDAT_DECLARE_FILE_IO_POLYMORPHIC( PropagatorProcessingSettings )

protected:
    bool clearNumericalSolutions_;
    bool setIntegratedResult_;
    bool createStateProcessors_;
    bool updateDependentVariableInterpolator_;
    bool setIntegratedVariationalResult_;
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

protected:
    virtual bool equals( const PropagatorProcessingSettings& rhs ) const
    {
        return clearNumericalSolutions_ == rhs.clearNumericalSolutions_ && setIntegratedResult_ == rhs.setIntegratedResult_ &&
                createStateProcessors_ == rhs.createStateProcessors_ &&
                updateDependentVariableInterpolator_ == rhs.updateDependentVariableInterpolator_ &&
                setIntegratedVariationalResult_ == rhs.setIntegratedVariationalResult_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( clearNumericalSolutions_ ),
            CEREAL_NVP( setIntegratedResult_ ),
            CEREAL_NVP( createStateProcessors_ ),
            CEREAL_NVP( updateDependentVariableInterpolator_ ),
            CEREAL_NVP( setIntegratedVariationalResult_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( clearNumericalSolutions_ ),
            CEREAL_NVP( setIntegratedResult_ ),
            CEREAL_NVP( createStateProcessors_ ),
            CEREAL_NVP( updateDependentVariableInterpolator_ ),
            CEREAL_NVP( setIntegratedVariationalResult_ ) );
    }
};

//! Base class for defining output and processing settings for single-arc propagation.
//! In addition to implementing base class functionality, it defines the output
//! that is to b printed to a terminal during a single-arc propagation (in the printSettings_ member)
class SingleArcPropagatorProcessingSettings : public PropagatorProcessingSettings
{
public:
    SingleArcPropagatorProcessingSettings(
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const int resultsSaveFrequencyInSteps = 1,
            const double resultsSaveFrequencyInSeconds = TUDAT_NAN,
            const std::shared_ptr< PropagationPrintSettings > printSettings = std::make_shared< PropagationPrintSettings >( ),
            const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        resultsSaveFrequencyInSteps_( resultsSaveFrequencyInSteps ), resultsSaveFrequencyInSeconds_( resultsSaveFrequencyInSeconds ),
        printSettings_( printSettings ), isPartOfMultiArc_( false ), arcIndex_( -1 )
    {}
    virtual ~SingleArcPropagatorProcessingSettings( ) {}

    std::shared_ptr< PropagationPrintSettings > getPrintSettings( )
    {
        return printSettings_;
    }

    void setResultsSaveFrequencyInSteps( const int resultsSaveFrequencyInSteps )
    {
        saveWarningPrinted_ = false;
        resultsSaveFrequencyInSteps_ = resultsSaveFrequencyInSteps;
    }

    void setResultsSaveFrequencyInSeconds( const double resultsSaveFrequencyInSeconds )
    {
        saveWarningPrinted_ = false;
        resultsSaveFrequencyInSeconds_ = resultsSaveFrequencyInSeconds;
    }

    int getResultsSaveFrequencyInSteps( )
    {
        return resultsSaveFrequencyInSteps_;
    }

    double getResultsSaveFrequencyInSeconds( )
    {
        return resultsSaveFrequencyInSeconds_;
    }

    void setClearNumericalSolutions( const bool clearNumericalSolutions ) override
    {
        if( isPartOfMultiArc_ )
        {
            throw std::runtime_error(
                    "Error, resetting setClearNumericalSolutions of single-arc propagation setting after it has been made part of a "
                    "multi-arc propagation setting, this is not permitted as it breaks multi-arc settings consistency" );
        }
        PropagatorProcessingSettings::setClearNumericalSolutions( clearNumericalSolutions );
    }

    void setIntegratedResult( const bool setIntegratedResult ) override
    {
        if( isPartOfMultiArc_ )
        {
            throw std::runtime_error(
                    "Error, resetting setIntegratedResult of single-arc propagation setting after it has been made part of a multi-arc "
                    "propagation setting" );
        }
        PropagatorProcessingSettings::setIntegratedResult( setIntegratedResult );
    }

    void setIntegratedVariationalResult( const bool setIntegratedVariationalResult ) override
    {
        if( isPartOfMultiArc_ )
        {
            std::cerr << "Warning, resetting setIntegratedVariationalResult of single-arc propagation setting after it has been made part "
                         "of a multi-arc "
                         "propagation setting has no effect"
                      << std::endl;
        }
        PropagatorProcessingSettings::setIntegratedVariationalResult( setIntegratedVariationalResult );
    }

    void setCreateStateProcessors( const bool createStateProcessors ) override
    {
        if( isPartOfMultiArc_ )
        {
            throw std::runtime_error(
                    "Error, resetting setCreateStateProcessors of single-arc propagation setting after it has been made part of a "
                    "multi-arc propagation setting, this is not permitted as it breaks multi-arc settings consistency" );
        }
        PropagatorProcessingSettings::setCreateStateProcessors( createStateProcessors );
    }

    bool saveCurrentStep( const int stepsSinceLastSave, const double timeSinceLastSave )
    {
        if( !saveWarningPrinted_ &&
            ( resultsSaveFrequencyInSeconds_ == resultsSaveFrequencyInSeconds_ && resultsSaveFrequencyInSteps_ == 1 ) )
        {
            std::cerr << "Warning when saving propagation step. "
                      << "Attribute processing_settings.results_save_frequency_in_steps of propagator settings is set to default value of "
                         "1, while processing_settings.results_save_frequency_in_seconds is also defined. "
                      << "This will result in the results still being saved every time step. To make the "
                         "processing_settings.results_save_frequency_in_seconds active, set "
                         "processing_settings.results_save_frequency_in_steps to 0"
                      << std::endl;
            saveWarningPrinted_ = true;
        }
        bool saveCurrentStep = false;
        if( stepsSinceLastSave >= resultsSaveFrequencyInSteps_ && resultsSaveFrequencyInSteps_ > 0 )
        {
            saveCurrentStep = true;
        }
        else if( timeSinceLastSave >= resultsSaveFrequencyInSeconds_ )
        {
            saveCurrentStep = true;
        }
        return saveCurrentStep;
    }

    bool printAnyOutput( ) override
    {
        return printSettings_->printAnyOutput( );
    }

    std::string getPropagationStartHeader( ) override
    {
        if( isPartOfMultiArc_ )
        {
            return "---------------  STARTING PROPAGATION FOR ARC " + std::to_string( arcIndex_ ) + "  ----------------";
        }
        else
        {
            return "===============  STARTING SINGLE-ARC PROPAGATION  ===============";
        }
    }

    std::string getPropagationEndHeader( ) override
    {
        if( isPartOfMultiArc_ )
        {
            return "-----------------------------------------------------------------";
        }
        else
        {
            return "=================================================================";
        }
    }

private:
    int resultsSaveFrequencyInSteps_;

    double resultsSaveFrequencyInSeconds_;

    std::shared_ptr< PropagationPrintSettings > printSettings_;

    void setAsMultiArc( const unsigned int arcIndex, const bool printArcIndex )
    {
        isPartOfMultiArc_ = true;
        arcIndex_ = arcIndex;
        printSettings_->setPrintArcIndex( printArcIndex );
    }

    void setClearNumericalSolutionsFromMultiArc( const bool clearNumericalSolutions )
    {
        PropagatorProcessingSettings::setClearNumericalSolutions( clearNumericalSolutions );
    }

    void setIntegratedResultFromMultiArc( const bool setIntegratedResult )
    {
        PropagatorProcessingSettings::setIntegratedResult( setIntegratedResult );
    }

    void setIntegratedVariationalResultFromMultiArc( const bool setIntegratedVariationalResult )
    {
        PropagatorProcessingSettings::setIntegratedVariationalResult( setIntegratedVariationalResult );
    }

    void setCreateStateProcessorsFromMultiArc( const bool createStateProcessors )
    {
        PropagatorProcessingSettings::setCreateStateProcessors( createStateProcessors );
    }

    bool isPartOfMultiArc_;
    int arcIndex_;

    bool saveWarningPrinted_ = false;

    bool equals( const PropagatorProcessingSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const SingleArcPropagatorProcessingSettings* >( &other );
        if( !rhs )
        {
            return false;
        }

        if( !PropagatorProcessingSettings::equals( other ) )
        {
            return false;
        }

        const bool printSettingsAreEqual = ( printSettings_ == rhs->printSettings_ ) ||
                ( printSettings_ != nullptr && rhs->printSettings_ != nullptr && ( *printSettings_ == *rhs->printSettings_ ) );

        return resultsSaveFrequencyInSteps_ == rhs->resultsSaveFrequencyInSteps_ &&
                resultsSaveFrequencyInSeconds_ == rhs->resultsSaveFrequencyInSeconds_ && printSettingsAreEqual &&
                isPartOfMultiArc_ == rhs->isPartOfMultiArc_ && arcIndex_ == rhs->arcIndex_ &&
                saveWarningPrinted_ == rhs->saveWarningPrinted_;
    }

    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( resultsSaveFrequencyInSteps_ ),
            CEREAL_NVP( resultsSaveFrequencyInSeconds_ ),
            CEREAL_NVP( printSettings_ ),
            CEREAL_NVP( isPartOfMultiArc_ ),
            CEREAL_NVP( arcIndex_ ),
            CEREAL_NVP( saveWarningPrinted_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( resultsSaveFrequencyInSteps_ ),
            CEREAL_NVP( resultsSaveFrequencyInSeconds_ ),
            CEREAL_NVP( printSettings_ ),
            CEREAL_NVP( isPartOfMultiArc_ ),
            CEREAL_NVP( arcIndex_ ),
            CEREAL_NVP( saveWarningPrinted_ ) );
    }

    friend class MultiArcPropagatorProcessingSettings;
};

template< typename StateScalarType, typename TimeType >
class MultiArcPropagatorSettings;

class MultiArcPropagatorProcessingSettings : public PropagatorProcessingSettings
{
public:
    MultiArcPropagatorProcessingSettings( const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings,
                                          const bool clearNumericalSolutions = false,
                                          const bool setIntegratedResult = false,
                                          const bool printFirstArcOnly = false,
                                          const bool printCurrentArcIndex = false,
                                          const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentSingleArcPrintSettings_( consistentSingleArcPrintSettings ), useIdenticalSettings_( true ),
        printFirstArcOnly_( printFirstArcOnly ), printCurrentArcIndex_( printCurrentArcIndex ), areSingleArcSettingsSet_( false ),
        isPartOfHybridArc_( false )
    {}

    MultiArcPropagatorProcessingSettings( const bool clearNumericalSolutions = false,
                                          const bool setIntegratedResult = false,
                                          const bool printFirstArcOnly = false,
                                          const bool printCurrentArcIndex = false,
                                          const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentSingleArcPrintSettings_( nullptr ), useIdenticalSettings_( false ), printFirstArcOnly_( printFirstArcOnly ),
        printCurrentArcIndex_( printCurrentArcIndex ), areSingleArcSettingsSet_( false ), isPartOfHybridArc_( false )
    {}

    virtual ~MultiArcPropagatorProcessingSettings( ) {}

    void resetSingleArcSettings( const bool printWarning = false )
    {
        if( !areSingleArcSettingsSet_ )
        {
            throw std::runtime_error( "Error in multi-arc output settings, single arc settings not yet defined when resetting" );
        }

        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            singleArcSettings_.at( i )->setClearNumericalSolutionsFromMultiArc( false );
            singleArcSettings_.at( i )->setIntegratedResultFromMultiArc( false );
            singleArcSettings_.at( i )->setIntegratedVariationalResultFromMultiArc( false );
            singleArcSettings_.at( i )->setAsMultiArc( i, printCurrentArcIndex_ );
            singleArcSettings_.at( i )->setCreateStateProcessorsFromMultiArc( setIntegratedResult_ );
            singleArcSettings_.at( i )->setInterpolatorSettings( interpolatorSettings_ );
            if( useIdenticalSettings_ )
            {
                if( consistentSingleArcPrintSettings_ == nullptr )
                {
                    throw std::runtime_error( "Error in multi-arc output settings, no consistent single arc print settings defined" );
                }
                singleArcSettings_.at( i )->getPrintSettings( )->reset( consistentSingleArcPrintSettings_ );
            }

            if( printFirstArcOnly_ && i > 0 )
            {
                singleArcSettings_.at( i )->getPrintSettings( )->disableAllPrinting( );
            }
        }
    }

    void resetConsistentSingleArcPrintSettings( const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings )
    {
        if( useIdenticalSettings_ )
        {
            consistentSingleArcPrintSettings_ = consistentSingleArcPrintSettings;
            resetSingleArcSettings( );
        }
    }

    void resetAndApplyConsistentSingleArcPrintSettings( const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings )
    {
        useIdenticalSettings_ = true;
        resetConsistentSingleArcPrintSettings( consistentSingleArcPrintSettings );
    }

    bool useIdenticalSettings( )
    {
        return useIdenticalSettings_;
    }

    void resetUseIdenticalSettings( const bool useIdenticalSettings )
    {
        useIdenticalSettings_ = useIdenticalSettings;
    }

    void resetPrintCurrentArcIndex( const bool printCurrentArcIndex )
    {
        printCurrentArcIndex_ = printCurrentArcIndex;
        resetSingleArcSettings( );
    }

    bool printAnyOutput( ) override
    {
        bool printOutput = false;
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            if( singleArcSettings_.at( i )->printAnyOutput( ) )
            {
                printOutput = true;
            }
        }
        return printOutput;
    }

    std::string getPropagationStartHeader( ) override
    {
        if( isPartOfHybridArc_ )
        {
            return "- - - - - - - -  STARTING MULTI-ARC PROPAGATION  - - - - - - - -";
        }
        else
        {
            return "===============  STARTING MULTI-ARC PROPAGATION  ===============";
        }
    }

    std::string getPropagationEndHeader( ) override
    {
        if( isPartOfHybridArc_ )
        {
            return "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -";
        }
        else
        {
            return "=================================================================";
        }
    }

    bool getPrintFirstArcOnly( )
    {
        return printFirstArcOnly_;
    }

    void resetPrintFirstArcOnly( const bool printFirstArcOnly )
    {
        printFirstArcOnly_ = printFirstArcOnly;
    }

    std::vector< std::shared_ptr< SingleArcPropagatorProcessingSettings > > getSingleArcSettings( )
    {
        return singleArcSettings_;
    }

    void setIntegratedResult( const bool setIntegratedResult ) override
    {
        setIntegratedResult_ = setIntegratedResult;
        createStateProcessors_ = setIntegratedResult;
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            singleArcSettings_.at( i )->setCreateStateProcessorsFromMultiArc( setIntegratedResult_ );
        }
    }

    void setInterpolatorSettings( const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings ) override
    {
        PropagatorProcessingSettings::setInterpolatorSettings( interpolatorSettings );
        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            singleArcSettings_.at( i )->setInterpolatorSettings( interpolatorSettings );
        }
    }

protected:
    std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings_;

    bool useIdenticalSettings_;

    bool printFirstArcOnly_;

    bool printCurrentArcIndex_;

    std::vector< std::shared_ptr< SingleArcPropagatorProcessingSettings > > singleArcSettings_;

    bool areSingleArcSettingsSet_;

    bool isPartOfHybridArc_;

private:
    void setSingleArcSettings( const std::vector< std::shared_ptr< SingleArcPropagatorProcessingSettings > >& singleArcSettings )
    {
        if( !areSingleArcSettingsSet_ )
        {
            singleArcSettings_ = singleArcSettings;
            areSingleArcSettingsSet_ = true;
            resetSingleArcSettings( );
        }
        else
        {
            throw std::runtime_error(
                    "Error, cannot set constituent single-arc output settings more than once in multi-arc output settings" );
        }
    }

    void setPartOfHybridArc( )
    {
        isPartOfHybridArc_ = true;
    }

    bool equals( const PropagatorProcessingSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const MultiArcPropagatorProcessingSettings* >( &other );
        if( !rhs )
        {
            return false;
        }

        if( !PropagatorProcessingSettings::equals( other ) )
        {
            return false;
        }

        const bool consistentPrintSettingsAreEqual = ( consistentSingleArcPrintSettings_ == rhs->consistentSingleArcPrintSettings_ ) ||
                ( consistentSingleArcPrintSettings_ != nullptr && rhs->consistentSingleArcPrintSettings_ != nullptr &&
                  ( *consistentSingleArcPrintSettings_ == *rhs->consistentSingleArcPrintSettings_ ) );

        if( !consistentPrintSettingsAreEqual || useIdenticalSettings_ != rhs->useIdenticalSettings_ ||
            printFirstArcOnly_ != rhs->printFirstArcOnly_ || printCurrentArcIndex_ != rhs->printCurrentArcIndex_ ||
            areSingleArcSettingsSet_ != rhs->areSingleArcSettingsSet_ || isPartOfHybridArc_ != rhs->isPartOfHybridArc_ ||
            singleArcSettings_ != rhs->singleArcSettings_ )
        {
            return false;
        }

        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            const auto& leftSettings = singleArcSettings_.at( i );
            const auto& rightSettings = rhs->singleArcSettings_.at( i );
            const bool singleArcSettingsAreEqual = ( leftSettings == rightSettings ) ||
                    ( leftSettings != nullptr && rightSettings != nullptr && ( *leftSettings == *rightSettings ) );
            if( !singleArcSettingsAreEqual )
            {
                return false;
            }
        }

        return true;
    }

    template< typename StateScalarType, typename TimeType >
    friend class MultiArcPropagatorSettings;

    friend class HybridArcPropagatorProcessingSettings;

    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( consistentSingleArcPrintSettings_ ),
            CEREAL_NVP( useIdenticalSettings_ ),
            CEREAL_NVP( printFirstArcOnly_ ),
            CEREAL_NVP( printCurrentArcIndex_ ),
            CEREAL_NVP( singleArcSettings_ ),
            CEREAL_NVP( areSingleArcSettingsSet_ ),
            CEREAL_NVP( isPartOfHybridArc_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( consistentSingleArcPrintSettings_ ),
            CEREAL_NVP( useIdenticalSettings_ ),
            CEREAL_NVP( printFirstArcOnly_ ),
            CEREAL_NVP( printCurrentArcIndex_ ),
            CEREAL_NVP( singleArcSettings_ ),
            CEREAL_NVP( areSingleArcSettingsSet_ ),
            CEREAL_NVP( isPartOfHybridArc_ ) );
    }
};

class HybridArcPropagatorProcessingSettings : public PropagatorProcessingSettings
{
public:
    HybridArcPropagatorProcessingSettings( const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings,
                                           const bool clearNumericalSolutions = false,
                                           const bool setIntegratedResult = false,
                                           const bool printStateTypeStart = false,
                                           const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentArcPrintSettings_( consistentArcPrintSettings ), useIdenticalSettings_( true ),
        printStateTypeStart_( printStateTypeStart )
    {}

    HybridArcPropagatorProcessingSettings( const bool clearNumericalSolutions = false,
                                           const bool setIntegratedResult = false,
                                           const bool printStateTypeStart = false,
                                           const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        useIdenticalSettings_( false ), printStateTypeStart_( printStateTypeStart )
    {}

    virtual ~HybridArcPropagatorProcessingSettings( ) {}

    void setClearNumericalSolutions( const bool clearNumericalSolutions ) override
    {
        this->clearNumericalSolutions_ = clearNumericalSolutions;
        singleArcSettings_->setClearNumericalSolutions( clearNumericalSolutions );
        multiArcSettings_->setClearNumericalSolutions( clearNumericalSolutions );
    }

    void setIntegratedResult( const bool setIntegratedResult ) override
    {
        setIntegratedResult_ = setIntegratedResult;
        this->setCreateStateProcessors( setIntegratedResult );
        singleArcSettings_->setIntegratedResult( setIntegratedResult );
        multiArcSettings_->setIntegratedResult( setIntegratedResult );
    }

    void setIntegratedVariationalResult( const bool setIntegratedVariationalResult ) override
    {
        setIntegratedVariationalResult_ = setIntegratedVariationalResult;
        singleArcSettings_->setIntegratedVariationalResult( setIntegratedVariationalResult );
        multiArcSettings_->setIntegratedVariationalResult( setIntegratedVariationalResult );
    }

    void setUpdateDependentVariableInterpolator( const bool updateDependentVariableInterpolator ) override
    {
        this->updateDependentVariableInterpolator_ = updateDependentVariableInterpolator;
        singleArcSettings_->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );
        multiArcSettings_->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );
    }

    void setInterpolatorSettings( const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings ) override
    {
        PropagatorProcessingSettings::setInterpolatorSettings( interpolatorSettings );
        if( singleArcSettings_ != nullptr )
        {
            singleArcSettings_->setInterpolatorSettings( interpolatorSettings );
        }
        if( multiArcSettings_ != nullptr )
        {
            multiArcSettings_->setInterpolatorSettings( interpolatorSettings );
        }
    }

    void resetArcSettings( const bool printWarning = false )
    {
        if( !areArcSettingsSet_ )
        {
            throw std::runtime_error( "Error in hybrid-arc output settings, constitunt arc settings not yet defined when resetting" );
        }

        singleArcSettings_->setClearNumericalSolutions( clearNumericalSolutions_ );
        singleArcSettings_->setIntegratedResult( setIntegratedResult_ );
        singleArcSettings_->setIntegratedVariationalResult( setIntegratedVariationalResult_ );
        singleArcSettings_->setInterpolatorSettings( interpolatorSettings_ );

        multiArcSettings_->setClearNumericalSolutions( clearNumericalSolutions_ );
        multiArcSettings_->setIntegratedResult( setIntegratedResult_ );
        multiArcSettings_->setIntegratedVariationalResult( setIntegratedVariationalResult_ );
        multiArcSettings_->setInterpolatorSettings( interpolatorSettings_ );

        if( useIdenticalSettings_ )
        {
            if( consistentArcPrintSettings_ == nullptr )
            {
                throw std::runtime_error( "Error in bybrid-arc output settings, no consistent arc print settings defined" );
            }
            singleArcSettings_->getPrintSettings( )->reset( consistentArcPrintSettings_ );
            if( !multiArcSettings_->useIdenticalSettings( ) )
            {
                multiArcSettings_->resetUseIdenticalSettings( true );
            }
            multiArcSettings_->resetConsistentSingleArcPrintSettings( consistentArcPrintSettings_ );
        }
    }

    bool printAnyOutput( ) override
    {
        return singleArcSettings_->printAnyOutput( ) || multiArcSettings_->printAnyOutput( );
    }

    std::string getPropagationStartHeader( ) override
    {
        return "==============  STARTING HYBRID-ARC PROPAGATION  ===============";
    }

    std::string getPropagationEndHeader( ) override
    {
        return "=================================================================";
    }

    std::shared_ptr< SingleArcPropagatorProcessingSettings > getSingleArcSettings( )
    {
        return singleArcSettings_;
    }

    std::shared_ptr< MultiArcPropagatorProcessingSettings > getMultiArcSettings( )
    {
        return multiArcSettings_;
    }

    void resetConsistentPrintSettings( const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings )
    {
        if( useIdenticalSettings_ )
        {
            consistentArcPrintSettings_ = consistentArcPrintSettings;
            multiArcSettings_->resetSingleArcSettings( );
            singleArcSettings_->getPrintSettings( )->reset( consistentArcPrintSettings );
        }
    }

    void resetAndApplyConsistentPrintSettings( const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings )
    {
        useIdenticalSettings_ = true;
        resetConsistentPrintSettings( consistentArcPrintSettings );
    }

protected:
    std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings_;

    bool useIdenticalSettings_;

    bool printStateTypeStart_;

    std::shared_ptr< SingleArcPropagatorProcessingSettings > singleArcSettings_ = nullptr;

    std::shared_ptr< MultiArcPropagatorProcessingSettings > multiArcSettings_ = nullptr;

    bool areArcSettingsSet_ = false;

private:
    void setSingleArcSettings( const std::shared_ptr< SingleArcPropagatorProcessingSettings > singleArcSettings,
                               const std::shared_ptr< MultiArcPropagatorProcessingSettings > multiArcSettings )
    {
        if( !areArcSettingsSet_ )
        {
            singleArcSettings_ = singleArcSettings;
            multiArcSettings_ = multiArcSettings;
            multiArcSettings_->setPartOfHybridArc( );
            areArcSettingsSet_ = true;
            resetArcSettings( );
        }
        else
        {
            throw std::runtime_error(
                    "Error, cannot set constituent single-arc output settings more than once in multi-arc output settings" );
        }
    }

    bool equals( const PropagatorProcessingSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const HybridArcPropagatorProcessingSettings* >( &other );
        if( !rhs )
        {
            return false;
        }

        if( !PropagatorProcessingSettings::equals( other ) )
        {
            return false;
        }

        const bool consistentPrintSettingsAreEqual = ( consistentArcPrintSettings_ == rhs->consistentArcPrintSettings_ ) ||
                ( consistentArcPrintSettings_ != nullptr && rhs->consistentArcPrintSettings_ != nullptr &&
                  ( *consistentArcPrintSettings_ == *rhs->consistentArcPrintSettings_ ) );
        const bool singleArcSettingsAreEqual = ( singleArcSettings_ == rhs->singleArcSettings_ ) ||
                ( singleArcSettings_ != nullptr && rhs->singleArcSettings_ != nullptr &&
                  ( *singleArcSettings_ == *rhs->singleArcSettings_ ) );
        const bool multiArcSettingsAreEqual = ( multiArcSettings_ == rhs->multiArcSettings_ ) ||
                ( multiArcSettings_ != nullptr && rhs->multiArcSettings_ != nullptr && ( *multiArcSettings_ == *rhs->multiArcSettings_ ) );

        return consistentPrintSettingsAreEqual && useIdenticalSettings_ == rhs->useIdenticalSettings_ &&
                printStateTypeStart_ == rhs->printStateTypeStart_ && singleArcSettingsAreEqual && multiArcSettingsAreEqual &&
                areArcSettingsSet_ == rhs->areArcSettingsSet_;
    }

    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( consistentArcPrintSettings_ ),
            CEREAL_NVP( useIdenticalSettings_ ),
            CEREAL_NVP( printStateTypeStart_ ),
            CEREAL_NVP( singleArcSettings_ ),
            CEREAL_NVP( multiArcSettings_ ),
            CEREAL_NVP( areArcSettingsSet_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< PropagatorProcessingSettings >( this ) );
        ar( CEREAL_NVP( consistentArcPrintSettings_ ),
            CEREAL_NVP( useIdenticalSettings_ ),
            CEREAL_NVP( printStateTypeStart_ ),
            CEREAL_NVP( singleArcSettings_ ),
            CEREAL_NVP( multiArcSettings_ ),
            CEREAL_NVP( areArcSettingsSet_ ) );
    }

    template< typename StateScalarType, typename TimeType >
    friend class HybridArcPropagatorSettings;
};

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_PROPAGATIONPROCESSINGSETTINGS_H
