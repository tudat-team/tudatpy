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
    PropagatorProcessingSettings(
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool updateDependentVariableInterpolator = false ):
        clearNumericalSolutions_( clearNumericalSolutions ),
        setIntegratedResult_( setIntegratedResult ),
        updateDependentVariableInterpolator_( updateDependentVariableInterpolator )
    { }

    virtual ~PropagatorProcessingSettings( ){ }

    bool getClearNumericalSolutions( )
    {
        return clearNumericalSolutions_;
    }

    bool getSetIntegratedResult( )
    {
        return setIntegratedResult_;
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
    }

    virtual void setUpdateDependentVariableInterpolator( const bool updateDependentVariableInterpolator )
    {
         updateDependentVariableInterpolator_ = updateDependentVariableInterpolator;
    }

    virtual bool printAnyOutput( ) = 0;

    virtual std::string getPropagationStartHeader( ) = 0;

    virtual std::string getPropagationEndHeader( ) = 0;

protected:

    bool clearNumericalSolutions_;
    bool setIntegratedResult_;
    bool updateDependentVariableInterpolator_;
};

//! Base class for defining output and processing settings for single-arc propagation.
//! In addition to implementing base class functionality, it defines the output
//! that is to b printed to a terminal during a single-arc propagation (in the printSettings_ member)
class SingleArcPropagatorProcessingSettings: public PropagatorProcessingSettings
{
public:

    SingleArcPropagatorProcessingSettings(
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const int resultsSaveFrequencyInSteps = 1,
            const double resultsSaveFrequencyInSeconds = TUDAT_NAN,
            const std::shared_ptr< PropagationPrintSettings > printSettings =
            std::make_shared< PropagationPrintSettings >( ),
            const bool updateDependentVariableInterpolator = false ):
            PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
            resultsSaveFrequencyInSteps_( resultsSaveFrequencyInSteps ),
            resultsSaveFrequencyInSeconds_( resultsSaveFrequencyInSeconds ),
            printSettings_( printSettings ),
        isPartOfMultiArc_( false ), arcIndex_( -1 ){ }
    virtual ~SingleArcPropagatorProcessingSettings( ){ }

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

    bool saveCurrentStep(
            const int stepsSinceLastSave, const double timeSinceLastSave )
    {
        if( !saveWarningPrinted_ && ( resultsSaveFrequencyInSeconds_ == resultsSaveFrequencyInSeconds_ && resultsSaveFrequencyInSteps_ == 1 ) )
        {
            std::cerr<<"Warning when saving propagation step. "<<
            "Attribute processing_settings.results_save_frequency_in_steps of propagator settings is set to default value of 1, while processing_settings.results_save_frequency_in_seconds is also defined. "<<
            "This will result in the results still being saved every time step. To make the processing_settings.results_save_frequency_in_seconds active, set processing_settings.results_save_frequency_in_steps to 0"<<std::endl;
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



    bool printAnyOutput( )
    {
        return printSettings_->printAnyOutput( );
    }


    std::string getPropagationStartHeader( )
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

    std::string getPropagationEndHeader( )
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

    const std::shared_ptr< PropagationPrintSettings > printSettings_;

    void setAsMultiArc( const unsigned int arcIndex, const bool printArcIndex )
    {
        isPartOfMultiArc_ = true;
        arcIndex_ = arcIndex;
        printSettings_->setPrintArcIndex( printArcIndex );
    }

    bool isPartOfMultiArc_;
    int arcIndex_;

    bool saveWarningPrinted_ = false;

    friend class MultiArcPropagatorProcessingSettings;
};

template< typename StateScalarType, typename TimeType >
class MultiArcPropagatorSettings;

class MultiArcPropagatorProcessingSettings: public PropagatorProcessingSettings
{
public:
    MultiArcPropagatorProcessingSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printFirstArcOnly = false,
            const bool printCurrentArcIndex = false,
            const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentSingleArcPrintSettings_( consistentSingleArcPrintSettings ),
        useIdenticalSettings_( true ),
        printFirstArcOnly_( printFirstArcOnly ),
        printCurrentArcIndex_( printCurrentArcIndex ),
        areSingleArcSettingsSet_( false ),
        isPartOfHybridArc_( false )
    {
    }

    MultiArcPropagatorProcessingSettings(
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printFirstArcOnly = false,
            const bool printCurrentArcIndex = false,
            const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentSingleArcPrintSettings_( nullptr ),
        useIdenticalSettings_( false ),
        printFirstArcOnly_( printFirstArcOnly ),
        printCurrentArcIndex_( printCurrentArcIndex ),
        areSingleArcSettingsSet_( false ),
        isPartOfHybridArc_( false )
    {
    }

    virtual ~MultiArcPropagatorProcessingSettings( ){ }


    void resetSingleArcSettings( const bool printWarning = false )
    {
        if( !areSingleArcSettingsSet_ )
        {
            throw std::runtime_error( "Error in multi-arc output settings, single arc settings not yet defined when resetting" );
        }

        for( unsigned int i = 0; i < singleArcSettings_.size( ); i++ )
        {
            singleArcSettings_.at( i )->setClearNumericalSolutions( false );
            singleArcSettings_.at( i )->setIntegratedResult( false );
            singleArcSettings_.at( i )->setAsMultiArc( i, printCurrentArcIndex_ );

            if( useIdenticalSettings_ )
            {
                if( consistentSingleArcPrintSettings_ == nullptr )
                {
                    throw std::runtime_error( "Error in multi-arc output settings, no consistent single arc print settings defined" );
                }
                singleArcSettings_.at( i )->getPrintSettings( )->reset(
                            consistentSingleArcPrintSettings_ );
            }

            if( printFirstArcOnly_ && i > 0 )
            {
                singleArcSettings_.at( i )->getPrintSettings( )->disableAllPrinting( );
            }

        }
    }



    void resetConsistentSingleArcPrintSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings )
    {
        if( useIdenticalSettings_ )
        {
            consistentSingleArcPrintSettings_ = consistentSingleArcPrintSettings;
            resetSingleArcSettings( );
        }
    }

    void resetAndApplyConsistentSingleArcPrintSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings )
    {
        useIdenticalSettings_ = true;
        resetConsistentSingleArcPrintSettings( consistentSingleArcPrintSettings );
    }

    bool useIdenticalSettings( )
    {
        return useIdenticalSettings_;
    }

    void resetUseIdenticalSettings(
            const bool useIdenticalSettings )
    {
        useIdenticalSettings_ = useIdenticalSettings;
    }



    void resetPrintCurrentArcIndex(
            const bool printCurrentArcIndex )
    {
        printCurrentArcIndex_ = printCurrentArcIndex;
        resetSingleArcSettings( );
    }

    bool printAnyOutput( )
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

    std::string getPropagationStartHeader( )
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

    std::string getPropagationEndHeader( )
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


protected:

    std::shared_ptr< PropagationPrintSettings > consistentSingleArcPrintSettings_;

    bool useIdenticalSettings_;

    bool printFirstArcOnly_;

    bool printCurrentArcIndex_;

    std::vector< std::shared_ptr< SingleArcPropagatorProcessingSettings > > singleArcSettings_;

    bool areSingleArcSettingsSet_;

    bool isPartOfHybridArc_;

private:


    void setSingleArcSettings(
            const std::vector< std::shared_ptr< SingleArcPropagatorProcessingSettings > >& singleArcSettings )
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

    template< typename StateScalarType, typename TimeType >
    friend class MultiArcPropagatorSettings;

    friend class HybridArcPropagatorProcessingSettings;

};



class HybridArcPropagatorProcessingSettings: public PropagatorProcessingSettings
{
public:
    HybridArcPropagatorProcessingSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printStateTypeStart = false,
            const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        consistentArcPrintSettings_( consistentArcPrintSettings ),
        useIdenticalSettings_( true ),
        printStateTypeStart_( printStateTypeStart ){ }

    HybridArcPropagatorProcessingSettings(
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printStateTypeStart = false,
            const bool updateDependentVariableInterpolator = false ):
        PropagatorProcessingSettings( clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
        useIdenticalSettings_( false ),
        printStateTypeStart_( printStateTypeStart ){ }

    virtual ~HybridArcPropagatorProcessingSettings( ){ }

    virtual void setClearNumericalSolutions( const bool clearNumericalSolutions )
    {
        this->clearNumericalSolutions_ = clearNumericalSolutions;
        singleArcSettings_->setClearNumericalSolutions( clearNumericalSolutions );
        multiArcSettings_->setClearNumericalSolutions( clearNumericalSolutions );
    }

    virtual void setIntegratedResult( const bool setIntegratedResult )
    {
        this->setIntegratedResult_ = setIntegratedResult;
        singleArcSettings_->setIntegratedResult( setIntegratedResult );
        multiArcSettings_->setIntegratedResult( setIntegratedResult );
    }

    virtual void setUpdateDependentVariableInterpolator( const bool updateDependentVariableInterpolator )
    {
        this->updateDependentVariableInterpolator_ = updateDependentVariableInterpolator;
        singleArcSettings_->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );
        multiArcSettings_->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );
    }

    void resetArcSettings( const bool printWarning = false )
    {
        if( !areArcSettingsSet_ )
        {
            throw std::runtime_error( "Error in hybrid-arc output settings, constitunt arc settings not yet defined when resetting" );
        }

        singleArcSettings_->setClearNumericalSolutions( clearNumericalSolutions_ );
        singleArcSettings_->setIntegratedResult( setIntegratedResult_ );

        multiArcSettings_->setClearNumericalSolutions( clearNumericalSolutions_ );
        multiArcSettings_->setIntegratedResult( setIntegratedResult_ );

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

    bool printAnyOutput( )
    {
        return singleArcSettings_->printAnyOutput( ) || multiArcSettings_->printAnyOutput( );
    }

    std::string getPropagationStartHeader( )
    {
        return "==============  STARTING HYBRID-ARC PROPAGATION  ===============";
    }

    std::string getPropagationEndHeader( )
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

    void resetConsistentPrintSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings )
    {
        if( useIdenticalSettings_ )
        {
            consistentArcPrintSettings_ = consistentArcPrintSettings;
            multiArcSettings_->resetSingleArcSettings( );
            singleArcSettings_->getPrintSettings( )->reset( consistentArcPrintSettings );
        }
    }

    void resetAndApplyConsistentPrintSettings(
            const std::shared_ptr< PropagationPrintSettings > consistentArcPrintSettings )
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


    void setSingleArcSettings(
            const std::shared_ptr< SingleArcPropagatorProcessingSettings > singleArcSettings,
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

    template< typename StateScalarType, typename TimeType >
    friend class HybridArcPropagatorSettings;

};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONPROCESSINGSETTINGS_H
