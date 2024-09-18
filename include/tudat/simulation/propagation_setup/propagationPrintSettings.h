/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONPRINTSETTINGS_H
#define TUDAT_PROPAGATIONPRINTSETTINGS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace propagators
{


//! Class defining settings for output written to cout (terminal) during a propagation
class PropagationPrintSettings
{
public:
    PropagationPrintSettings(
            const bool printNumberOfFunctionEvaluations = false,
            const bool printDependentVariableData = false,
            const double resultsPrintFrequencyInSeconds = TUDAT_NAN,
            const int resultsPrintFrequencyInSteps = 0,
            const bool printTerminationReason = false,
            const bool printPropagationTime = false,
            const bool printPropagatedStateData = false,
            const bool printInitialAndFinalConditions = false,
            const bool printDependentVariableDuringPropagation = false,
            const bool printProcessedStateData = false ): printArcIndex_( false )
    {
        reset( printNumberOfFunctionEvaluations,
               printDependentVariableData, resultsPrintFrequencyInSeconds, resultsPrintFrequencyInSteps,
               printTerminationReason, printPropagationTime, printPropagatedStateData,
               printInitialAndFinalConditions, printDependentVariableDuringPropagation, printProcessedStateData );
    }

    bool getPrintNumberOfFunctionEvaluations( ) { return printNumberOfFunctionEvaluations_; }

    void setPrintNumberOfFunctionEvaluations( const bool printNumberOfFunctionEvaluations )
    { printNumberOfFunctionEvaluations_ = printNumberOfFunctionEvaluations; }


    bool getPrintDependentVariableData( ) {  return printDependentVariableData_;  }

    void setPrintDependentVariableData( const bool printDependentVariableData )
    {  printDependentVariableData_ = printDependentVariableData; }


    bool getPrintProcessedStateData( ) {  return printProcessedStateData_;  }

    void setPrintProcessedStateData( const bool printProcessedStateData )
    {  printProcessedStateData_ = printProcessedStateData; }


    bool getPrintPropagatedStateData( ){ return printPropagatedStateData_; }

    void setPrintPropagatedStateData( const bool printPropagatedStateData )
    { printPropagatedStateData_ = printPropagatedStateData; }
    

    double getResultsPrintFrequencyInSeconds( ){ return resultsPrintFrequencyInSeconds_; }

    void setResultsPrintFrequencyInSeconds( const double resultsPrintFrequencyInSeconds )
    { resultsPrintFrequencyInSeconds_ = resultsPrintFrequencyInSeconds; }


    double getResultsPrintFrequencyInSteps( ){ return resultsPrintFrequencyInSteps_; }

    void setResultsPrintFrequencyInSteps( const double resultsPrintFrequencyInSteps )
    { resultsPrintFrequencyInSteps_ = resultsPrintFrequencyInSteps; }


    bool getPrintTerminationReason( ) { return printTerminationReason_;  }

    void setPrintTerminationReason( const bool printTerminationReason ) 
    { printTerminationReason_ = printTerminationReason; }


    bool getPrintPropagationTime( ) { return printPropagationTime_; }

    void setPrintPropagationTime( const bool printPropagationTime ) 
    { printPropagationTime_ = printPropagationTime; }


    bool getPrintInitialAndFinalConditions( ){ return printInitialAndFinalConditions_; }

    void setPrintInitialAndFinalConditions( const bool printInitialAndFinalConditions )
    { printInitialAndFinalConditions_ = printInitialAndFinalConditions; }


    bool getPrintDependentVariableDuringPropagation( ){ return printDependentVariableDuringPropagation_; }

    void setPrintDependentVariableDuringPropagation( const bool printDependentVariableDuringPropagation )
    { printDependentVariableDuringPropagation_ = printDependentVariableDuringPropagation; }


    bool printCurrentStep(
            const int stepsSinceLastPrint, const double timeSinceLastPrint )
    {
        bool printCurrentStep = false;
//        std::cout<<"Settings: "<<resultsPrintFrequencyInSeconds_<<" "<<resultsPrintFrequencyInSteps_<<std::endl;
//        std::cout<<"Input: "<<timeSinceLastPrint<<" "<<stepsSinceLastPrint<<std::endl;
        if( stepsSinceLastPrint >= resultsPrintFrequencyInSteps_ && resultsPrintFrequencyInSteps_ > 0 )
        {
            printCurrentStep = true;
        }
        else if( timeSinceLastPrint >= resultsPrintFrequencyInSeconds_ )
        {
            printCurrentStep = true;
        }
        else if( resultsPrintFrequencyInSeconds_ == resultsPrintFrequencyInSteps_ && !( timeSinceLastPrint == timeSinceLastPrint ) )
        {
            printCurrentStep = true;
        }
        return printCurrentStep;
    }
    
    void setPrintArcIndex( const bool printArcIndex )
    {
        printArcIndex_ = printArcIndex;
    }

    // Check if any output is to be printed before propagation
    bool printPostPropagation( )
    {
        return ( printNumberOfFunctionEvaluations_ || printTerminationReason_ || printPropagationTime_ || printInitialAndFinalConditions_ || printProcessedStateData_ );
    }

    // Check if any output is to be printed during propagation
    bool printDuringPropagation( )
    {
        return ( ( resultsPrintFrequencyInSeconds_ == resultsPrintFrequencyInSeconds_ ) || ( resultsPrintFrequencyInSteps_ > 0 ) );
    }

    // Check if any output is to be printed after propagation
    bool printBeforePropagation( )
    {
        return ( printPropagatedStateData_ || printDependentVariableData_ || printArcIndex_ );
    }

    bool printAnyOutput( )
    {
        return ( printPostPropagation( ) ||
                 printDuringPropagation( ) ||
                 printBeforePropagation( ) );
    }

    void reset( const bool printNumberOfFunctionEvaluations,
                const bool printDependentVariableData,
                const double resultsPrintFrequencyInSeconds,
                const int resultsPrintFrequencyInSteps,
                const bool printTerminationReason,
                const bool printPropagationTime,
                const bool printPropagatedStateData,
                const bool printInitialAndFinalConditions,
                const bool printDependentVariableDuringPropagation,
                const bool printProcessedStateData )
    {
        printNumberOfFunctionEvaluations_ =  printNumberOfFunctionEvaluations;
        printDependentVariableData_ =  printDependentVariableData;
        resultsPrintFrequencyInSeconds_ = resultsPrintFrequencyInSeconds;
        resultsPrintFrequencyInSteps_ = resultsPrintFrequencyInSteps;
        printTerminationReason_ = printTerminationReason;
        printPropagationTime_ = printPropagationTime;
        printPropagatedStateData_ = printPropagatedStateData;
        printInitialAndFinalConditions_ = printInitialAndFinalConditions;
        printDependentVariableDuringPropagation_ = printDependentVariableDuringPropagation;
        printProcessedStateData_ = printProcessedStateData;

    }

    void reset( const std::shared_ptr< PropagationPrintSettings > printSettings )
    {
        printNumberOfFunctionEvaluations_ =  printSettings->getPrintNumberOfFunctionEvaluations( );
        printDependentVariableData_ =  printSettings->getPrintDependentVariableData( );
        resultsPrintFrequencyInSeconds_ = printSettings->getResultsPrintFrequencyInSeconds( );
        resultsPrintFrequencyInSteps_ = printSettings->getResultsPrintFrequencyInSteps( );
        printTerminationReason_ = printSettings->getPrintTerminationReason( );
        printPropagationTime_ = printSettings->getPrintPropagationTime( );
        printPropagatedStateData_ = printSettings->getPrintPropagatedStateData( );
        printInitialAndFinalConditions_ = printSettings->getPrintInitialAndFinalConditions( );
        printProcessedStateData_ = printSettings->getPrintProcessedStateData( );

    }

    // Print nothing
    void disableAllPrinting( )
    {
        reset( false, false, TUDAT_NAN, 0, false, false, false, false, false, false );
    }

    // Print everything, but keep print interval during propagation the same
    void enableAllPrinting( )
    {
        reset( true, true, resultsPrintFrequencyInSeconds_, resultsPrintFrequencyInSteps_, true, true, true, true, true, true );
    }

    // Print everything, and reset print interval during propagation
    void enableAllPrinting( const double resultsPrintFrequencyInSeconds, const int resultsPrintFrequencyInSteps )
    {
        reset( true, true, resultsPrintFrequencyInSeconds, resultsPrintFrequencyInSteps, true, true, true, true, true, true );
    }

private:

    bool printNumberOfFunctionEvaluations_;
    bool printDependentVariableData_;
    double resultsPrintFrequencyInSeconds_;
    int resultsPrintFrequencyInSteps_;
    bool printTerminationReason_;
    bool printPropagationTime_;
    bool printPropagatedStateData_;
    bool printInitialAndFinalConditions_;
    bool printDependentVariableDuringPropagation_;
    bool printProcessedStateData_;

    bool printArcIndex_;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONPRINTSETTINGS_H
