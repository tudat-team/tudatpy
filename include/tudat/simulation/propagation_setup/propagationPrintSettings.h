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
            const bool printStateData = false,
            const double statePrintInterval = TUDAT_NAN,
            const bool printTerminationReason = false,
            const bool printPropagationTime = false,
            const bool printPropagatedStateData = false,
            const bool printInitialAndFinalConditions = false,
            const bool printDependentVariableDuringPropagation = false ): printArcIndex_( false )
    {
        reset( printNumberOfFunctionEvaluations,
               printDependentVariableData, printStateData, statePrintInterval,
               printTerminationReason, printPropagationTime, printPropagatedStateData,
               printInitialAndFinalConditions, printDependentVariableDuringPropagation );
    }

    bool getPrintNumberOfFunctionEvaluations( ) { return printNumberOfFunctionEvaluations_; }

    void setPrintNumberOfFunctionEvaluations( const bool printNumberOfFunctionEvaluations )
    { printNumberOfFunctionEvaluations_ = printNumberOfFunctionEvaluations; }


    bool getPrintDependentVariableData( ) {  return printDependentVariableData_;  }

    void setPrintDependentVariableData( const bool printDependentVariableData )
    {  printDependentVariableData_ = printDependentVariableData; }

    bool getPrintStateData( ) { return printStateData_; }

    void setPrintStateData( const bool printStateData ) { printStateData_ = printStateData; }


    bool getPrintPropagatedStateData( ){ return printPropagatedStateData_; }

    void setPrintPropagatedStateData( const bool printPropagatedStateData )
    { printPropagatedStateData_ = printPropagatedStateData; }

    double getStatePrintInterval( ){ return statePrintInterval_; }

    void setStatePrintInterval( const double statePrintInterval )
    { statePrintInterval_ = statePrintInterval; }


    bool getPrintTerminationReason( ) { return printTerminationReason_;  }

    void setPrintTerminationReason( const bool printTerminationReason ) { printTerminationReason_ = printTerminationReason; }



    bool getPrintPropagationTime( ) { return printPropagationTime_; }

    void setPrintPropagationTime( const bool printPropagationTime ) { printPropagationTime_ = printPropagationTime; }


    bool getPrintInitialAndFinalConditions( ){ return printInitialAndFinalConditions_; }

    void setPrintInitialAndFinalConditions( const bool printInitialAndFinalConditions )
    { printInitialAndFinalConditions_ = printInitialAndFinalConditions; }


    bool getPrintDependentVariableDuringPropagation( ){ return printDependentVariableDuringPropagation_; }

    void setPrintDependentVariableDuringPropagation_( const bool printDependentVariableDuringPropagation )
    { printDependentVariableDuringPropagation_ = printDependentVariableDuringPropagation; }




    void setPrintArcIndex( const bool printArcIndex )
    {
        printArcIndex_ = printArcIndex;
    }

    // Check if any output is to be printed before propagation
    bool printPostPropagation( )
    {
        return ( printNumberOfFunctionEvaluations_ || printTerminationReason_ || printPropagationTime_ || printInitialAndFinalConditions_ );
    }

    // Check if any output is to be printed during propagation
    bool printDuringPropagation( )
    {
        return ( ( statePrintInterval_ == statePrintInterval_ ) );
    }

    // Check if any output is to be printed after propagation
    bool printBeforePropagation( )
    {
        return ( printStateData_ || printPropagatedStateData_ || printDependentVariableData_ || printArcIndex_ );
    }

    void reset( const bool printNumberOfFunctionEvaluations,
                const bool printDependentVariableData,
                const bool printStateData,
                const double statePrintInterval,
                const bool printTerminationReason,
                const bool printPropagationTime,
                const bool printPropagatedStateData,
                const bool printInitialAndFinalConditions,
                const bool printDependentVariableDuringPropagation )
    {
        printNumberOfFunctionEvaluations_ =  printNumberOfFunctionEvaluations;
        printDependentVariableData_ =  printDependentVariableData;
        printStateData_ = printStateData;
        statePrintInterval_ = statePrintInterval;
        printTerminationReason_ = printTerminationReason;
        printPropagationTime_ = printPropagationTime;
        printPropagatedStateData_ = printPropagatedStateData;
        printInitialAndFinalConditions_ = printInitialAndFinalConditions;
        printDependentVariableDuringPropagation_ = printDependentVariableDuringPropagation;

    }

    void reset( const std::shared_ptr< PropagationPrintSettings > printSettings )
    {
        printNumberOfFunctionEvaluations_ =  printSettings->getPrintNumberOfFunctionEvaluations( );
        printDependentVariableData_ =  printSettings->getPrintDependentVariableData( );
        printStateData_ = printSettings->getPrintStateData( );
        statePrintInterval_ = printSettings->getStatePrintInterval( );
        printTerminationReason_ = printSettings->getPrintTerminationReason( );
        printPropagationTime_ = printSettings->getPrintPropagationTime( );
        printPropagatedStateData_ = printSettings->getPrintPropagatedStateData( );
        printInitialAndFinalConditions_ = printSettings->getPrintInitialAndFinalConditions( );

    }

    // Print nothing
    void disableAllPrinting( )
    {
        reset( false, false, false, TUDAT_NAN, false, false, false, false, false );
    }

    // Print everything, but keep print interval during propagation the same
    void enableAllPrinting( )
    {
        reset( true, true, true, statePrintInterval_, true, true, true, true, true );
    }

    // Print everything, and reset print interval during propagation
    void enableAllPrinting( const double statePrintInterval )
    {
        reset( true, true, true, statePrintInterval, true, true, true, true, true );
    }

private:

    bool printNumberOfFunctionEvaluations_;
    bool printDependentVariableData_;
    bool printStateData_;
    double statePrintInterval_;
    bool printTerminationReason_;
    bool printPropagationTime_;
    bool printPropagatedStateData_;
    bool printInitialAndFinalConditions_;
    bool printDependentVariableDuringPropagation_;

    bool printArcIndex_;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONPRINTSETTINGS_H
