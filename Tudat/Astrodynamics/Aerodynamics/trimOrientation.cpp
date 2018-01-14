/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Aerodynamics/trimOrientation.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"
#include "Tudat/Mathematics/RootFinders/secantRootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"
namespace tudat
{

namespace aerodynamics
{

//! Constructor
TrimOrientationCalculator::TrimOrientationCalculator(
        const boost::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface,
        const boost::shared_ptr< root_finders::RootFinderCore< double > > rootFinder ):
    coefficientInterface_( coefficientInterface ), rootFinder_( rootFinder )
{
    // Find index of angle of attack in aerodynamic coefficient interface (throw error if not found)
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariables =
            coefficientInterface->getIndependentVariableNames( );
    std::vector< AerodynamicCoefficientsIndependentVariables >::iterator variableIterator =
            std::find( independentVariables.begin( ), independentVariables.end( ), angle_of_attack_dependent );

    if( variableIterator == independentVariables.end( ) )
    {
        throw std::runtime_error( "Error when getting trim angle of attack, no angle of attack dependency is found" );
    }
    variableIndex_ = std::distance( independentVariables.begin( ), variableIterator );

    // Find indices of angle of attack in control surface interfaces.
    std::map< std::string, std::vector< AerodynamicCoefficientsIndependentVariables > > controlSurfaceIndependentVariables =
            coefficientInterface->getControlSurfaceIndependentVariables( );
    for( std::map< std::string, std::vector< AerodynamicCoefficientsIndependentVariables > >::iterator
         controlSurfaceIterator = controlSurfaceIndependentVariables.begin( );
         controlSurfaceIterator != controlSurfaceIndependentVariables.end( ); controlSurfaceIterator++ )
    {
        std::vector< AerodynamicCoefficientsIndependentVariables > currentIndependentVariables =
                controlSurfaceIterator->second;
        variableIterator =
                std::find(
                    currentIndependentVariables.begin( ), currentIndependentVariables.end( ), angle_of_attack_dependent );
        if( ( variableIterator == currentIndependentVariables.end( ) ) )
        {
            std::cerr << "Warning when getting trim angle of attack, no angle of attack dependency is found for control surface "
                       + controlSurfaceIterator->first << std::endl;
        }
        else
        {
            controlSurfaceVariableIndex_[ controlSurfaceIterator->first ] =
                    std::distance( currentIndependentVariables.begin( ), variableIterator );
        }
    }

    // If no root finder provided, use default value.
    if ( !rootFinder_.get( ) )
    {
        rootFinder_ = boost::make_shared< root_finders::SecantRootFinderCore< double > >(
                    boost::bind(
                        &root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition< double >::
                        checkTerminationCondition,
                        boost::make_shared< root_finders::termination_conditions::RootAbsoluteToleranceTerminationCondition
                        < double > >( 1.0E-15, 1000 ), _1, _2, _3, _4, _5 ), 0.5 );
    }
}

//! Function to find the trimmed angle of attack for a given set of independent  variables
double TrimOrientationCalculator::findTrimAngleOfAttack(
        const std::vector< double > untrimmedIndependentVariables,
        const std::map< std::string, std::vector< double > > untrimmedControlSurfaceIndependentVariables )
{
    // Determine function for which the root is to be determined.
    boost::function< double( const double ) > coefficientFunction =
            boost::bind( &TrimOrientationCalculator::getPerturbedMomentCoefficient,
                         this, _1, untrimmedIndependentVariables, untrimmedControlSurfaceIndependentVariables );

    double trimmedAngleOfAttack = TUDAT_NAN;

    // Find root of pitch moment function
    try
    {
        trimmedAngleOfAttack = rootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >( coefficientFunction ),
                    untrimmedIndependentVariables.at( variableIndex_ ) );
    }
    // Throw error if not converged
    catch( std::runtime_error )
    {
        throw std::runtime_error( "Error when getting trim angle of attack, root finder did not converge." );

    }

    return trimmedAngleOfAttack;
}


//! Function to get the moment coefficient for a given angle of attack
double TrimOrientationCalculator::getPerturbedMomentCoefficient(
        const double perturbedAngleOfAttack,
        const std::vector< double >& unperturbedConditions,
        const std::map< std::string, std::vector< double > > unperturbedControlSurfaceIndependentVariables )
{
    // Update coefficients to perturbed independent variables
    std::vector< double > perturbedConditions = unperturbedConditions;
    std::map< std::string, std::vector< double > > perturbedControlSurfaceConditions =
            unperturbedControlSurfaceIndependentVariables;

    perturbedConditions[ variableIndex_ ] = perturbedAngleOfAttack;
    for( std::map< std::string, int >::const_iterator controlSurfaceIterator = controlSurfaceVariableIndex_.begin( );
         controlSurfaceIterator != controlSurfaceVariableIndex_.end( ); controlSurfaceIterator++ )
    {
        if( perturbedControlSurfaceConditions.count( controlSurfaceIterator->first ) > 0 )
        {
            perturbedControlSurfaceConditions[ controlSurfaceIterator->first ][ controlSurfaceIterator->second ] =
                    perturbedAngleOfAttack;
        }
    }

    coefficientInterface_->updateFullCurrentCoefficients( perturbedConditions, perturbedControlSurfaceConditions );

    // Get pitch moment coefficient
    return coefficientInterface_->getCurrentMomentCoefficients( )( 1 );
}

}

}


