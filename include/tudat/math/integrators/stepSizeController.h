/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *
 */

#ifndef TUDAT_STEPSIZE_CONTROLLER_H
#define TUDAT_STEPSIZE_CONTROLLER_H




#include <functional>
#include <memory>

#include <Eigen/Core>

#include <limits>
#include <vector>

#include "tudat/basics/utilityMacros.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{

namespace numerical_integrators
{

template< typename TimeStepType >
class IntegratorStepSizeValidator
{
public:
    IntegratorStepSizeValidator(
        const bool acceptInfinityStep = false,
        const bool acceptNanStep = false
        ):acceptInfinityStep_( acceptInfinityStep ), acceptNanStep_( acceptNanStep ){ }

    virtual ~IntegratorStepSizeValidator( ){ }

    virtual std::pair< TimeStepType, bool > validateStep(
        const std::pair< TimeStepType, bool > recommendedStep, const double currentStep ) = 0;

    void resetAcceptInfinityStep( const bool acceptInfinityStep )
    {
        acceptInfinityStep_ = acceptInfinityStep;
    }

    void resetAcceptNanStep( const bool acceptNanStep )
    {
        acceptNanStep_ = acceptNanStep;
    }

protected:
    bool acceptInfinityStep_;

    bool acceptNanStep_;
};

enum MinimumIntegrationTimeStepHandling
{
    throw_exception_below_minimum,
    set_to_minimum_step_silently,
    set_to_minimum_step_single_warning,
    set_to_minimum_step_every_time_warning
};

template< typename TimeStepType >
class BasicIntegratorStepSizeValidator: public IntegratorStepSizeValidator< TimeStepType >
{
public:
    BasicIntegratorStepSizeValidator(
        const TimeStepType minimumStep,
        const TimeStepType maximumStep,
        const MinimumIntegrationTimeStepHandling minimumIntegrationTimeStepHandling = throw_exception_below_minimum,
        const bool acceptInfinityStep = false,
        const bool acceptNanStep = false ):
        IntegratorStepSizeValidator< TimeStepType >( acceptInfinityStep, acceptNanStep ),
        minimumStep_( minimumStep ), maximumStep_( maximumStep ),
        minimumIntegrationTimeStepHandling_( minimumIntegrationTimeStepHandling ){ }

    virtual ~BasicIntegratorStepSizeValidator( ){ }

    std::pair< TimeStepType, bool > validateStep(
        const std::pair< TimeStepType, bool > recommendedStep, const double currentStep )
    {
        bool acceptStep = recommendedStep.second;
        double newStepSize = TUDAT_NAN;

        if ( std::fabs( recommendedStep.first ) < std::fabs( minimumStep_ ) )
        {
            if( minimumIntegrationTimeStepHandling_ == throw_exception_below_minimum )
            {
                throw std::runtime_error( "Error in step-size control, minimum step size " + std::to_string( minimumStep_ ) +
                " is higher than required time step " + std::to_string( recommendedStep.first ) );
            }
            else
            {
                if( ( minimumIntegrationTimeStepHandling_ == set_to_minimum_step_every_time_warning ) ||
                    ( ( minimumIntegrationTimeStepHandling_ == set_to_minimum_step_single_warning ) &&!minimumStepWarningIsPrinted_ ) )
                {
                    std::cerr<<"Warning in step-size control, minimum step size " + std::to_string( minimumStep_ ) +
                                              " is higher than required time step " + std::to_string( recommendedStep.first ) + ", minimum step will be used."<<std::endl;
                }
                minimumStepWarningIsPrinted_ = true;
                newStepSize = currentStep / std::fabs( currentStep ) * std::fabs( maximumStep_ );

            }
        }
        else if( std::fabs( recommendedStep.first ) > std::fabs( maximumStep_ ) )
        {
            newStepSize = currentStep / std::fabs( currentStep ) * std::fabs( maximumStep_ );
        }
        else
        {
            newStepSize = recommendedStep.first;
        }

        if( newStepSize * recommendedStep.first < 0 )
        {
            throw std::runtime_error( "Error during step size control, step size flipped sign" );
        }

        if( std::isnan( newStepSize ) && !this->acceptNanStep_ )
        {
            throw std::runtime_error( "Error, recommended step is NaN" );
        }
        else if( std::isinf( newStepSize ) && !this->acceptInfinityStep_ )
        {
            throw std::runtime_error( "Error, recommended step is NaN" );
        }

        // Check if computed error in state is too large and reject step if true.
        return std::make_pair(  newStepSize, acceptStep );
    }

    void resetMinimumIntegrationTimeStepHandling( const MinimumIntegrationTimeStepHandling minimumIntegrationTimeStepHandling )
    {
        minimumIntegrationTimeStepHandling_ = minimumIntegrationTimeStepHandling;
    }

    void restartPropagation( )
    {
        minimumStepWarningIsPrinted_ = false;
    }

protected:

    const TimeStepType minimumStep_;

    const TimeStepType maximumStep_;

    MinimumIntegrationTimeStepHandling minimumIntegrationTimeStepHandling_;

    bool minimumStepWarningIsPrinted_;
};

template< typename TimeStepType, typename StateType = Eigen::VectorXd >
class IntegratorStepSizeController
{
public:
    IntegratorStepSizeController(
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        integratorOrder_( static_cast< double >( integratorOrder ) ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ),
        maximumFactorDecreaseForNextStepSize_( maximumFactorDecreaseForNextStepSize )
    { }

    virtual ~IntegratorStepSizeController( )
    { }

    virtual void initialize( const StateType& state ){ }

    virtual std::pair< TimeStepType, bool > computeNewStepSize(
        const StateType &firstStateEstimate,
        const StateType &secondStateEstimate,
        const TimeStepType &currentStep ) = 0;


protected:

    std::pair< TimeStepType, bool > computeTimeStepFromErrorEstimate(
        const TimeStepType& maximumErrorInState,
        const TimeStepType& currentStep )
    {
        // Compute the new step size. This is based off of the equation given in
        // (Montenbruck and Gill, 2005).
        const TimeStepType timeStepRatio = safetyFactorForNextStepSize_
                                           * std::pow( 1.0 / static_cast< double >( maximumErrorInState ),
                                                       1.0 / static_cast< double >( integratorOrder_ ) );

        bool tolerancesMet = maximumErrorInState <= 1.0;
        if ( timeStepRatio <= minimumFactorDecreaseForNextStepSize_ )
        {
            return std::make_pair( currentStep * minimumFactorDecreaseForNextStepSize_, tolerancesMet );
        }
        else if ( timeStepRatio >= maximumFactorDecreaseForNextStepSize_ )
        {
            return std::make_pair( currentStep * maximumFactorDecreaseForNextStepSize_, tolerancesMet );
        }
        else
        {
            return std::make_pair( currentStep * timeStepRatio, tolerancesMet );
        }
    }

    const double safetyFactorForNextStepSize_;

    const double integratorOrder_;

    const double minimumFactorDecreaseForNextStepSize_;

    const double maximumFactorDecreaseForNextStepSize_;

};

template< typename TimeStepType, typename StateType = Eigen::VectorXd >
class PerElementIntegratorStepSizeController: public IntegratorStepSizeController< TimeStepType, StateType >
{
public:
    PerElementIntegratorStepSizeController(
        const StateType relativeErrorTolerance,
        const StateType absoluteErrorTolerance,
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        IntegratorStepSizeController< TimeStepType, StateType >(
            safetyFactorForNextStepSize, integratorOrder, minimumFactorDecreaseForNextStepSize, maximumFactorDecreaseForNextStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ),
        absoluteErrorTolerance_( absoluteErrorTolerance ),
        tolerancesSet_( true ){ }

    PerElementIntegratorStepSizeController(
        const double relativeErrorTolerance,
        const double absoluteErrorTolerance,
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        IntegratorStepSizeController< TimeStepType, StateType >(
            safetyFactorForNextStepSize, integratorOrder, minimumFactorDecreaseForNextStepSize, maximumFactorDecreaseForNextStepSize ),
        scalarRelativeErrorTolerance_( relativeErrorTolerance ),
        scalarAbsoluteErrorTolerance_( absoluteErrorTolerance ),
        tolerancesSet_( false ){ }

    virtual ~PerElementIntegratorStepSizeController( ){ }

    void initialize( const StateType& state )
    {
        if( !tolerancesSet_ )
        {
            relativeErrorTolerance_ = StateType::Constant( state.rows( ), state.cols( ),
                                                          std::fabs( scalarRelativeErrorTolerance_ ) );
            absoluteErrorTolerance_ = StateType::Constant( state.rows( ), state.cols( ),
                                                           std::fabs( scalarAbsoluteErrorTolerance_ ) );
            tolerancesSet_ = true;
        }
        else
        {
            if( ( relativeErrorTolerance_.rows( ) != state.rows( ) ) || ( relativeErrorTolerance_.cols( ) != state.cols( ) ) )
            {
                throw std::runtime_error( "Error in step size controller, size of toleranes is incompatible with state" );
            }
        }
    }

    std::pair< TimeStepType, bool > computeNewStepSize(
        const StateType& firstStateEstimate,
        const StateType& secondStateEstimate,
        const TimeStepType& currentStep )
    {
        if( !tolerancesSet_ )
        {
            throw std::runtime_error( "Error in per-element step size control; tolerances not initialized" );
        }
        // Compute the truncation error based on the higher and lower order estimates.
        const StateType truncationError_ =
            ( firstStateEstimate - secondStateEstimate ).array( ).abs( );

        // Compute error tolerance based on relative and absolute error tolerances.
        const StateType errorTolerance_ =
            ( firstStateEstimate.array( ).abs( ) *
              relativeErrorTolerance_.array( ) ).matrix( )
            + absoluteErrorTolerance_;

        // Compute relative truncation error. This will indicate if the current step satisfies the
        // required tolerances.
        const StateType relativeTruncationError_ = truncationError_.array( ) /
                                                   errorTolerance_.array( );

        // Compute the maximum error based on the largest coefficient in the relative truncation error
        // matrix.
        const typename StateType::Scalar maximumErrorInState_
            = relativeTruncationError_.array( ).abs( ).maxCoeff( );
        
        return this->computeTimeStepFromErrorEstimate( maximumErrorInState_, currentStep );

    }

protected:

    StateType relativeErrorTolerance_;

    StateType absoluteErrorTolerance_;

    const double scalarRelativeErrorTolerance_ = TUDAT_NAN;

    const double scalarAbsoluteErrorTolerance_ = TUDAT_NAN;

    bool tolerancesSet_;

};


template< typename TimeStepType, typename StateType = Eigen::VectorXd >
class PerBlockIntegratorStepSizeController: public IntegratorStepSizeController< TimeStepType, StateType >
{
public:


    PerBlockIntegratorStepSizeController(
        const std::function< std::vector< std::tuple< int, int, int, int > >( const int, const int ) > blocksToCheckFunction,
        const StateType relativeErrorTolerance,
        const StateType absoluteErrorTolerance,
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        IntegratorStepSizeController< TimeStepType, StateType >(
            safetyFactorForNextStepSize, integratorOrder, minimumFactorDecreaseForNextStepSize, maximumFactorDecreaseForNextStepSize ),
        blocksToCheckFunction_( blocksToCheckFunction ),
        relativeErrorTolerance_( relativeErrorTolerance ),
        absoluteErrorTolerance_( absoluteErrorTolerance ),
        tolerancesSet_( true )
    {
    }

    PerBlockIntegratorStepSizeController(
        const std::function< std::vector< std::tuple< int, int, int, int > >( const int, const int ) > blocksToCheckFunction,
        const double relativeErrorTolerance,
        const double absoluteErrorTolerance,
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        IntegratorStepSizeController< TimeStepType, StateType >(
            safetyFactorForNextStepSize, integratorOrder, minimumFactorDecreaseForNextStepSize, maximumFactorDecreaseForNextStepSize ),
        blocksToCheckFunction_( blocksToCheckFunction ),
        scalarRelativeErrorTolerance_( relativeErrorTolerance ),
        scalarAbsoluteErrorTolerance_( absoluteErrorTolerance ),
        tolerancesSet_( false )
    {
    }

    void initialize( const StateType& state )
    {
        if( blocksToCheck_.size( ) == 0 )
        {
            if( blocksToCheckFunction_ == nullptr )
            {
                throw std::runtime_error( "Error when setting per-segment step-size control, no blocks are provided" );
            }
            blocksToCheck_ = blocksToCheckFunction_( state.rows( ), state.cols( ) );
            if( blocksToCheck_.size( ) == 0 )
            {
                throw std::runtime_error( "Error when setting per-segment step-size control, no blocks are provided by blocks-to-check function" );
            }
        }

        relativeTruncationError_.resize( blocksToCheck_.size( ), 1 );

        for( unsigned int i = 0; i < blocksToCheck_.size( ); i++ )
        {
            int maximumRow = std::get< 0 >( blocksToCheck_.at( i ) ) + std::get< 2 >( blocksToCheck_.at( i ) );
            if( maximumRow > state.rows( ) )
            {
                throw std::runtime_error( "Error when setting per-segment step-size control, block to check is out of bounds. Number of rows is " +
                std::to_string( state.rows( ) ) + ", but control is requested on row segment" +
                std::to_string( std::get< 0 >( blocksToCheck_.at( i ) )  ) + ", " +
                std::to_string( std::get< 2 >( blocksToCheck_.at( i ) ) ) );
            }

            int maximumColumn = std::get< 1 >( blocksToCheck_.at( i ) ) + std::get< 3 >( blocksToCheck_.at( i ) );
            if( maximumColumn > state.cols( ) )
            {
                throw std::runtime_error( "Error when setting per-segment step-size control, block to check is out of bounds. Number of columns is " +
                                          std::to_string( state.cols( ) ) + ", but control is requested on column segment " +
                                          std::to_string( std::get< 1 >( blocksToCheck_.at( i ) )  ) + ", " +
                                          std::to_string( std::get< 3 >( blocksToCheck_.at( i ) ) ) );
            }
        }

        if( !tolerancesSet_ )
        {
            relativeErrorTolerance_ = StateType::Constant( blocksToCheck_.size( ), 1,
                                                           std::fabs( scalarRelativeErrorTolerance_ ) );
            absoluteErrorTolerance_ = StateType::Constant( blocksToCheck_.size( ), 1,
                                                           std::fabs( scalarAbsoluteErrorTolerance_ ) );
            tolerancesSet_ = true;
        }
        else
        {
            if( ( relativeErrorTolerance_.rows( ) != static_cast< int >( blocksToCheck_.size( ) ) ) || ( relativeErrorTolerance_.cols( ) != 1 ) )
            {
                throw std::runtime_error( "Error in per-segment step size controller, size of tolerances is incompatible with state blocks" );
            }
        }
    }


    virtual ~PerBlockIntegratorStepSizeController( ){ }

    std::pair< TimeStepType, bool > computeNewStepSize(
        const StateType& firstStateEstimate,
        const StateType& secondStateEstimate,
        const TimeStepType& currentStep )
    {
        if( blocksToCheck_.size( ) == 0 )
        {
            throw std::runtime_error( "Error when using per-segment step-size control, no blocks are provided" );
        }

        // Compute the truncation error based on the higher and lower order estimates.
        const StateType truncationError_ =
            ( firstStateEstimate - secondStateEstimate ).array( ).abs( );

        for( unsigned int i = 0; i < blocksToCheck_.size( ); i++ )
        {
            relativeTruncationError_( i ) = truncationError_.block(
                std::get< 0 >( blocksToCheck_.at( i ) ),
                std::get< 1 >( blocksToCheck_.at( i ) ),
                std::get< 2 >( blocksToCheck_.at( i ) ),
                std::get< 3 >( blocksToCheck_.at( i ) ) ).norm( ) /
                ( firstStateEstimate.block(
                    std::get< 0 >( blocksToCheck_.at( i ) ),
                    std::get< 1 >( blocksToCheck_.at( i ) ),
                    std::get< 2 >( blocksToCheck_.at( i ) ),
                    std::get< 3 >( blocksToCheck_.at( i ) ) ).norm( ) *
                  relativeErrorTolerance_( i ) + absoluteErrorTolerance_( i ) );
        }

        // Compute the maximum error based on the largest coefficient in the relative truncation error
        // matrix.
        const typename StateType::Scalar maximumErrorInState_
            = relativeTruncationError_.array( ).abs( ).maxCoeff( );

        return this->computeTimeStepFromErrorEstimate( maximumErrorInState_, currentStep );

    }

protected:

    std::vector< std::tuple< int, int, int, int > > blocksToCheck_;

    std::function< std::vector< std::tuple< int, int, int, int > >( const int, const int ) > blocksToCheckFunction_;

    StateType relativeErrorTolerance_;

    StateType absoluteErrorTolerance_;

    const double scalarRelativeErrorTolerance_ = TUDAT_NAN;

    const double scalarAbsoluteErrorTolerance_ = TUDAT_NAN;

    bool tolerancesSet_;

    StateType relativeTruncationError_;
};


template< typename TimeStepType, typename StateType = Eigen::VectorXd >
class CustomIntegratorStepSizeController: public IntegratorStepSizeController< TimeStepType, StateType >
{
public:
    CustomIntegratorStepSizeController(
        const std::function< TimeStepType( const StateType&, const StateType& ) > customErrorFunction,
        const double safetyFactorForNextStepSize,
        const int integratorOrder,
        const double minimumFactorDecreaseForNextStepSize,
        const double maximumFactorDecreaseForNextStepSize ):
        IntegratorStepSizeController< TimeStepType, StateType >(
            safetyFactorForNextStepSize, integratorOrder, minimumFactorDecreaseForNextStepSize, maximumFactorDecreaseForNextStepSize ),
        customErrorFunction_( customErrorFunction )
    { }

    virtual ~CustomIntegratorStepSizeController( ){ }

    std::pair< TimeStepType, bool > computeNewStepSize(
        const StateType& firstStateEstimate,
        const StateType& secondStateEstimate,
        const TimeStepType& currentStep )
    {
          return computeTimeStepFromErrorEstimate( customErrorFunction_( firstStateEstimate, secondStateEstimate ), currentStep );
    }

protected:

    std::function< TimeStepType( const StateType&, const StateType& ) > customErrorFunction_;
};

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_STEPSIZE_CONTROLLER_H
