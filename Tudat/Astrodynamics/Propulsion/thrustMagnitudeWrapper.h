/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THRUSTMAGNITUDEWRAPPER_H
#define TUDAT_THRUSTMAGNITUDEWRAPPER_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"

namespace tudat
{

namespace propulsion
{

//! Base class for determining the current thrust magnitude from method defined in derived class.
class ThrustMagnitudeWrapper
{
public:

    //! Constructor
    ThrustMagnitudeWrapper( ):
        currentTime_( TUDAT_NAN ){ }

    //! Destructor.
    virtual ~ThrustMagnitudeWrapper( ){ }

    //! Pure virtual function to update the thrust magnitude to the current time.
    /*!
     * Pure virtual function to update the thrust magnitude to the current time. Method of computation is to be
     * implemented in the derived class.
     * \param time Time at which the magitude is to be computed.
     */
    virtual void update( const double time ) = 0;

    //! Pure virtual function to get the current thrust magnitude.
    /*!
     * Pure virtual function to get the current thrust magnitude.
     * \return Current thrust magnitude.
     */
    virtual double getCurrentThrustMagnitude( ) = 0;

    //! Pure virtual function to get the current mass rate.
    /*!
     * Pure virtual function to get the current mass rate.
     * \return Current mass rate.
     */
    virtual double getCurrentMassRate( ) = 0;

    //! Function to reset the current time of the thrust model
    /*!
     *  Function to reset the current time of the thrust model. Function is typically used to reset the time to NaN,
     *  signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        resetDerivedClassCurrentTime( currentTime );
    }

    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {

    }


protected:

    //! Current time for model.
    double currentTime_;
};

//! Class for custom computations of thrust magnitude and mass rate.
/*!
 *  Class for custom computations of thrust magnitude and mass rate. Both the magnitude and specific impulse
 *  are parameterized as a function of time and provided as input, providing a flexible interface for user-defined
 *  thrust magnitude profiles.
 */
class CustomThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning whether the function is on (returns true if so) at a given time.
     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
     * being started (empty by default)
     */
    CustomThrustMagnitudeWrapper(
            const boost::function< double( const double ) > thrustMagnitudeFunction,
            const boost::function< double( const double ) > specificImpulseFunction,
            const boost::function< bool( const double ) > isEngineOnFunction = boost::lambda::constant( true ) ,
            const boost::function< void( const double ) > customThrustResetFunction = boost::function< void( const double ) >( ) ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        isEngineOnFunction_( isEngineOnFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        customThrustResetFunction_( customThrustResetFunction ){ }

    //! Destructor.
    ~CustomThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            // If engine is one, update engine.
            if( isEngineOnFunction_( time ) )
            {
                currentThrustMagnitude_ = thrustMagnitudeFunction_( time );
                currentSpecificImpulse_ = specificImpulseFunction_( time );
            }
            // If engine is off, set thrust to zero.
            else
            {
                currentThrustMagnitude_ = 0.0;
                currentSpecificImpulse_ = TUDAT_NAN;
            }
            currentTime_ = time;
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( )
    {
        if( currentThrustMagnitude_ != 0.0 )
        {
            return propulsion::computePropellantMassRateFromSpecificImpulse(
                        currentThrustMagnitude_, currentSpecificImpulse_ );
        }
        else
        {
            return 0.0;
        }
    }

    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( customThrustResetFunction_.empty( ) ) )
        {
            customThrustResetFunction_( currentTime );\
        }
    }

private:

    //! Function returning thrust as a function of time..
    boost::function< double( const double ) > thrustMagnitudeFunction_;

    //! Function returning specific impulse as a function of time.
    boost::function< double( const double ) > specificImpulseFunction_;

    //! Function returning whether the function is on (returns true if so) at a given time.
    boost::function< bool( const double ) > isEngineOnFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    boost::function< void( const double ) > customThrustResetFunction_;

};


//! Class to compute the engine thrust and mass rate from EngineModel object(s).
class ThrustMagnitudeFromEngineWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor for single engine used for thrust.
    /*!
     * Constructor for single engine used for thrust
     * \param engineModel Engine model used to compute thrust/mass rate.
     */
    ThrustMagnitudeFromEngineWrapper(
            const boost::shared_ptr< system_models::EngineModel > engineModel ):
        currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    {
        engineModels_.push_back( engineModel );
    }

    //! Constructor for multiple engines used for thrust.
    /*!
     * Constructor for multiple engines used for thrust
     * \param engineModels List of engine models used to compute thrust/mass rate.
     */
    ThrustMagnitudeFromEngineWrapper(
            const std::vector< boost::shared_ptr< system_models::EngineModel > > engineModels ):
        engineModels_( engineModels ), currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    { }

    //! Destructor.
    ~ThrustMagnitudeFromEngineWrapper( ){ }


    //! Function to update the thrust magnitude to the current time from the engine models.
    /*!
     *  Function to update the thrust magnitude to the current time from the engine models
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ = time ) )
        {
            currentThrust_ = 0.0;
            currentMassRate_ = 0.0;

            // Update engine models
            for( unsigned int i = 0; i < engineModels_.size( ); i++ )
            {
                engineModels_.at( i )->updateEngineModel( time );
            }

            // Add thrust and mass rate.
            for( unsigned int i = 0; i < engineModels_.size( ); i++ )
            {
                currentThrust_ += engineModels_.at( i )->getCurrentThrust( );
                currentMassRate_ += engineModels_.at( i )->getCurrentMassRate( );

            }
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrust_;
    }

    //! Function to return the current mass rate
    /*!
     * Function to return the current mass rate, as computed by last call to update member function.
     * \return Current mass rate
     */
    double getCurrentMassRate( )
    {
        return currentMassRate_;
    }

protected:

    //! List of engine models used to compute thrust/mass rate.
    std::vector< boost::shared_ptr< system_models::EngineModel > > engineModels_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrust_;

    //! Current mass rate, as computed by last call to update member function.
    double currentMassRate_;
};

//! Variables on which parameterized thrust can depend.
enum ThrustDependentVariables
{
    time_dependent_thrust,
    mach_number_dependent_thrust,
    altitude_dependent_thrust,
    density_dependent_thrust,
    dynamic_pressure_dependent_thrust,
    pressure_dependent_thrust,
    guidance_input_dependent_thrust,
    maximum_thrust_multiplier
};

//! Class to compute the engine thrust and specific impulse as a parameterized function of any number of indepedent
//! variables.
class ParameterizedThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning whether the function is on (returns true if so) at a given time.
     */
    ParameterizedThrustMagnitudeWrapper(
            const boost::function< double( const std::vector< double >& ) > thrustMagnitudeFunction,
            const boost::function< double( const std::vector< double >& ) > specificImpulseFunction,
            const std::vector< boost::function< double( const double ) > > thrustInputVariableFunctions,
            const std::vector< boost::function< double( const double ) > > specificImpulseInputVariableFunctions,
            const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
            const std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        thrustInputVariableFunctions_( thrustInputVariableFunctions ),
        specificImpulseInputVariableFunctions_( specificImpulseInputVariableFunctions ),
        thrustDependentVariables_( thrustDependentVariables ),
        specificImpulseDependentVariables_( specificImpulseDependentVariables ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN )
    {
        if( thrustInputVariableFunctions_.size( ) != thrustDependentVariables_.size( ) )
        {
            throw std::runtime_error( "Error in parameterized thrust, inconsistent number of user-defined input variables for thrust" );
        }

        if( specificImpulseInputVariableFunctions_.size( ) != specificImpulseDependentVariables_.size( ) )
        {
            throw std::runtime_error( "Error in parameterized thrust, inconsistent number of user-defined input variables for Isp" );
        }

        currentThrustInputVariables_.resize( thrustInputVariableFunctions.size( ) );
        currentSpecificImpulseInputVariables_.resize( specificImpulseInputVariableFunctions.size( ) );
    }

    //! Destructor.
    ~ParameterizedThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            for( unsigned int i = 0; i < thrustInputVariableFunctions_.size( ); i++ )
            {
                currentThrustInputVariables_[ i ] = thrustInputVariableFunctions_.at( i )( time );
            }

            currentThrustMagnitude_ = thrustMagnitudeFunction_( currentThrustInputVariables_ );

            for( unsigned int i = 0; i < specificImpulseInputVariableFunctions_.size( ); i++ )
            {
                currentSpecificImpulseInputVariables_[ i ] = specificImpulseInputVariableFunctions_.at( i )( time );
            }

            currentSpecificImpulse_ = specificImpulseFunction_( currentSpecificImpulseInputVariables_ );
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( )
    {
        return propulsion::computePropellantMassRateFromSpecificImpulse(
                    currentThrustMagnitude_, currentSpecificImpulse_ );
    }

private:


    boost::function< double( const std::vector< double >& ) > thrustMagnitudeFunction_;

    boost::function< double( const std::vector< double >& ) > specificImpulseFunction_;

    std::vector< boost::function< double( const double ) > > thrustInputVariableFunctions_;

    std::vector< boost::function< double( const double ) > > specificImpulseInputVariableFunctions_;

    std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables_;

    std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables_;

    std::vector< double > currentThrustInputVariables_;

    std::vector< double > currentSpecificImpulseInputVariables_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_THRUSTMAGNITUDEWRAPPER_H
