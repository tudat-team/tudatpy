/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <memory>
#include <functional>
#include <boost/lambda/lambda.hpp>
#include <iostream>

#include "tudat/math/interpolators/interpolator.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"
#include "tudat/astro/propulsion/thrustFunctions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"

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
    virtual double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN ) = 0;

    virtual double getCurrentThrustAccelerationMagnitude( const double currentMass = TUDAT_NAN ) = 0;


    //! Pure virtual function to get the current mass rate.
    /*!
     * Pure virtual function to get the current mass rate.
     * \return Current mass rate.
     */
    virtual double getCurrentMassRate( const double currentMass = TUDAT_NAN ) = 0;

    virtual double getCurrentSpecificImpulse( ) = 0;

    virtual bool modelIsForceBased( ) = 0;

    //! Function to reset the current time of the thrust model
    /*!
     *  Function to reset the current time of the thrust model. Function is typically used to reset the time to NaN,
     *  signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
        resetDerivedClassCurrentTime( );
    }

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities. This function can be redefined in
     *  derived class
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( )
    {

    }


protected:

    //! Current time for model.
    double currentTime_;
};

class ConstantThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    ConstantThrustMagnitudeWrapper(
            const double thrustMagnitude,
            const double specificImpulse ):
        thrustMagnitude_( thrustMagnitude ),
        specificImpulse_( specificImpulse )
    {
        massRate_ = computePropellantMassRateFromSpecificImpulse(
                    thrustMagnitude_, specificImpulse_ );
    }

    //! Destructor.
    ~ConstantThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            currentTime_ = time;
        }
    }

    double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN )
    {
        return thrustMagnitude_;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass )
    {
        return thrustMagnitude_ / currentMass;
    }


    void resetConstantThrustForceMagnitude( const double thrustMagnitude )
    {
        thrustMagnitude_ = thrustMagnitude;
        massRate_ = computePropellantMassRateFromSpecificImpulse(
                    thrustMagnitude_, specificImpulse_ );
    }


    double getCurrentMassRate( const double currentMass = TUDAT_NAN )
    {
        return massRate_;
    }

    double getCurrentSpecificImpulse( )
    {
        return specificImpulse_;
    }

    void resetConstantSpecificImpulse( const double specificImpulse )
    {
        specificImpulse_ = specificImpulse;
        massRate_ = computePropellantMassRateFromSpecificImpulse(
                    thrustMagnitude_, specificImpulse_ );
    }


    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }


    bool modelIsForceBased( )
    {
        return true;
    }

private:

    double thrustMagnitude_;

    double specificImpulse_;

    double massRate_;

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
            const std::function< double( const double ) > thrustMagnitudeFunction,
            const std::function< double( const double ) > specificImpulseFunction ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        isSpecificImpulseConstant_( false ){ }

    CustomThrustMagnitudeWrapper(
            const std::function< double( const double ) > thrustMagnitudeFunction,
            const double specificImpulse ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( [=](const double){return specificImpulse;} ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        isSpecificImpulseConstant_( true ){ }

    //! Destructor.
    ~CustomThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time );

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN )
    {
        return currentThrustMagnitude_;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass )
    {
        return currentThrustMagnitude_ / currentMass;
    }

    bool modelIsForceBased( )
    {
        return true;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( const double currentMass = TUDAT_NAN )
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

    double getCurrentSpecificImpulse( )
    {
        return currentSpecificImpulse_;
    }

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( )
    {
        thrustMagnitudeFunction_( TUDAT_NAN );
        specificImpulseFunction_( TUDAT_NAN );
        currentThrustMagnitude_ = TUDAT_NAN;
        currentSpecificImpulse_ = TUDAT_NAN;
        currentTime_ = TUDAT_NAN;
    }

    bool isSpecificImpulseConstant( )
    {
        return isSpecificImpulseConstant_;
    }

private:

    //! Function returning thrust as a function of time..
    std::function< double( const double ) > thrustMagnitudeFunction_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    bool isSpecificImpulseConstant_;

};


class CustomThrustAccelerationMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    CustomThrustAccelerationMagnitudeWrapper(
            const std::function< double( const double ) > thrustAccelerationMagnitudeFunction,
            const std::function< double( const double ) > specificImpulseFunction ):
        thrustAccelerationMagnitudeFunction_( thrustAccelerationMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        currentThrustAccelerationMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        isSpecificImpulseConstant_( false ){ }

    CustomThrustAccelerationMagnitudeWrapper(
            const std::function< double( const double ) > thrustAccelerationMagnitudeFunction,
            const double specificImpulse ):
        thrustAccelerationMagnitudeFunction_( thrustAccelerationMagnitudeFunction ),
        specificImpulseFunction_( [=](const double){return specificImpulse;} ),
        currentThrustAccelerationMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        isSpecificImpulseConstant_( true ){ }

    ~CustomThrustAccelerationMagnitudeWrapper( ){ }

    void update( const double time );

    double getCurrentThrustForceMagnitude( const double currentMass )
    {
        return currentThrustAccelerationMagnitude_ * currentMass;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass = TUDAT_NAN )
    {
        return currentThrustAccelerationMagnitude_;
    }

    bool modelIsForceBased( )
    {
        return false;
    }

    double getCurrentMassRate( const double currentMass )
    {
        if( currentThrustAccelerationMagnitude_ != 0.0 )
        {
            return propulsion::computePropellantMassRateFromSpecificImpulse(
                        getCurrentThrustForceMagnitude( currentMass ), currentSpecificImpulse_ );
        }
        else
        {
            return 0.0;
        }
    }

    double getCurrentSpecificImpulse( )
    {
        return currentSpecificImpulse_;
    }

    virtual void resetDerivedClassCurrentTime( )
    {
        thrustAccelerationMagnitudeFunction_( TUDAT_NAN );
        specificImpulseFunction_( TUDAT_NAN );
        currentThrustAccelerationMagnitude_ = TUDAT_NAN;
        currentSpecificImpulse_ = TUDAT_NAN;
        currentTime_ = TUDAT_NAN;
    }

    bool isSpecificImpulseConstant( )
    {
        return isSpecificImpulseConstant_;
    }

private:

    std::function< double( const double ) > thrustAccelerationMagnitudeFunction_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    double currentThrustAccelerationMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    bool isSpecificImpulseConstant_;

};


//! Class for bang-bang thrust magnitude from MEE co-states (optimal control theory).
/*!
 *  Class for bang-bang thrust magnitude from MEE co-states (optimal control theory).
 */
class MeeCostatesBangBangThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
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
    MeeCostatesBangBangThrustMagnitudeWrapper(
            const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
            const std::function< Eigen::Vector6d( ) > centralBodyStateFunction,
            const std::function< double( ) > centralBodyGravitationalParameterFunction,
            const std::function< Eigen::VectorXd( const double ) > costateFunction,
            const double thrustMagnitude,
            const std::function< double( const double ) > specificImpulseFunction,
            const std::function< double( ) > thrustingBodyMassFunction,
            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
        thrustingBodyStateFunction_( thrustingBodyStateFunction ),
        centralBodyStateFunction_( centralBodyStateFunction ),
        centralBodyGravitationalParameterFunction_( centralBodyGravitationalParameterFunction ),
        costateFunction_( costateFunction ),
        maximumThrustMagnitude_( thrustMagnitude ),
        specificImpulseFunction_( specificImpulseFunction ),
        thrustingBodyMassFunction_( thrustingBodyMassFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        customThrustResetFunction_( customThrustResetFunction ){ }

    //! Destructor.
    ~MeeCostatesBangBangThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time );

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN  )
    {
        return currentThrustMagnitude_;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass )
    {
        return currentThrustMagnitude_ / currentMass;
    }

    bool modelIsForceBased( )
    {
        return true;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( const double currentMass = TUDAT_NAN )
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

    double getCurrentSpecificImpulse( )
    {
        return currentSpecificImpulse_;
    }

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( )
    {
        if( !( customThrustResetFunction_ == nullptr ) )
        {
            customThrustResetFunction_( TUDAT_NAN );
        }
    }

private:

    //!  Function returning the state of the body under thrust as a function of time.
    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction_;

    //!  Function returning the state of the central body as a function of time.
    std::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //!  Function returning the gravitational parameter of the central body as a function of time.
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! General function which gives the costates vector as a function of time.
    std::function< Eigen::VectorXd( const double ) > costateFunction_;

    //! Thrust magnitude when the engine is on.
    double maximumThrustMagnitude_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    //! Function returning the mass of the body under thrust as a function of time.
    std::function< double( ) > thrustingBodyMassFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    std::function< void( const double ) > customThrustResetFunction_;

};

//! Variables on which parameterized thrust can depend.
enum ThrustIndependentVariables
{
    mach_number_dependent_thrust,
    altitude_dependent_thrust,
    density_dependent_thrust,
    dynamic_pressure_dependent_thrust,
    pressure_dependent_thrust,
    guidance_input_dependent_thrust,
    throttle_dependent_thrust
};

//! Class to compute the engine thrust and specific impulse as a parameterized function of any number of indepedent
//! variables.
/*!
 *  Class to compute the engine thrust and specific impulse as a parameterized function of any number of indepedent
 *  variables.  The physical meaning of the variables must be defined here, selecting from the options in the
 *  ThrustIndependentVariables enum, and they are automatically retrieved from the relevant environment models during the
 *  propagation.
 *  Note that any number of user-specific functions may be included, as a  guidance_input_dependent_thrust type or
 *  throttle_dependent_thrust. Note that a throttle_dependent_thrust may not be used as one of the independent variables
 *  of the thrust magnitude for this class. This setting is parsed through the ParameterizedThrustMagnitudeSettings
 *  class, which is the class that is typically used to create this class.
 */
class ParameterizedThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor for parameterized thrust and specific impulse.
    /*!
     * Constructor, defines the functions for thrust and specific impulse, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeFunction Function returning the current thrust as a function of the independent variables.
     * \param specificImpulseFunction  Function returning the current specific impulse as a function of the
     * independent variables.
     * \param thrustInputVariableFunctions List of functions returning input variables for the thrust
     * The order of the functions in this vector is passed to the thrustMagnitudeFunction in the same order as
     * entries of this vector.
     * \param specificImpulseInputVariableFunctions List of functions returning input variables for the
     * specific impulse. The order of the functions in this vector is passed to the specificImpulseFunction in the same order
     * as entries of this vector.
     * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the thrustMagnitudeFunction function.
     * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
     * input to the specificImpulseDependentVariables function.
     * \param inputUpdateFunction Function that is called to update the user-defined guidance to the current time
     * (empty by default).
     */
    ParameterizedThrustMagnitudeWrapper(
            const std::function< double( const std::vector< double >& ) > thrustMagnitudeFunction,
            const std::function< double( const std::vector< double >& ) > specificImpulseFunction,
            const std::vector< std::function< double( ) > > thrustInputVariableFunctions,
            const std::vector< std::function< double( ) > > specificImpulseInputVariableFunctions,
            const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
            const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
            const std::function< void( const double) > inputUpdateFunction =
            std::function< void( const double) >( ) ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        thrustInputVariableFunctions_( thrustInputVariableFunctions ),
        specificImpulseInputVariableFunctions_( specificImpulseInputVariableFunctions ),
        thrustIndependentVariables_( thrustIndependentVariables ),
        specificImpulseDependentVariables_( specificImpulseDependentVariables ),
        inputUpdateFunction_( inputUpdateFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN )
    {
        if( thrustInputVariableFunctions_.size( ) != thrustIndependentVariables_.size( ) )
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
    void update( const double time );

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime(  )
    {
        if( !( inputUpdateFunction_ == nullptr ) )
        {
            inputUpdateFunction_( TUDAT_NAN );
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN )
    {
        return currentThrustMagnitude_;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass )
    {
        return currentThrustMagnitude_ / currentMass;
    }

    bool modelIsForceBased( )
    {
        return true;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( const double currentMass = TUDAT_NAN  )
    {
        return propulsion::computePropellantMassRateFromSpecificImpulse(
                    currentThrustMagnitude_, currentSpecificImpulse_ );
    }

    double getCurrentSpecificImpulse( )
    {
        return currentSpecificImpulse_;
    }

private:

    //! Function returning the current thrust as a function of its independent variables
    std::function< double( const std::vector< double >& ) > thrustMagnitudeFunction_;

    //! Function returning the current specific impulse as a function of its independent variables
    std::function< double( const std::vector< double >& ) > specificImpulseFunction_;

    //!  List of functions returning input variables for thrustMagnitudeFunction_.
    std::vector< std::function< double( ) > > thrustInputVariableFunctions_;

    //!  List of functions returning input variables for specificImpulseFunction_.
    std::vector< std::function< double( ) > > specificImpulseInputVariableFunctions_;

    //! List of identifiers for the physical meaning of each of the entries of the input to the thrustMagnitudeFunction
    //! function.
    std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables_;

    //! List of identifiers for the physical meaning of each of the entries of the input to the
    //! specificImpulseDependentVariables function.
    std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables_;

    //! List of current input data to thrust function
    std::vector< double > currentThrustInputVariables_;

    //! List of current input data to specific impulse function.
    std::vector< double > currentSpecificImpulseInputVariables_;


    //! Function that is called to update the user-defined guidance to the current time
    const std::function< void( const double) > inputUpdateFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

};

class CustomThrustVectorWrapper: public ThrustMagnitudeWrapper, public ephemerides::InertialBodyFixedDirectionCalculator
{
    CustomThrustVectorWrapper(
            const std::function< Eigen::Vector3d( const double ) > thrustVectorFunction,
            const std::function< double( const double ) > specificImpulseFunction ):
        ThrustMagnitudeWrapper( ), ephemerides::InertialBodyFixedDirectionCalculator( ),
        thrustVectorFunction_( thrustVectorFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        isSpecificImpulseConstant_( false ){ }

    void update( const double time );

    double getCurrentThrustForceMagnitude( const double currentMass = TUDAT_NAN )
    {
        return currentThrustMagnitude_;
    }

    double getCurrentThrustAccelerationMagnitude( const double currentMass )
    {
        return currentThrustMagnitude_ / currentMass;
    }


    Eigen::Vector3d getDirection( const double time )
    {
        if( time != currentTime_ )
        {
            update( time );
        }
        return currentThrustVector_.normalized( );
    }


    bool modelIsForceBased( )
    {
        return true;
    }

    double getCurrentMassRate( const double currentMass = TUDAT_NAN )
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

    double getCurrentSpecificImpulse( )
    {
        return currentSpecificImpulse_;
    }

    virtual void resetDerivedClassCurrentTime( )
    {
        thrustVectorFunction_( NAN );
        specificImpulseFunction_( TUDAT_NAN );
        currentThrustMagnitude_ = TUDAT_NAN;
        currentSpecificImpulse_ = TUDAT_NAN;
        currentTime_ = TUDAT_NAN;
    }

    bool isSpecificImpulseConstant( )
    {
        return isSpecificImpulseConstant_;
    }

private:

    //! Function returning thrust as a function of time..
    std::function< Eigen::Vector3d( const double ) > thrustVectorFunction_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    Eigen::Vector3d currentThrustVector_;

    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    bool isSpecificImpulseConstant_;
};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_THRUSTMAGNITUDEWRAPPER_H
