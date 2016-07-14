#ifndef THRUSTMAGNITUDEWRAPPER_H
#define THRUSTMAGNITUDEWRAPPER_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace propulsion
{


class ThrustMagnitudeWrapper
{
public:

    ThrustMagnitudeWrapper( ){ }

    virtual ~ThrustMagnitudeWrapper( ){ }

    virtual void update( const double time ) = 0;

    virtual double getCurrentThrust( ) = 0;

    virtual double getCurrentMassRate( ) = 0;

};


class CustomThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    CustomThrustMagnitudeWrapper(
            const boost::function< double( const double ) > thrustMagnitudeFunction,
            const boost::function< double( const double ) > specificImpulseFunction,
            const boost::function< bool( const double ) > isEngineOnFunction = boost::lambda::constant( true ) ):
    thrustMagnitudeFunction_( thrustMagnitudeFunction ),
    specificImpulseFunction_( specificImpulseFunction ),
    isEngineOnFunction_( isEngineOnFunction ),
    currentThrustMagnitude_( TUDAT_NAN ),
    currentSpecificImpulse_( TUDAT_NAN ){ }

    void update( const double time )
    {
        if( isEngineOnFunction_( time ) )
        {
            currentThrustMagnitude_ = thrustMagnitudeFunction_( time );
            currentSpecificImpulse_ = specificImpulseFunction_( time );
        }
        else
        {
            currentThrustMagnitude_ = 0.0;
            currentSpecificImpulse_ = TUDAT_NAN;
        }
    }

    double getCurrentThrust( )
    {
        return currentThrustMagnitude_;
    }

    double getCurrentMassRate( )
    {
        return propulsion::computePropellantMassRateFromSpecificImpulse(
                    currentThrustMagnitude_, currentSpecificImpulse_ );
    }

private:
    boost::function< double( const double ) > thrustMagnitudeFunction_;

    boost::function< double( const double ) > specificImpulseFunction_;

    boost::function< bool( const double ) > isEngineOnFunction_;

    double currentThrustMagnitude_;

    double currentSpecificImpulse_;

};

class ThrustMagnitudeFromEngineWrapper: public ThrustMagnitudeWrapper
{
public:

    ThrustMagnitudeFromEngineWrapper(
            const boost::shared_ptr< system_models::EngineModel > engineModel ):
        currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    {
        engineModels_.push_back( engineModel );
    }

    ThrustMagnitudeFromEngineWrapper(
            const std::vector< boost::shared_ptr< system_models::EngineModel > > engineModels ):
        engineModels_( engineModels ), currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    { }

    ~ThrustMagnitudeFromEngineWrapper( ){ }

    void update( const double time )
    {
        currentThrust_ = 0.0;
        currentMassRate_ = 0.0;

        for( unsigned int i = 0; i < engineModels_.size( ); i++ )
        {
            engineModels_.at( i )->updateEngineModel( time );
        }

        for( unsigned int i = 0; i < engineModels_.size( ); i++ )
        {
            currentThrust_ += engineModels_.at( i )->getCurrentThrust( );
            currentMassRate_ += engineModels_.at( i )->getCurrentMassRate( );

        }
    }

    double getCurrentThrust( )
    {
        return currentThrust_;
    }

    double getCurrentMassRate( )
    {
        return currentMassRate_;
    }

protected:

    std::vector< boost::shared_ptr< system_models::EngineModel > > engineModels_;

    double currentThrust_;

    double currentMassRate_;
};

}

}

#endif // THRUSTMAGNITUDEWRAPPER_H
