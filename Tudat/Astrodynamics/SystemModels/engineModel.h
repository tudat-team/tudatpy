#ifndef ENGINEMODEL_H
#define ENGINEMODEL_H

#include <boost/function.hpp>

namespace tudat
{

namespace system_models
{

class EngineModel
{
public:

    double getCurrentThrust( )
    {
        return currentThrust_;
    }

    virtual ~EngineModel( ){ }

    virtual void updateEngineModel( const double currentTime ) = 0;

    virtual double getCurrentMassRate( ) = 0;

protected:

    double currentThrust_;

};

class DirectEngineModel: public EngineModel
{
public:

    DirectEngineModel(
            const double specificImpulse,
            const boost::function< double( ) > massFlowFunction ):
        specificImpulse_( specificImpulse ),
        massFlowFunction_( massFlowFunction )
    { }

    void updateEngineModel( const double currentTime )
    {
        currentThrust_ = specificImpulse_ * massFlowFunction_( );
    }

    double getCurrentMassRate( )
    {
        return massFlowFunction_( );
    }


protected:

    double specificImpulse_;

    boost::function< double( ) > massFlowFunction_;
};

}

}

#endif // ENGINEMODEL_H
