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

    void updateEngineModel( const double currentTime ) = 0;

protected:

    double currentThrust_;

};

class DirectEngineModel: public EngineModel
{
public:

    DirectEngineModel(
            const boost::function< double( ) > specificImpulseFunction,
            const boost::function< double( ) > massFlowFunction ):
        specificImpulseFunction_( specificImpulseFunction ),
        massFlowFunction_( massFlowFunction )
    {

    }
\    void updateEngineModel( const double currentTime )
     {
         currentThrust_ = specificImpulseFunction_( ) * massFlowFunction_( );
     }

protected:

    boost::function< double( ) > specificImpulseFunction_;

    boost::function< double( ) > massFlowFunction_;
};

}

}

#endif // ENGINEMODEL_H
