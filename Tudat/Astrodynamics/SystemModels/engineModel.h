#ifndef TUDAT_ENGINEMODEL_H
#define TUDAT_ENGINEMODEL_H

#include <Eigen/Core>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Propulsion/thrustFunctions.h"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace system_models
{

class EngineModel
{
public:

    EngineModel(
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
       currentThrust_( TUDAT_NAN ),
       bodyFixedThrustDirection_( bodyFixedThrustDirection ){ }

    double getCurrentThrust( )
    {
        return currentThrust_;
    }

    virtual ~EngineModel( ){ }

    virtual void updateEngineModel( const double currentTime ) = 0;

    virtual double getCurrentMassRate( ) = 0;

    Eigen::Vector3d getBodyFixedThrustDirection( )
    {
        return bodyFixedThrustDirection_;
    }

protected:

    double currentThrust_;

    Eigen::Vector3d bodyFixedThrustDirection_;

};

class DirectEngineModel: public EngineModel
{
public:

    DirectEngineModel(
            const double specificImpulse,
            const boost::function< double( ) > massFlowFunction,
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        EngineModel( bodyFixedThrustDirection ),
        specificImpulse_( specificImpulse ),
        massFlowFunction_( massFlowFunction )
    { }

    void updateEngineModel( const double currentTime )
    {
        currentThrust_ = propulsion::computeThrustFromSpecificImpulse( massFlowFunction_( ), specificImpulse_ );
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

#endif // TUDAT_ENGINEMODEL_H
