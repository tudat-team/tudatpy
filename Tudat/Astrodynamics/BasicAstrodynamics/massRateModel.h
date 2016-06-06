#ifndef MASSRATEMODEL_H
#define MASSRATEMODEL_H

#include <map>
#include <vector>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace basic_astrodynamics
{

class MassRateModel
{
public:

    virtual ~MassRateModel( ) { }

    virtual double getMassRate( )
    {
        return currentMassRate_;
    }

    virtual void updateMembers( const double currentTime = TUDAT_NAN ) = 0;

    virtual void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
    }

protected:

    double currentTime_;

    double currentMassRate_;

protected:

private:
};

class CustomMassRateModel: public MassRateModel
{
public:

    CustomMassRateModel(
            const boost::function< double( const double ) > massRateFunction ){ }

    ~CustomMassRateModel( ){ }

    virtual void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime == currentTime ) )
        {
            currentMassRate_ = massRateFunction( currentTime );
        }
    }

protected:

    boost::function< double( const double ) > massRateFunction ;

protected:

private:
};


}

}
#endif // MASSRATEMODEL_H
