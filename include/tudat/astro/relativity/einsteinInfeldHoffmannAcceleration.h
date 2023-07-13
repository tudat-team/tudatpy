#ifndef EINSTEININFELDHOFFMANNACCELERATION_H
#define EINSTEININFELDHOFFMANNACCELERATION_H

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"

namespace tudat
{

namespace relativity
{

class EinsteinInfeldHoffmannAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    EinsteinInfeldHoffmannAcceleration( const std::shared_ptr< EinsteinInfeldHoffmannEquations > eihEquations,
                                        const std::string bodyUndergoingAcceleration,
                                        const std::vector< std::string > bodiesExertingAcceleration ):
        eihEquations_( eihEquations ), bodyUndergoingAcceleration_( bodyUndergoingAcceleration ),
        bodiesExertingAcceleration_( bodiesExertingAcceleration ){ }


    virtual void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( ! ( currentTime == currentTime_ ) )
        {
            eihEquations_->update( currentTime );
            currentTime_ = currentTime;
            if( currentTime == currentTime_ )
            {
                currentAcceleration_ = eihEquations_->getAccelerationOfBody( bodyUndergoingAcceleration_ );
            }
        }
    }

    std::vector< std::string > getBodiesExertingAcceleration( )
    {
        return bodiesExertingAcceleration_;
    }

    std::shared_ptr< EinsteinInfeldHoffmannEquations > getEihEquations( )
    {
        return eihEquations_;
    }

    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
        eihEquations_->resetCurrentTime( );
    }

private:

    std::shared_ptr< EinsteinInfeldHoffmannEquations > eihEquations_;
    std::string bodyUndergoingAcceleration_;
    std::vector< std::string > bodiesExertingAcceleration_;
};

}

}
#endif // EINSTEININFELDHOFFMANNACCELERATION_H
