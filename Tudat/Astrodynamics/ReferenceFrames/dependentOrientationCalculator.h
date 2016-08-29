#ifndef TUDAT_DEPENDENTORIENTATIONCALCULATOR_H
#define TUDAT_DEPENDENTORIENTATIONCALCULATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{

namespace reference_frames
{


class DependentOrientationCalculator
{
public:

    DependentOrientationCalculator( ): currentTime_( TUDAT_NAN ){ }

    virtual ~DependentOrientationCalculator( ){ }

    virtual Eigen::Quaterniond getRotationToLocalFrame( ) = 0;

    virtual Eigen::Quaterniond getRotationToGlobalFrame( ) = 0;

    virtual void updateCalculator( const double currentTime ) = 0;

    Eigen::Quaterniond getRotationToLocalFrame( const double currentTime )
    {
        updateCalculator( currentTime );
        return getRotationToLocalFrame( );
    }

    Eigen::Quaterniond getRotationToGlobalFrame( const double currentTime )
    {
        updateCalculator( currentTime );
        return getRotationToGlobalFrame( );
    }

    virtual void resetDerivedClassTime( const double currentTime = TUDAT_NAN ){ }

    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        resetDerivedClassTime( currentTime );
    }

protected:

    double currentTime_;
};

} // namespace reference_frames

} // namespace tudat

#endif // TUDAT_DEPENDENTORIENTATIONCALCULATOR_H
