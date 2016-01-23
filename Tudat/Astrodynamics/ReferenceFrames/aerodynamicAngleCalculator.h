#ifndef TUDAT_AERODYNAMICANGLECALCULATOR_H
#define TUDAT_AERODYNAMICANGLECALCULATOR_H

#include <vector>
#include <map>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace reference_frames
{

enum AerodynamicsReferenceFrames
{
    inertial_frame = -1,
    corotating_frame = 0,
    vertical_frame = 1,
    trajectory_frame = 2,
    aerodynamic_frame = 3,
    body_frame = 4
};

enum AerodynamicsReferenceFrameAngles
{
    latitude_angle,
    longitude_angle,
    heading_angle,
    flight_path_angle,
    angle_of_attack,
    angle_of_sideslip,
    bank_angle
};

class AerodynamicAngleCalculator
{
public:
    AerodynamicAngleCalculator(
            const boost::function< basic_mathematics::Vector6d( ) > bodyFixedStateFunction,
            const boost::function< double( ) > angleOfAttackFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > angleOfSideslipFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > bankAngleFunction =
            boost::lambda::constant ( 0.0 ),
            const bool calculateVerticalToAerodynamicFrame = 0 ):
        bodyFixedStateFunction_( bodyFixedStateFunction ),
        angleOfAttackFunction_( angleOfAttackFunction ),
        angleOfSideslipFunction_( angleOfSideslipFunction ),
        bankAngleFunction_( bankAngleFunction ),
        calculateVerticalToAerodynamicFrame_( calculateVerticalToAerodynamicFrame ){ }

    void update( );

    Eigen::Quaterniond getRotationQuaternionBetweenFrames(
            const AerodynamicsReferenceFrames originalFrame,
            const AerodynamicsReferenceFrames targetFrame );

    double getAerodynamicAngle( const AerodynamicsReferenceFrameAngles angleId );

private:
    std::map< AerodynamicsReferenceFrameAngles, double > currentAerodynamicAngles_;

    std::map< std::pair< AerodynamicsReferenceFrames, AerodynamicsReferenceFrames >, Eigen::Quaterniond > currentRotationMatrices_;

    basic_mathematics::Vector6d currentBodyFixedState_;

    boost::function< basic_mathematics::Vector6d( ) > bodyFixedStateFunction_;

    boost::function< double( ) > angleOfAttackFunction_;

    boost::function< double( ) > angleOfSideslipFunction_;

    boost::function< double( ) > bankAngleFunction_;

    bool calculateVerticalToAerodynamicFrame_;
};

}

}
#endif // TUDAT_AERODYNAMICANGLECALCULATOR_H
