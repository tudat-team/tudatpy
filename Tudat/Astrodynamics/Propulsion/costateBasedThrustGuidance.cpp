
#include <cmath>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/Propulsion/costateBasedThrustGuidance.h"

namespace tudat
{

namespace propulsion
{

MeeCostateBasedThrustGuidance::MeeCostateBasedThrustGuidance(
        const boost::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
        const boost::function< Eigen::Vector6d( ) > centralBodyStateFunction,
        const boost::function< double( ) > centralBodyGravitationalParameterFunction,
        boost::function< Eigen::VectorXd( const double ) > costateFunction,
        const boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection )
    : BodyFixedForceDirectionGuidance( bodyFixedForceDirection ),
      thrustingBodyStateFunction_( thrustingBodyStateFunction ),
      centralBodyStateFunction_( centralBodyStateFunction ),
      centralBodyGravitationalParameterFunction_( centralBodyGravitationalParameterFunction ),
      costateFunction_( costateFunction ){ }

void MeeCostateBasedThrustGuidance::updateForceDirection( const double time )
{
    Eigen::VectorXd costates_ = costateFunction_( time );

    // Get the current state in cartesian coordinates and keplerian elements, and some convenient parameters
    Eigen::Vector6d currentState = thrustingBodyStateFunction_( ) - centralBodyStateFunction_( );
    double centralBodyGravitationalParameter = centralBodyGravitationalParameterFunction_( );

    // Obtain ModifiedEquinoctial elements, flag of 0 indicates that singularity occurs at 180 deg inclination.
    Eigen::Vector6d modifiedEquinoctialElements =
            orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                currentState, centralBodyGravitationalParameter, 0 );

    // Optimal control laws local variables declared for clarity
    double w = ( 1.0 + modifiedEquinoctialElements( 1 ) * cos( modifiedEquinoctialElements( 5 ) )
                 + modifiedEquinoctialElements( 2 ) * sin( modifiedEquinoctialElements( 5 ) ) );
    double sSquared = 1.0 + modifiedEquinoctialElements( 3 ) * modifiedEquinoctialElements( 3 )
            + modifiedEquinoctialElements( 4 ) * modifiedEquinoctialElements( 4 );

    // Local variables for al constant terms for the calculation of pitch angle
    double Lap = costates_( 0 ) * 2.0 * modifiedEquinoctialElements( 0 ) / w;
    double Laf1 = costates_( 1 )  * sin( modifiedEquinoctialElements( 5 ) ) ;
    double Laf2 = costates_( 1 ) / w * ( ( w + 1.0 ) * cos( modifiedEquinoctialElements( 5 ) )
                                         + modifiedEquinoctialElements( 1 ) ) ;
    double Lag1 = costates_( 2 ) * cos( modifiedEquinoctialElements( 5 ) );
    double Lag2 = costates_( 2 ) / w * ( ( w + 1.0 ) * sin( modifiedEquinoctialElements( 5 ) )
                                         + modifiedEquinoctialElements( 2 ) );

    // Calculate pitch angle, NOTE: denomitator ommitted since it is not relevant for the atan2 function,
    // since both denominators are the same.
    double alpha = std::atan2( -Laf1+Lag1, -Lap-Laf2-Lag2);

    // Local variables for al constant terms for the calculation of yaw angle
    double Lbp = costates_( 0 ) * 2.0 * modifiedEquinoctialElements( 0 ) * cos( alpha) / w;
    double Lbf1 = costates_( 1 )  * sin( modifiedEquinoctialElements( 5 ) ) * sin( alpha );
    double Lbf2 = costates_( 1 ) / w * ( ( w + 1.0 ) * cos( modifiedEquinoctialElements( 5 ) )
                                         + modifiedEquinoctialElements( 1 ) ) * cos( alpha );
    double Lbf3 = costates_( 1 ) / w * ( modifiedEquinoctialElements( 2 ) *(
                                             modifiedEquinoctialElements( 3 ) * sin( modifiedEquinoctialElements( 5 ) )
                                             - modifiedEquinoctialElements( 4 ) * cos( modifiedEquinoctialElements( 5 ) ) ) );

    double Lbg1 = costates_( 2 ) * cos( modifiedEquinoctialElements( 5 ) ) * sin( alpha);
    double Lbg2 = costates_( 2 ) / w * ( ( w + 1.0 ) * sin( modifiedEquinoctialElements( 5 ) )
                                         + modifiedEquinoctialElements( 2 ) ) * cos( alpha);
    double Lbg3 = costates_( 2 ) / w * ( modifiedEquinoctialElements( 1 ) *(
                                             modifiedEquinoctialElements( 3 ) * sin( modifiedEquinoctialElements( 5 ) )
                                             - modifiedEquinoctialElements( 4 ) * cos( modifiedEquinoctialElements( 5 ) ) ) );
    double Lbh = costates_( 3 ) * sSquared  * cos( modifiedEquinoctialElements( 5 ) ) / ( 2.0 * w );;
    double Lbk = costates_( 4 ) * sSquared  * sin( modifiedEquinoctialElements( 5 ) ) / ( 2.0 * w );;

    // Calculate yaw angle, NOTE: denomitator ommitted since it is not relevant for the atan2 function,
    // since both denominators are the same.
    double beta = std::atan2( Lbf3 - Lbg3 - Lbh - Lbk, - Lbp - Lbf1 - Lbf2 + Lbg1 - Lbg2 );

    // Calculate thrust direction
    currentForceDirection_ = reference_frames::getVelocityBasedLvlhToInertialRotation(
                currentState, Eigen::Vector6d::Zero( ), false, true ) *
            ( ( Eigen::Vector3d( ) << cos( alpha ) * cos( beta ), -sin( beta ) ,
                -sin( alpha ) * cos( beta ) ).finished( ).normalized( ) );

}

} // namespace propulsion

} // namespace tudat
