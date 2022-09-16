
#include <cmath>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/propulsion/costateBasedThrustGuidance.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"

namespace tudat
{

namespace propulsion
{

////! Constructor
//MeeCostateBasedThrustGuidance::MeeCostateBasedThrustGuidance(
//        const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
//        const std::function< Eigen::Vector6d( ) > centralBodyStateFunction,
//        const std::function< double( ) > centralBodyGravitationalParameterFunction,
//        std::function< Eigen::VectorXd( const double ) > costateFunction,
//        const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection )
//    : BodyFixedForceDirectionGuidance( bodyFixedForceDirection ),
//      thrustingBodyStateFunction_( thrustingBodyStateFunction ),
//      centralBodyStateFunction_( centralBodyStateFunction ),
//      centralBodyGravitationalParameterFunction_( centralBodyGravitationalParameterFunction ),
//      costateFunction_( costateFunction ){ }

////! Function to update the force direction to the current time.
//void MeeCostateBasedThrustGuidance::updateForceDirection( const double time )
//{
//    if( !( time == currentTime_ ) )
//    {
//        Eigen::VectorXd costates_ = costateFunction_( time );

//        // Get the current state in cartesian coordinates and keplerian elements, and some convenient parameters
//        Eigen::Vector6d currentState = thrustingBodyStateFunction_( ) - centralBodyStateFunction_( );
//        double centralBodyGravitationalParameter = centralBodyGravitationalParameterFunction_( );

//        // Obtain ModifiedEquinoctial elements, flag of 0 indicates that singularity occurs at 180 deg inclination.
//        Eigen::Vector6d modifiedEquinoctialElements =
//                orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
//                    currentState, centralBodyGravitationalParameter, 0 );

//        // Optimal control laws local variables declared for clarity
//        double auxiliaryParameterW = ( 1.0 + modifiedEquinoctialElements( 1 ) * cos( modifiedEquinoctialElements( 5 ) )
//                     + modifiedEquinoctialElements( 2 ) * sin( modifiedEquinoctialElements( 5 ) ) );
//        double auxiliaryParameterSSquared = 1.0 + modifiedEquinoctialElements( 3 ) * modifiedEquinoctialElements( 3 )
//                + modifiedEquinoctialElements( 4 ) * modifiedEquinoctialElements( 4 );

//        // Local variables for al constant terms for the calculation of pitch angle
//        double Lap = costates_( 0 ) * 2.0 * modifiedEquinoctialElements( 0 ) / auxiliaryParameterW;
//        double Laf1 = costates_( 1 )  * sin( modifiedEquinoctialElements( 5 ) ) ;
//        double Laf2 = costates_( 1 ) / auxiliaryParameterW *
//                ( ( auxiliaryParameterW + 1.0 ) * cos( modifiedEquinoctialElements( 5 ) )
//                                             + modifiedEquinoctialElements( 1 ) ) ;
//        double Lag1 = costates_( 2 ) * cos( modifiedEquinoctialElements( 5 ) );
//        double Lag2 = costates_( 2 ) / auxiliaryParameterW *
//                ( ( auxiliaryParameterW + 1.0 ) * sin( modifiedEquinoctialElements( 5 ) )
//                                             + modifiedEquinoctialElements( 2 ) );

//        // Calculate pitch angle, NOTE: denomitator ommitted since it is not relevant for the atan2 function,
//        // since both denominators are the same.
//        double thrustAngleAlpha = std::atan2( -Laf1+Lag1, -Lap-Laf2-Lag2);

//        // Local variables for al constant terms for the calculation of yaw angle
//        double Lbp = costates_( 0 ) * 2.0 * modifiedEquinoctialElements( 0 ) * cos( thrustAngleAlpha) / auxiliaryParameterW;
//        double Lbf1 = costates_( 1 )  * sin( modifiedEquinoctialElements( 5 ) ) * sin( thrustAngleAlpha );
//        double Lbf2 = costates_( 1 ) / auxiliaryParameterW *
//                ( ( auxiliaryParameterW + 1.0 ) * cos( modifiedEquinoctialElements( 5 ) )
//                                             + modifiedEquinoctialElements( 1 ) ) * cos( thrustAngleAlpha );
//        double Lbf3 = costates_( 1 ) / auxiliaryParameterW * ( modifiedEquinoctialElements( 2 ) *(
//                                                 modifiedEquinoctialElements( 3 ) * sin( modifiedEquinoctialElements( 5 ) )
//                                                 - modifiedEquinoctialElements( 4 ) * cos( modifiedEquinoctialElements( 5 ) ) ) );

//        double Lbg1 = costates_( 2 ) * cos( modifiedEquinoctialElements( 5 ) ) * sin( thrustAngleAlpha);
//        double Lbg2 = costates_( 2 ) / auxiliaryParameterW *
//                ( ( auxiliaryParameterW + 1.0 ) * sin( modifiedEquinoctialElements( 5 ) )
//                                             + modifiedEquinoctialElements( 2 ) ) * cos( thrustAngleAlpha);
//        double Lbg3 = costates_( 2 ) / auxiliaryParameterW * ( modifiedEquinoctialElements( 1 ) *(
//                                                 modifiedEquinoctialElements( 3 ) * sin( modifiedEquinoctialElements( 5 ) )
//                                                 - modifiedEquinoctialElements( 4 ) * cos( modifiedEquinoctialElements( 5 ) ) ) );
//        double Lbh = costates_( 3 ) * auxiliaryParameterSSquared *
//                cos( modifiedEquinoctialElements( 5 ) ) / ( 2.0 * auxiliaryParameterW );
//        double Lbk = costates_( 4 ) * auxiliaryParameterSSquared *
//                sin( modifiedEquinoctialElements( 5 ) ) / ( 2.0 * auxiliaryParameterW );

//        // Calculate yaw angle, NOTE: denomitator ommitted since it is not relevant for the atan2 function,
//        // since both denominators are the same.
//        double thrustAngleBeta = std::atan2( Lbf3 - Lbg3 - Lbh - Lbk, - Lbp - Lbf1 - Lbf2 + Lbg1 - Lbg2 );

//        // Calculate thrust direction
//        currentForceDirection_ = reference_frames::getTnwToInertialRotation(
//                    currentState, Eigen::Vector6d::Zero( ), false ) *
//                ( ( Eigen::Vector3d( ) <<
//                    cos( thrustAngleAlpha ) * cos( thrustAngleBeta ), sin( thrustAngleAlpha ) * cos( thrustAngleBeta ) ,
//                    sin( thrustAngleBeta )  ).finished( ).normalized( ) );
//        currentTime_ = time;


////        // Switching function for the thrust magnitude.
////        double thrustMagnitudeSwitchingCondition = /*( 1.0 / thrustingBodyMassFunction_( ) ) **/
////                ( Lbp * cos( thrustAngleBeta ) + Lbh * sin( thrustAngleBeta ) + Lbk * sin( thrustAngleBeta )
////                + Lbf1 * cos( thrustAngleBeta ) + Lbf2 * cos( thrustAngleBeta ) - Lbf3 * sin( thrustAngleBeta )
////                - Lbg1 * cos( thrustAngleBeta ) + Lbg2 * cos( thrustAngleBeta ) + Lbg3 * sin( thrustAngleBeta ) );
////        if ( thrustMagnitudeSwitchingCondition <= 0.0 )
////        {
////            std::cout << "INSIDE THRUST DIRECTION FUNCTION, THRUST ON. " << "\n\n";
////        }
////        else
////        {
////            std::cout << "INSIDE THRUST DIRECTION FUNCTION, THRUST OFF. " << "\n\n";
////        }
//    }

//}

} // namespace propulsion

} // namespace tudat
