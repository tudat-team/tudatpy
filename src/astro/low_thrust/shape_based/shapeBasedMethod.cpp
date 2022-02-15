/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/low_thrust/shape_based/shapeBasedMethod.h"

namespace tudat
{

namespace shape_based_methods
{




//! Returns state history.
void ShapeBasedMethod::getTrajectory(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    propagatedTrajectory.clear( );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the mass profile of a shape-based trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        double currentIndependentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );
        propagatedTrajectory[ epochsVector[ i ] ] = computeCurrentStateVector( currentIndependentVariable );
    }
}


////! Return thrust profile.
//void ShapeBasedMethod::getThrustForceProfile(
//        std::vector< double >& epochsVector,
//        std::map< double, Eigen::VectorXd >& thrustProfile,
//        std::function< double ( const double ) > specificImpulseFunction,
//        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )

//{
//    thrustProfile.clear( );
//    std::map< double, Eigen::VectorXd > massProfile;

//    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

//    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
//    {
//        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
//        {
//            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectory, "
//                                      "epochs are not provided in increasing order." );
//        }

//        double independentVariable = convertTimeToAzimuth( epochsVector[ i ] );

//        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];
//        thrustProfile[ epochsVector[ i ] ] = currentMass * computeCurrentThrustAccelerationMagnitudeFromAzimuth(
//                    independentVariable )
//                * computeCurrentThrustAccelerationDirectionFromAzimuth( independentVariable );
//    }
//}




} // namespace transfer_trajectories

} // namespace tudat
