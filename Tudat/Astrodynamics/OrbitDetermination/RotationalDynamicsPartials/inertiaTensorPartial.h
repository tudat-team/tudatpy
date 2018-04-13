/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
#define TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H

#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/torquePartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class InertiaTensorPartial
{
public:

    InertiaTensorPartial(
            boost::function< double( ) > inertiaTensorNormalizationFunction ):
        inertiaTensorNormalizationFunction_( inertiaTensorNormalizationFunction )
    {
        unscaledPartialWrtC20_.setZero( );
        unscaledPartialWrtC20_( 0, 0 ) = 1.0 / 3.0;
        unscaledPartialWrtC20_( 1, 1 ) = 1.0 / 3.0;
        unscaledPartialWrtC20_( 2, 2 ) = 2.0 / 3.0;

        unscaledPartialWrtC21_.setZero( );
        unscaledPartialWrtC21_( 2, 0 ) = -1.0;
        unscaledPartialWrtC21_( 0, 2 ) = -1.0;

        unscaledPartialWrtC22_.setZero( );
        unscaledPartialWrtC20_( 0, 0 ) = -2.0;
        unscaledPartialWrtC20_( 1, 1 ) = -2.0;

        unscaledPartialWrtS21_.setZero( );
        unscaledPartialWrtS21_( 2, 1 ) = -1.0;
        unscaledPartialWrtS21_( 1, 2 ) = -1.0;

        unscaledPartialWrtS22_.setZero( );
        unscaledPartialWrtS22_( 0, 1 ) = -2.0;
        unscaledPartialWrtS22_( 1, 0 ) = -2.0;
        unscaledPartialWrtMeanInertia_ = Eigen::Matrix3d::Identity( );
    }

    ~InertiaTensorPartial( ){ }


    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        boost::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    std::pair< boost::function< void( std::vector< Eigen::MatrixXd >& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< void( std::vector< Eigen::MatrixXd >& ) > partialFunction;
        return std::macke_pair( partialFunction, 0 );
    }

    void update( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentNormalizationFactor_ = inertiaTensorNormalizationFunction_( );
            currentTime_ = currentTime;
        }
    }

protected:

    double currentTime_;

    double currentNormalizationFactor_;

    Eigen::Matrix3d unscaledPartialWrtC20_;

    Eigen::Matrix3d unscaledPartialWrtC21_;

    Eigen::Matrix3d unscaledPartialWrtC22_;

    Eigen::Matrix3d unscaledPartialWrtS21_;

    Eigen::Matrix3d unscaledPartialWrtS22_;

    Eigen::Matrix3d unscaledPartialWrtMeanInertia_;


    boost::function< double( ) > inertiaTensorNormalizationFunction_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUEPARTIALS_H
