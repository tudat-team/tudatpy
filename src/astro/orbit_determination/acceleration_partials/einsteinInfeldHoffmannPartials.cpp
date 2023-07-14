/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/einsteinInfeldHoffmannPartials.h"

namespace tudat
{

namespace acceleration_partials
{


    void EihEquationsPartials::getAccelerationWrtGravitationalParameter( Eigen::MatrixXd& muPartial, const int bodyIndex, const int muIndex )
    {
        muPartial.setZero( 3, 1 );
        double currentMu = eihEquations_->getGravitationalParameter( muIndex     );
        if( currentMu == 0.0 )
        {
            throw std::runtime_error( "Error when computing EIH partial w.r.t. mu, value of mu is 0" );
        }

        if( bodyIndex != muIndex )
        {
            muPartial += eihEquations_->getSinglePointMassAccelerations( bodyIndex, muIndex )  *
                ( 1.0 + physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getTotalScalarTermCorrection( bodyIndex, muIndex ) );
            muPartial += eihEquations_->getSingleSourceLocalPotential( bodyIndex, muIndex ) *
                physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * eihEquations_->getTotalVectorTermCorrection( bodyIndex, muIndex );
        }

        double scalarMultiplier0 = eihEquations_->getScalarTermMultiplier( 0 );
        double scalarMultiplier1 = eihEquations_->getScalarTermMultiplier( 1 );
        double scalarMultiplier6 = eihEquations_->getScalarTermMultiplier( 6 );
        double vectorMultiplier2 = eihEquations_->getVectorTermMultiplier( 2 );

        double scalarTerm0Partial = 0.0;
        if( bodyIndex != muIndex )
        {
            scalarTerm0Partial = eihEquations_->getSingleSourceLocalPotential( bodyIndex, muIndex );
        }

        for( int j = 0; j < numberOfExertingBodies_; j++ )
        {
            if ( j != bodyIndex && j != muIndex )
            {
                muPartial += (
                    eihEquations_->getSinglePointMassAccelerations( bodyIndex, j ) * (
                    scalarMultiplier0 * scalarTerm0Partial +
                    scalarMultiplier1 * eihEquations_->getSingleSourceLocalPotential( j, muIndex ) +
                    scalarMultiplier6 * eihEquations_->getRelativePositions( bodyIndex, j ).dot( eihEquations_->getSinglePointMassAccelerations( j, muIndex ) )
                    ) +
                    vectorMultiplier2 * eihEquations_->getSingleSourceLocalPotential( bodyIndex, j ) * eihEquations_->getSinglePointMassAccelerations( j, muIndex )
                                                                                               ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
            }
        }
        muPartial /= currentMu;
    }



} // namespace acceleration_partials

} // namespace tudat

