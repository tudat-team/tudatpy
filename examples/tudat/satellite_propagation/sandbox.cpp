/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <functional>

#include <tudat/simulation/simulation.h>

#include "tudat/io/applicationOutput.h"

using namespace tudat;

double findLocalMinimumOfTargetDistance(
        const double lowerBound,
        const double upperBound,
        const double tolerance,
        const std::function< Eigen::Vector3d( const double )> positionFunction )
{
    double currentUpperBound = upperBound;
    double currentLowerBound = lowerBound;
    double newTestValue = ( currentUpperBound + currentLowerBound ) / 2.0;
    double currentTestValue;

    double currentUpperDistance, currentLowerDistance, currentTestDistance ;
    int currentUpperDerivative, currentLowerDerivative, currentTestDerivative;

    int counter = 0;
    do
    {

        currentTestValue = newTestValue;

        currentUpperDistance = positionFunction( currentUpperBound ).norm( );
        currentLowerDistance  = positionFunction( currentLowerBound ).norm( );
        currentTestDistance = positionFunction( currentTestValue ).norm( );

        currentUpperDerivative = utilities::sgn< double >( ( currentUpperDistance - positionFunction( currentUpperBound - tolerance ).norm( ) ) / tolerance );
        currentLowerDerivative  = utilities::sgn< double >( ( currentLowerDistance - positionFunction( currentLowerBound - tolerance ).norm( ) ) / tolerance );
        currentTestDerivative = utilities::sgn< double >( ( currentTestDistance - positionFunction( currentTestValue - tolerance ).norm( ) ) / tolerance );

        //        if( currentUpperDerivative == currentLowerDerivative )
        //        {
        //            std::cerr<<"Error 1 when finding minium target distance"<<std::endl;
        //            std::cout<<std::setprecision( 12 )<<currentLowerBound<<" "<<currentTestValue<<" "<<currentUpperBound<<" "<<
        //                       currentLowerDistance<<" "<<currentTestDistance<<" "<<currentUpperDistance<<" "<<
        //                       ( ( currentLowerDistance - positionFunction( currentLowerBound - tolerance ).norm( ) ) / tolerance )<<" "<<
        //                       ( ( currentTestDistance - positionFunction( currentTestValue - tolerance ).norm( ) ) / tolerance )<<" "<<
        //                       ( ( currentUpperDistance - positionFunction( currentUpperBound - tolerance ).norm( ) ) / tolerance )<<" "<<counter<<std::endl<<std::endl;
        //        }

        if(  currentUpperDerivative > 0 && currentTestDerivative < 0 )
        {

            newTestValue = ( currentUpperBound + currentTestValue ) / 2.0;
            currentLowerBound = currentTestValue;
        }
        else if( currentLowerDerivative < 0 && currentTestDerivative > 0 )
        {
            newTestValue = ( currentTestValue + currentLowerBound ) / 2.0;
            currentUpperBound = currentTestValue;

        }
        else
        {

//            std::cout<<std::setprecision( 12 )<<currentLowerBound<<" "<<currentTestValue<<" "<<currentUpperBound<<" "<<
//                       currentLowerDistance<<" "<<currentTestDistance<<" "<<currentUpperDistance<<" "<<
//                       ( ( currentLowerDistance - positionFunction( currentLowerBound - tolerance ).norm( ) ) / tolerance )<<" "<<
//                       ( ( currentTestDistance - positionFunction( currentTestValue - tolerance ).norm( ) ) / tolerance )<<" "<<
//                       ( ( currentUpperDistance - positionFunction( currentUpperBound - tolerance ).norm( ) ) / tolerance )<<" "<<counter<<std::endl<<std::endl;

//            std::cerr<<"Error 1 when finding minium target distance"<<std::endl;
        }

        counter ++;
        if( counter > 1E5 )
        {
            std::cerr<<"Error 2 when finding minium target distance"<<std::endl;
        }

    } while( std::fabs( newTestValue - currentTestValue ) > tolerance );

    return newTestValue;
}

std::map< double, double > findLocalMinimaOfTargetDistance(
        const double lowerBound,
        const double upperBound,
        const double threshold,
        const double tolerance,
        const double initialSearchTimeStep,
        const std::function< Eigen::Vector3d( const double )> positionFunction )
{
    std::map< double, double > minima;

    double upperTime = lowerBound + 2.0 * initialSearchTimeStep,
            middleTime = lowerBound + initialSearchTimeStep,
            lowerTime = lowerBound;

    double upperValue, middleValue, lowerValue;

    double candidateTime;

    while( upperTime <= upperBound )
    {
        upperValue = ( positionFunction( upperTime ) ).norm( );
        middleValue = ( positionFunction( middleTime ) ).norm( );
        lowerValue = ( positionFunction( lowerTime ) ).norm( );

        if( utilities::sgn( upperValue - middleValue ) > 0 &&
                utilities::sgn( middleValue - lowerValue ) < 0 )
        {
            candidateTime = findLocalMinimumOfTargetDistance( lowerTime, upperTime, tolerance, positionFunction );

            if( positionFunction( candidateTime ).norm( ) <= threshold )
            {
                minima[ candidateTime ] = positionFunction( candidateTime ).norm( );
            }
        }

        upperTime += initialSearchTimeStep;
        middleTime += initialSearchTimeStep;
        lowerTime += initialSearchTimeStep;


    }
    return minima;
}



void getCloseApproachTimes(
        const double initialTime, const double finalTime, const double approachThreshold )
{

        std::map< double, double > closeApproachTimes = findLocalMinimaOfTargetDistance(
                    initialTime, finalTime, approachThreshold, 5.0, 1800.0,
                    std::bind( spice_interface::getBodyCartesianPositionAtEpoch, "-28", " Callisto",
                               "ECLIPJ2000", "None", std::placeholders::_1 ) );

        for( auto it : closeApproachTimes )
        {
            std::cout<<std::setprecision( 16 )<<it.first<<" "<<it.second<<" "<<
                       spice_interface::getBodyCartesianStateAtEpoch(
                           "-28", " Callisto", "ECLIPJ2000", "None", it.first ).segment( 3, 3 ).norm( )<<std::endl;
        }
        std::cout<<std::endl;
        for( auto it : closeApproachTimes )
        {
            std::cout<<std::setprecision( 16 )<<it.first<<" "<<spice_interface::getBodyCartesianStateAtEpoch(
                           "-28", " Callisto", "ECLIPJ2000", "None", it.first ).transpose( )<<std::endl;
        }

}
//! Execute propagation of orbit of Asterix around the Earth.
int main( )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    double initialTime = 30.0 * tudat::physical_constants::JULIAN_YEAR;
    double finalTime = 33.0 * tudat::physical_constants::JULIAN_YEAR;
    double approachThreshold = 2.0E7;
    getCloseApproachTimes( initialTime, finalTime, approachThreshold );



}

