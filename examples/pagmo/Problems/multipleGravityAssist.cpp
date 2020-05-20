/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "multipleGravityAssist.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::transfer_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
using namespace tudat;
using namespace pagmo;

MultipleGravityAssist::MultipleGravityAssist(std::vector< std::vector< double > > &bounds,
                                             std::vector< int > flybySequence,
                                             const bool useTripTime ):
    problemBounds_( bounds ), useTripTime_( useTripTime )
{

    // Specify required parameters
    // Specify the number of legs and type of legs.
    numberOfLegs_ = flybySequence.size( );
    legTypeVector_.resize( numberOfLegs_ );
    legTypeVector_[ 0 ] = mga_Departure;
    legTypeVector_[ numberOfLegs_ - 1 ] = capture;

    for(int i = 1; i < numberOfLegs_ - 1; i++){
        legTypeVector_[ i ] = mga_Swingby;
    }

    // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
    ephemerisVector_.resize( numberOfLegs_ );
    gravitationalParameterVector_.resize( numberOfLegs_ );
    minimumPericenterRadii_.resize( numberOfLegs_ );
    for(int i = 0; i < numberOfLegs_; i++)
    {
        switch(flybySequence[ i ])
        {
        case( 1 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );
            gravitationalParameterVector_[ i ] = 2.2032E13;
            minimumPericenterRadii_[ i ] = 2639.7E3;
            break;
        case( 2 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::venus );
            gravitationalParameterVector_[ i ] = 3.24859E14;
            minimumPericenterRadii_[ i ] = 6251.8E3;
            break;
        case( 3 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );
            gravitationalParameterVector_[ i ] = 3.986004418E14;
            minimumPericenterRadii_[ i ] = 6578.1E3;
            break;
        case( 4 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );
            gravitationalParameterVector_[ i ] = 4.282837E13;
            minimumPericenterRadii_[ i ] = 3596.2E3;
            break;
        case( 5 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::jupiter );
            gravitationalParameterVector_[ i ] = 1.26686534E17;
            minimumPericenterRadii_[ i ] = 72000.0E3;
            break;
        case( 6 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::saturn );
            gravitationalParameterVector_[ i ] = 3.7931187E16;
            minimumPericenterRadii_[ i ] = 61000.0E3;
            break;
        case( 7 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::uranus );
            gravitationalParameterVector_[ i ] = 5.793939E15;
            minimumPericenterRadii_[ i ] = 26000.0E3;
            break;
        case( 8 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::neptune );
            gravitationalParameterVector_[ i ] = 6.836529E15;
            minimumPericenterRadii_[ i ] = 25000.0E3;
            break;
        case( 9 ):
            ephemerisVector_[ i ] = std::make_shared< ephemerides::ApproximatePlanetPositions >
                    ( ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::pluto );
            gravitationalParameterVector_[ i ] = 8.71E11;
            minimumPericenterRadii_[ i ] = 1395.0E3;
            break;
        default:
            std::cerr<<"Planet in flyby sequence is not defined.";
        }
    }

    // Create departure and capture variables.
    semiMajorAxes_.resize( 2 );
    eccentricities_.resize( 2 );
    semiMajorAxes_ << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities_ << 0., 0.98;

}
//! Descriptive name of the problem
std::string MultipleGravityAssist::get_name() const {
    return "MGA transfer trajectory";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > MultipleGravityAssist::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const{
    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs_ + 1 );

    double TOF = 0;
    for(int i = 0; i < numberOfLegs_ ; i++){
        variableVector[ i ] = xv[ i ];
        if( i > 0 ){
            TOF += xv[i];
        }
    }
    variableVector[ numberOfLegs_ ] = 1;//dummy
    variableVector *= physical_constants::JULIAN_DAY;

    // Create the trajectory problem.
    Trajectory mgaTraj( numberOfLegs_, legTypeVector_, ephemerisVector_,
                          gravitationalParameterVector_, variableVector, sunGravitationalParameter,
                          minimumPericenterRadii_, semiMajorAxes_, eccentricities_ );

    // Start the deltaV vector.
    double resultingDeltaV;
    mgaTraj.calculateTrajectory( resultingDeltaV );

    if (std::isnan(resultingDeltaV))
    {
        resultingDeltaV = 1.0E10;
    }

    if ( useTripTime_ ){
        return { resultingDeltaV, TOF };
    }
    else {
        return { resultingDeltaV };
    }

}


