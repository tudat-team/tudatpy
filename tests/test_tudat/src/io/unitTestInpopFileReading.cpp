/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/io/readInpopEphemerisFile.h"

namespace tudat
{

using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_inpop_file_reader )

BOOST_AUTO_TEST_CASE( test_inpop_time_ephemeris )
{
    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string timeDifferenceFileName = spiceKernelsPath + "/inpop19a_TCB_m100_p100_asc/inpop19a_TCB_m100_p100_asc_pos_TCG.asc" ;
    std::map< double, std::vector< double > > timeEphemerisCoefficients = readInpopTimeEphemeris< double >( timeDifferenceFileName );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > timeEphemerisInterpolator =
            createInpopTimeEphemerisInterpolator( timeDifferenceFileName );

    BOOST_CHECK_SMALL( timeEphemerisInterpolator->interpolate( ( basic_astrodynamics::TAI_JULIAN_DAY_AT_TIME_SYNCHRONIZATION -
                                                         basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY ), 5.0E-13 );
}

BOOST_AUTO_TEST_CASE( test_inpop_time_ephemeris2 )
{
    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string timeDifferenceFileName = spiceKernelsPath + "/inpop19a_TDB_m100_p100_asc/inpop19a_TDB_m100_p100_asc_pos_TT.asc" ;
    std::map< double, std::vector< double > > timeEphemerisCoefficients = readInpopTimeEphemeris< double >( timeDifferenceFileName );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > timeEphemerisInterpolator =
            createInpopTimeEphemerisInterpolator( timeDifferenceFileName );
    
    BOOST_CHECK_SMALL( sofa_interface::getTDBminusTT( 1.0E7, 0.0, 0.0, 0.0, 0.0 )+timeEphemerisInterpolator->interpolate( 1.0E7 ), 5.0E-9 );
}

BOOST_AUTO_TEST_CASE( test_inpop_state_ephemeris )
{
    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string positionFileName = spiceKernelsPath + "/inpop19a_TCB_m100_p100_asc/inpop19a_TCB_m100_p100_asc_pos_Ear.asc" ;
    std::string velocityFileName = spiceKernelsPath + "/inpop19a_TCB_m100_p100_asc/inpop19a_TCB_m100_p100_asc_vel_Ear.asc" ;

    std::map< double, std::vector< Eigen::Vector6d > > stateCoefficients = readInpopStateEphemeris( positionFileName, velocityFileName );

    std::shared_ptr< ephemerides::Ephemeris > earthEphemeris = createInpopEphemerisFromFiles( positionFileName, velocityFileName );

    std::string kernelsPath = spiceKernelsPath;
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.bpc");

    Eigen::Vector6d spiceKernelState = spice_interface::getBodyCartesianStateAtEpoch( "Earth", "SSB", "J2000", "NONE",
                                                                                basic_astrodynamics::convertTcbToTdb( 0.0 ) );
    Eigen::Vector6d interpolatorState = earthEphemeris->getCartesianState( 0.0 );

    Eigen::Vector3d spicePosition = spiceKernelState.segment( 0, 3 );
    Eigen::Vector3d interpolatedPosition = ( 1.0 - physical_constants::LB_TIME_RATE_TERM ) * interpolatorState.segment( 0, 3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spicePosition, interpolatedPosition, 1.0E-14 );

    Eigen::Vector3d spiceVelocity = spiceKernelState.segment( 3, 3 );
    Eigen::Vector3d interpolatedVelocity = interpolatorState.segment( 3, 3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceVelocity, interpolatedVelocity, 1.0E-14 );

}

BOOST_AUTO_TEST_CASE( test_inpop_state_ephemeris2 )
{
    std::string spiceKernelsPath = paths::getSpiceKernelPath( );
    std::string positionFileName = spiceKernelsPath + "/inpop19a_TDB_m100_p100_asc/inpop19a_TDB_m100_p100_asc_pos_Ear.asc" ;
    std::string velocityFileName = spiceKernelsPath + "/inpop19a_TDB_m100_p100_asc/inpop19a_TDB_m100_p100_asc_vel_Ear.asc" ;

    std::map< double, std::vector< Eigen::Vector6d > > stateCoefficients = readInpopStateEphemeris( positionFileName, velocityFileName );

    std::shared_ptr< ephemerides::Ephemeris > earthEphemeris = createInpopEphemerisFromFiles( positionFileName, velocityFileName );

    std::string kernelsPath = spiceKernelsPath;
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "/inpop19a_TDB_m100_p100_spice.bpc");

    Eigen::Vector6d spiceKernelState = spice_interface::getBodyCartesianStateAtEpoch( "Earth", "SSB", "J2000", "NONE", 0.0 );
    Eigen::Vector6d interpolatorState = earthEphemeris->getCartesianState( 0.0 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceKernelState, interpolatorState, 1.0E-15 );
}
BOOST_AUTO_TEST_SUITE_END( )

}
