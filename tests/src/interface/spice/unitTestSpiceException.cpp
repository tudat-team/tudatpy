/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceException.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_spice_exception )

BOOST_AUTO_TEST_CASE( testImplementedSpiceExceptions )
{
    // The SPICE exceptions have been implemented from a list of exceptions
    // that are thrown by the SPICE library. This test checks that the
    // exceptions are thrown correctly and no off-by one errors are present.

    // First exception in the list
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(ADDRESSOUTOFBOUNDS)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceADDRESSOUTOFBOUNDS );

    // Arbitrary exceptions in the list
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(FILENOTFOUND)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceFILENOTFOUND );
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(INSUFFICIENTDATA)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceINSUFFICIENTDATA );
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(IDCODENOTFOUND)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceIDCODENOTFOUND );
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(NOLEAPSECONDS)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceNOLEAPSECONDS );

    // Last exception in the list
    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "SPICE(ZZHOLDDGETFAILED)", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceZZHOLDDGETFAILED );
}

BOOST_AUTO_TEST_CASE( testNotImplementedSpiceExceptions )
{
    // The

    BOOST_CHECK_THROW( tudat::exceptions::throwSpiceException( "GIBBERISH", "Explanation", "Long Message", "Traceback" ),
                       tudat::exceptions::SpiceError );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
