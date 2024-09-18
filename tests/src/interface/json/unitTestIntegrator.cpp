/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/interface/json/unitTestSupport.h"
#include "tudat/interface/json/math/integrator.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_integrator )

// Test 1: integrator types
BOOST_AUTO_TEST_CASE( test_json_integrator_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            numerical_integrators::integratorTypes,
                            numerical_integrators::unsupportedIntegratorTypes );
}

// Test 2: Runge-Kutta coefficient sets
BOOST_AUTO_TEST_CASE( test_json_integrator_rksets )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "rksets" ),
                            numerical_integrators::rungeKuttaCoefficientSets,
                            numerical_integrators::unsupportedRungeKuttaCoefficientSets );
}

// Test 3: euler
BOOST_AUTO_TEST_CASE( test_json_integrator_euler )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::json_interface;

    // Create IntegratorSettings from JSON file
    const std::shared_ptr< IntegratorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< IntegratorSettings< double > > >( INPUT( "euler" ) );

    // Create IntegratorSettings manually
    const AvailableIntegrators integratorType = euler;
    const double initialTime = 3.0;
    const double stepSize = 1.4;
    const std::shared_ptr< IntegratorSettings< double > > manualSettings =
            std::make_shared< IntegratorSettings< double > >( integratorType,
                                                                initialTime,
                                                                stepSize );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: rungeKutta4
BOOST_AUTO_TEST_CASE( test_json_integrator_rungeKutta4 )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::json_interface;

    // Create IntegratorSettings from JSON file
    const std::shared_ptr< IntegratorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< IntegratorSettings< double > > >( INPUT( "rungeKutta4" ) );

    // Create IntegratorSettings manually
    const AvailableIntegrators integratorType = rungeKutta4;
    const double initialTime = 3.0;
    const double stepSize = 1.4;
    const unsigned int saveFrequency = 2;
    const bool assessTerminationConditionDuringIntegrationSubsteps = true;
    const std::shared_ptr< IntegratorSettings< double > > manualSettings =
            std::make_shared< IntegratorSettings< double > >( integratorType,
                                                                initialTime,
                                                                stepSize,
                                                                saveFrequency,
                                                                assessTerminationConditionDuringIntegrationSubsteps );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}


// Test 5: rungeKuttaVariableStepSize
BOOST_AUTO_TEST_CASE( test_json_integrator_rungeKuttaVariableStepSize )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::json_interface;

    // Create IntegratorSettings from JSON file
    const std::shared_ptr< IntegratorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< IntegratorSettings< double > > >( INPUT( "rungeKuttaVariableStepSize" ) );

    // Create IntegratorSettings manually
    const AvailableIntegrators integratorType = rungeKuttaVariableStepSize;
    const double initialTime = -0.3;
    const double initialStepSize = 1.4;
    const CoefficientSets rungeKuttaCoefficientSet =
            rungeKuttaFehlberg78;
    const double minimumStepSize = 0.4;
    const double maximumStepSize = 2.4;
    const double relativeErrorTolerance = 1.0E-4;
    const double absoluteErrorTolerance = 1.0E-2;
    const double safetyFactorForNextStepSize = 2.0;
    const double maximumFactorIncreaseForNextStepSize = 10.0;
    const double minimumFactorDecreaseForNextStepSize = 0.1;

    const std::shared_ptr< IntegratorSettings< double > > manualSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                initialTime,
                initialStepSize,
                rungeKuttaCoefficientSet,
                minimumStepSize,
                maximumStepSize,
                relativeErrorTolerance,
                absoluteErrorTolerance,
                1,
                false,
                safetyFactorForNextStepSize,
                maximumFactorIncreaseForNextStepSize,
                minimumFactorDecreaseForNextStepSize );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: adamsBashforthMoulton
BOOST_AUTO_TEST_CASE( test_json_integrator_adamsBashforthMoulton )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::json_interface;

    // Create IntegratorSettings from JSON file
    const std::shared_ptr< IntegratorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< IntegratorSettings< double > > >( INPUT( "adamsBashforthMoulton" ) );

    // Create IntegratorSettings manually
    const AvailableIntegrators integratorType = adamsBashforthMoulton;
    const double initialTime = -0.3;
    const double initialStepSize = 1.4;
    const double minimumStepSize = 0.4;
    const double maximumStepSize = 2.4;
    const double relativeErrorTolerance = 1.0E-4;
    const double absoluteErrorTolerance = 1.0E-2;
    const double bandwidth = 200;
    const int minimumOrder = 6;
    const int maximumOrder = 11;
    const std::shared_ptr< IntegratorSettings< double > > manualSettings =
            std::make_shared< AdamsBashforthMoultonSettings< double > >(  initialTime,
                                                                initialStepSize,
                                                                minimumStepSize,
                                                                maximumStepSize,
                                                                relativeErrorTolerance,
                                                                absoluteErrorTolerance,
                                                                minimumOrder,
                                                                maximumOrder,
                                                                1,
                                                                false,
                                                                bandwidth );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: bulirschStoer
BOOST_AUTO_TEST_CASE( test_json_integrator_bulirschStoer )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::json_interface;

    // Create IntegratorSettings from JSON file
    const std::shared_ptr< IntegratorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< IntegratorSettings< double > > >( INPUT( "bulirschStoer" ) );

    // Create IntegratorSettings manually
    const double initialTime = -0.3;
    const double initialStepSize = 1.4;
    const double minimumStepSize = 0.4;
    const double maximumStepSize = 2.4;
    const double relativeErrorTolerance = 1.0E-4;
    const double absoluteErrorTolerance = 1.0E-2;
    const int maximumNumberOfSteps = 8;

    const std::shared_ptr< IntegratorSettings< double > > manualSettings =
            std::make_shared< BulirschStoerIntegratorSettings< double > >(
                initialTime, initialStepSize, bulirsch_stoer_sequence, maximumNumberOfSteps,
                minimumStepSize, maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
