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

#ifndef COMPAREESTIMATIONANDCOVARIANCERESULTSTESTCASE_H
#define COMPAREESTIMATIONANDCOVARIANCERESULTSTESTCASE_H

#include <boost/test/test_tools.hpp>

#include "tudat/simulation/estimation_setup/orbitDeterminationManagerForwardDeclarations.h"

namespace tudat
{
namespace unit_tests
{

template< typename TimeType = double, typename StateScalarType = double >
void compareEstimationAndCovarianceResults(
        const std::shared_ptr< simulation_setup::EstimationOutput< StateScalarType, TimeType > > estimationOutput,
        const std::shared_ptr< simulation_setup::CovarianceAnalysisOutput< StateScalarType, TimeType > > covarianceOutput )
{
    for( int i = 0; i < estimationOutput->getCorrelationMatrix( ).rows( ); i++ )
    {
        BOOST_CHECK_EQUAL( estimationOutput->getFormalErrorVector( )( i ), covarianceOutput->getFormalErrorVector( )( i ) );
        BOOST_CHECK_EQUAL( estimationOutput->getNormalizationTerms( )( i ), covarianceOutput->getNormalizationTerms( )( i ) );
        for( int j = 0; j < estimationOutput->getCorrelationMatrix( ).cols( ); j++ )
        {
            BOOST_CHECK_EQUAL( estimationOutput->getCorrelationMatrix( )( i, j ), covarianceOutput->getCorrelationMatrix( )( i, j ) );
        }
    }
}

}  // namespace unit_tests

}  // namespace tudat

#endif  // COMPAREESTIMATIONANDCOVARIANCERESULTSTESTCASE_H
