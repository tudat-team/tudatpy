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

#include <memory>
#include <sstream>
#include <stdexcept>

#include "tudat/simulation/estimation_setup/estimationInterfacesForwardDeclarations.h"

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
        if( estimationOutput->getFormalErrorVector( )( i ) != covarianceOutput->getFormalErrorVector( )( i ) )
        {
            std::ostringstream error;
            error << "Formal error mismatch at index " << i << ".";
            throw std::runtime_error( error.str( ) );
        }
        if( estimationOutput->getNormalizationTerms( )( i ) != covarianceOutput->getNormalizationTerms( )( i ) )
        {
            std::ostringstream error;
            error << "Normalization term mismatch at index " << i << ".";
            throw std::runtime_error( error.str( ) );
        }
        for( int j = 0; j < estimationOutput->getCorrelationMatrix( ).cols( ); j++ )
        {
            if( estimationOutput->getCorrelationMatrix( )( i, j ) != covarianceOutput->getCorrelationMatrix( )( i, j ) )
            {
                std::ostringstream error;
                error << "Correlation matrix mismatch at (" << i << ", " << j << ").";
                throw std::runtime_error( error.str( ) );
            }
        }
    }
}

}  // namespace unit_tests

}  // namespace tudat

#endif  // COMPAREESTIMATIONANDCOVARIANCERESULTSTESTCASE_H
