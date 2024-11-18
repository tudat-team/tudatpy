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

#include "tudat/math/integrators/bulirschStoerVariableStepsizeIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
std::vector< unsigned int > getBulirschStoerStepSequence(
        const ExtrapolationMethodStepSequences& extrapolationMethodStepSequenceType,
        const unsigned int lengthOfSequence )
{
    std::vector< unsigned int > stepSequence;
    switch ( extrapolationMethodStepSequenceType )
    {
    case bulirsch_stoer_sequence:
        using namespace boost::assign;
        stepSequence += 2, 4, 6;
        for ( unsigned int i = 3; i < lengthOfSequence; i++ )
        {
            stepSequence.push_back(
                        2 * stepSequence.at( i - 2 ) );
        }
        break;
    case deufelhard_sequence:
        for ( unsigned int i = 0; i < lengthOfSequence; i++ )
        {
            stepSequence.push_back( 2 * ( i + 1 ) );
        }
        break;
    default: // The default case will never occur because sequence is an enum
        throw std::runtime_error( "Error, did not recognize step sequence" );
    }

    return stepSequence;
}

//template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
//template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
//template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


} // namespace numerical_integrators

} // namespace tudat
