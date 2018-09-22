/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120316    K. Kumar          File created.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#include "Tudat/Mathematics/NumericalIntegrators/bulirschStoerVariableStepsizeIntegrator.h"

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

template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::VectorXd, Eigen::VectorXd >;
template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::Vector6d, Eigen::Vector6d >;
template class BulirschStoerVariableStepSizeIntegrator < double, Eigen::MatrixXd, Eigen::MatrixXd >;


} // namespace numerical_integrators

} // namespace tudat
