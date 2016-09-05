#ifndef CONTINUOUSRANDOMVARIABLE_H
#define CONTINUOUSRANDOMVARIABLE_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"


namespace tudat
{

namespace statistics
{

//! Enum defining implemented continuous boost distributions
enum ContinuousBoostStatisticalDistributions
{
    uniform_boost_distribution = 0,
    normal_boost_distribution = 1,
    exponential_boost_distribution = 2,
    gamma_boost_distribution = 3,
    lognormal_boost_distribution = 4,
    beta_boost_distribution = 5
};

//! Function to evaluate pdf of Gaussian distribution.
/*!
 *  Function to evaluate probability distribution function of Gaussian distribution
 *  \param independentVariable Value of indpendent variable at which pdf it is to be evaluated
 *  \param mean Mean of Gaussian distribution
 *  \param standardDeviation Standard deviation of Gaussian distribution
 *  \return Value of given pdf at requested independent variable value
 */
double evaluateGaussianPdf( const double independentVariable, const double mean, const double standardDeviation );

//! Function to evaluate cdf of Gaussian distribution.
/*!
 *  Function to evaluate cumulative distribution function of Gaussian distribution
 *  \param independentVariable Value of indpendent variable at which cdf it is to be evaluated
 *  \param mean Mean of Gaussian distribution
 *  \param standardDeviation Standard deviation of Gaussian distribution
 *  \return Value of given cdf at requested independent variable value
 */
double calculateGaussianCdf( const double independentVariable, const double mean, const double standardDeviation );

//! Base class for a continuous random variable
/*!
 *  Base class for a continuous random variable, providing interfaces for evaluating probability, cumulative and inverse cumulative
 *  distribution functions. Specific functions are to be implemented by derived class. Derived class may be used for random number
 *  generation through ContinuousVariableClassRandomVariableGenerator class.
 */
template< typename IndependentVariableType >
class ContinuousRandomVariable
{
public:

    //! Destructor
    virtual ~ContinuousRandomVariable( ){ }

    //! Function to evaluate pdf of distribution
    /*!
     *  Function to evaluate probability distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated pdf
     */
    virtual double evaluatePdf( const IndependentVariableType independentVariable ) = 0;

    //! Function to evaluate cdf of distribution
    /*!
     *  Function to evaluate cumulative distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated cdf
     */
    virtual double evaluateCdf( const IndependentVariableType independentVariable ) = 0;

    //! Function to evaluate inverse cdf of distribution
    /*!
     *  Function to evaluate inverse cumulative distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated inverse cdf
     */
    virtual double evaluateInverseCdf( const IndependentVariableType independentVariable ) = 0;
};

}

}

#endif // CONTINUOUSRANDOMVARIABLE_H
