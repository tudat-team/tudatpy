/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONTINUOUSPROBABILITYDISTRIBUTIONS_H
#define TUDAT_CONTINUOUSPROBABILITYDISTRIBUTIONS_H

#include <Eigen/Core>

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
 *  Base class for a continuous random variable, providing interfaces for evaluating probability and cumulative distribution
 *  functions (inverse cdf is only available for derived class InvertibleContinuousProbabilityDistribution).
 *  Specific functions are to be implemented by derived class. Derived class may be used for random number generation
 *  through ContinuousRandomVariableGenerator class.
 */
template< typename IndependentVariableType >
class ContinuousProbabilityDistribution
{
public:

    //! Destructor
    virtual ~ContinuousProbabilityDistribution( ){ }

    //! Function to evaluate pdf of distribution
    /*!
     *  Function to evaluate probability distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated pdf
     */
    virtual double evaluatePdf( const IndependentVariableType& independentVariable ) = 0;

    //! Function to evaluate cdf of distribution
    /*!
     *  Function to evaluate cumulative distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated cdf
     */
    virtual double evaluateCdf( const IndependentVariableType& independentVariable ) = 0;
};

//! Derived class of ContinuousProbabilityDistribution that includes inverse cdf computation.
template< typename IndependentVariableType >
class InvertibleContinuousProbabilityDistribution: public ContinuousProbabilityDistribution< IndependentVariableType >
{
public:

    using ContinuousProbabilityDistribution< IndependentVariableType >::evaluatePdf;
    using ContinuousProbabilityDistribution< IndependentVariableType >::evaluateCdf;

    //! Destructor
    virtual ~InvertibleContinuousProbabilityDistribution( ){ }

    //! Function to evaluate inverse cdf of distribution
    /*!
     *  Function to evaluate inverse cumulative distribution function at given probability value
     *  \param independentVariable Value of probability at which inverse cdf is to be computed (must be in the domain [0,1]).
     *  \return Evaluated inverse cdf
     */
    virtual double evaluateInverseCdf( const IndependentVariableType independentVariable ) = 0;
};

} // namespace statistics

} // namespace tudat

#endif // TUDAT_CONTINUOUSPROBABILITYDISTRIBUTIONS_H
