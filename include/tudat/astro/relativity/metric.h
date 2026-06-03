/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_METRIC_H
#define TUDAT_METRIC_H

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/relativity/relativisticPotentials.h"

namespace tudat
{

namespace relativity
{

//! Minkowski metric (-1,1,1,1 signature) represented as a Matrix
const static Eigen::Matrix4d minkowskiMetric =
        ( Eigen::Matrix4d( ) << -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 ).finished( );

//! Class that stores the PPN parameters used by relativity models.
class PPNParameterSet
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param parameterGamma Value of PPN parameter gamma.
     * \param parameterBeta Value of PPN parameter beta.
     */
    PPNParameterSet( const double parameterGamma,
                     const double parameterBeta,
                     const double parameterDelta = 0.0,
                     const double parameterEpsilon = 0.0 ):
        parameterGamma_( parameterGamma ), parameterBeta_( parameterBeta ), parameterDelta_( parameterDelta ),
        parameterEpsilon_( parameterEpsilon )
    {}

    //! Destructor
    ~PPNParameterSet( ) {}

    //! Function to retrieve value of PPN parameter gamma.
    /*!
     * Function to retrieve value of PPN parameter gamma.
     * \return Value of PPN parameter gamma.
     */
    double getParameterGamma( )
    {
        return parameterGamma_;
    }

    //! Function to retrieve value of PPN parameter beta.
    /*!
     * Function to retrieve value of PPN parameter beta.
     * \return Value of PPN parameter beta.
     */
    double getParameterBeta( )
    {
        return parameterBeta_;
    }

    //! Function to retrieve value of PPN parameter epsilon.
    /*!
     * Function to retrieve value of PPN parameter epsilon.
     * \return Value of PPN parameter epsilon.
     */
    double getParameterEpsilon( )
    {
        return parameterEpsilon_;
    }

    //! Function to retrieve value of PPN parameter delta.
    /*!
     * Function to retrieve value of PPN parameter delta.
     * \return Value of PPN parameter delta.
     */
    double getParameterDelta( )
    {
        return parameterDelta_;
    }

    //! Function to reset value of PPN parameter gamma.
    /*!
     * Function to reset value of PPN parameter gamma.
     * \param parameterGamma New value of PPN parameter gamma.
     */
    void setParameterGamma( const double parameterGamma )
    {
        parameterGamma_ = parameterGamma;
    }

    //! Function to reset value of PPN parameter beta.
    /*!
     * Function to reset value of PPN parameter beta.
     * \param parameterBeta New value of PPN parameter beta.
     */
    void setParameterBeta( const double parameterBeta )
    {
        parameterBeta_ = parameterBeta;
    }

    //! Function to reset value of PPN parameter epsilon.
    /*!
     * Function to reset value of PPN parameter epsilon.
     * \param parameterEpsilon New value of PPN parameter epsilon.
     */
    void setParameterEpsilon( const double parameterEpsilon )
    {
        parameterEpsilon_ = parameterEpsilon;
    }

    //! Function to reset value of PPN parameter delta.
    /*!
     * Function to reset value of PPN parameter delta.
     * \param parameterDelta New value of PPN parameter delta.
     */
    void setParameterDelta( const double parameterDelta )
    {
        parameterDelta_ = parameterDelta;
    }

protected:
    //! Value of PPN parameter gamma.
    double parameterGamma_;

    //! Value of PPN parameter beta.
    double parameterBeta_;

    double parameterDelta_;

    double parameterEpsilon_;
};

class Metric
{
public:
    Metric( )
    {
        currentChristoffelSymbols_.resize( 4 );
    }

    virtual ~Metric( ) {}

    virtual std::shared_ptr< Metric > Clone( ) = 0;

    Eigen::Matrix< double, 4, 4 > getCurrentCovariantMetric( )
    {
        return minkowskiMetric + getCurrentCovariantMetricPeturbation( );
    }

    Eigen::Matrix< double, 4, 4 > getCurrentContravariantMetric( )
    {
        return getCurrentCovariantMetric( ).inverse( );
    }

    Eigen::Matrix< double, 4, 4 > getCurrentCovariantMetricPeturbation( )
    {
        return currentCovariantMetricContribution_;
    }

    Eigen::Matrix< double, 4, 4 > getCurrentContravariantMetricPeturbation( )
    {
        return -( getCurrentContravariantMetric( ) -
                  minkowskiMetric );  // minkowskiMetric * currentCovariantMetricContribution_ * minkowskiMetric;
    }

    std::vector< Eigen::Matrix< double, 4, 4 > > getCurrentChristoffelSymbols( )
    {
        return currentChristoffelSymbols_;
    }

    Eigen::Matrix< double, 6, 1 > getCurrentEvaluationState( )
    {
        return currentEvaluationState_;
    }

    double getCurrentTime( )
    {
        return currentTime_;
    }

    virtual void update( const Eigen::Matrix< double, 6, 1 >& state,
                         const double time,
                         const bool updateCurrentMetric,
                         const bool updateCurrentChristoffelSymbols ) = 0;

protected:
    Eigen::Matrix< double, 4, 4 > currentCovariantMetricContribution_;

    std::vector< Eigen::Matrix< double, 4, 4 > > currentChristoffelSymbols_;

    Eigen::Matrix< double, 6, 1 > currentEvaluationState_;

    double currentTime_;
};

}  // namespace relativity

}  // namespace tudat

#endif  // TUDAT_METRIC_H
