/*    Copyright (c) 2010-2018, Delft University of Technology
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

namespace tudat
{

namespace relativity
{


//! Minkowski metric (-1,1,1,1 signature) represented as a Matrix
const static Eigen::Matrix4d minkowskiMetric = ( Eigen::Matrix4d( ) <<
                                                 -1.0, 0.0, 0.0, 0.0,
                                                 0.0, 1.0, 0.0, 0.0,
                                                 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 0.0, 1.0 ).finished( );

//! Class that stores the PPN parameters, typically used as a 'global' environment property stored in ppnParameterSet variable
class PPNParameterSet
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterGamma Value of PPN parameter gamma.
     * \param parameterBeta Value of PPN parameter beta.
     */
    PPNParameterSet( const double parameterGamma, const double parameterBeta ):
        parameterGamma_( parameterGamma ), parameterBeta_( parameterBeta )
    { }

    //! Destructor
    ~PPNParameterSet( ){ }

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

protected:

    //! Value of PPN parameter gamma.
    double parameterGamma_;

    //! Value of PPN parameter beta.
    double parameterBeta_;

};

//! Global PPN parameter set, initialized upon compilation (with values equal to GR).
extern std::shared_ptr< PPNParameterSet > ppnParameterSet;

//! Global parameter denoting EP violation in proper time rate, initialized to GR value of 0 upon compilation.
extern double equivalencePrincipleLpiViolationParameter;

}

}

#endif // TUDAT_METRIC_H
