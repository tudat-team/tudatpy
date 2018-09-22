/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PPNPARAMETERS_H
#define TUDAT_PPNPARAMETERS_H

#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class used to estimate PPN parameter gamma
class PPNParameterGamma: public EstimatableParameter< double >
{

public:

    //! Constuctor
    /*!
     * Constuctor
     * \param ppnParameterSet Object used to store PPN parameters
     */
    PPNParameterGamma( const std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet
                       = relativity::ppnParameterSet ):
        EstimatableParameter< double >( ppn_parameter_gamma, "global_metric" ),
      ppnParameterSet_( ppnParameterSet ){ }

    //! Destructor
    ~PPNParameterGamma( ) { }

    //! Function to get the current value of the PPN parameter gamma.
    /*!
     * Function to get the current value of the PPN parameter gamma.
     * \return Current value of the PPN parameter gamma.
     */
    double getParameterValue( )
    {
        return ppnParameterSet_->getParameterGamma( );
    }

    //! Function to reset the value of the PPN parameter gamma.
    /*!
     * Function to reset the value of the PPN parameter gamma.
     * \param parameterValue New value of the PPN parameter gamma.
     */
    void setParameterValue( double parameterValue )
    {
        ppnParameterSet_->setParameterGamma( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Object used to store PPN parameters
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

};

//! Class used to estimate PPN parameter beta
class PPNParameterBeta: public EstimatableParameter< double >
{

public:

    //! Constuctor
    /*!
     * Constuctor
     * \param ppnParameterSet Object used to store PPN parameters
     */
    PPNParameterBeta( const std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet
                      = relativity::ppnParameterSet  ):
        EstimatableParameter< double >( ppn_parameter_beta, "global_metric" ),
      ppnParameterSet_( ppnParameterSet ){ }

    //! Destructor
    ~PPNParameterBeta( ) { }

    //! Function to get the current value of the PPN parameter beta.
    /*!
     * Function to get the current value of the PPN parameter beta.
     * \return Current value of the PPN parameter beta.
     */
    double getParameterValue( )
    {
        return ppnParameterSet_->getParameterBeta( );
    }

    //! Function to reset the value of the PPN parameter beta.
    /*!
     * Function to reset the value of the PPN parameter beta.
     * \param parameterValue New value of the PPN parameter beta.
     */
    void setParameterValue( double parameterValue )
    {
        ppnParameterSet_->setParameterBeta( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Object used to store PPN parameters
    std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

};

}

}

#endif // TUDAT_PPNPARAMETERS_H
