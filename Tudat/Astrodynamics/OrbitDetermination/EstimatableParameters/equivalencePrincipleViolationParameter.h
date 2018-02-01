/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EQUIVALENCEPRINCIPLEVIOLATIONPARAMETER_H
#define TUDAT_EQUIVALENCEPRINCIPLEVIOLATIONPARAMETER_H

#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class used to estimate equivalence principle LPI violation parameter, e.g. putative effect on observers proper time rate
class EquivalencePrincipleLpiViolationParameter: public EstimatableParameter< double >
{

public:

    //! Constuctor
    /*!
     * Constuctor
     */
    EquivalencePrincipleLpiViolationParameter( ):
        EstimatableParameter< double >( equivalence_principle_lpi_violation_parameter, "global_metric" ){ }

    //! Destructor
    ~EquivalencePrincipleLpiViolationParameter( ) { }

    //! Function to get the current value of the equivalence principle LPI violation parametera.
    /*!
     * Function to get the current value of the equivalence principle LPI violation parametera.
     * \return Current value of the equivalence principle LPI violation parametera.
     */
    double getParameterValue( )
    {
        return relativity::equivalencePrincipleLpiViolationParameter;
    }

    //! Function to reset the value of the equivalence principle LPI violation parametera.
    /*!
     * Function to reset the value of the equivalence principle LPI violation parametera.
     * \param parameterValue New value of the equivalence principle LPI violation parametera.
     */
    void setParameterValue( double parameterValue )
    {
        relativity::equivalencePrincipleLpiViolationParameter = parameterValue;
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

};

}

}

#endif // TUDAT_EQUIVALENCEPRINCIPLEVIOLATIONPARAMETER_H
