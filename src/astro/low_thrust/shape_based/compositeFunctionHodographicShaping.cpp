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


#include <cmath>
#include <iostream>
#include "tudat/astro/low_thrust/shape_based/compositeFunctionHodographicShaping.h"

namespace tudat
{
namespace shape_based_methods
{


void CompositeFunctionHodographicShaping::resetCompositeFunctionCoefficients(
        const Eigen::VectorXd& compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.rows() == static_cast< int >( compositeFunctionComponents_.size( ) ) )
    {
        compositeFunctionCoefficients_ = compositeFunctionCoefficients;
    }
    else
    {
        std::cerr << "The number of composite function coefficients is incorrect!" << std::endl
             << compositeFunctionCoefficients.size() << " coefficients were put in, however, "
             << compositeFunctionComponents_.size() << " are required. The composite function coefficients are not reset!\n";
    }
}


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionCurrentValue (const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionDerivativeCurrentValue( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateDerivative( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionIntegralCurrentValue (const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateIntegral( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::getComponentFunctionDerivativeCurrentValue(
        const int componentIndex, const double independentVariable )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateDerivative(independentVariable );
}


double CompositeFunctionHodographicShaping::getComponentFunctionCurrentValue(
        const int componentIndex, const double independentVariable )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction(independentVariable );
}


double CompositeFunctionHodographicShaping::getComponentFunctionIntegralCurrentValue(
        const int componentIndex, const double independentVariable )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateIntegral(independentVariable );
}

} // namespace shape_based_methods
} // namespace tudat
