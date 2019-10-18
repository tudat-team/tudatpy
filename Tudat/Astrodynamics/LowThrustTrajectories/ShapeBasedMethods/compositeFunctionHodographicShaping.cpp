/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "compositeFunctionHodographicShaping.h"
#include <math.h>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{


void CompositeFunctionHodographicShaping::resetCompositeFunctionCoefficients(
        const Eigen::VectorXd& compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.size() == compositeFunctionComponents_.size() )
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


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionCurrentTime( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionDerivativeCurrentTime( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateDerivative( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::evaluateCompositeFunctionIntegralCurrentTime( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateIntegral( independentVariable );
    }
    return functionValue;
}


double CompositeFunctionHodographicShaping::getComponentFunctionDerivativeCurrentTime(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateDerivative( currentTime );
}


double CompositeFunctionHodographicShaping::getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction( currentTime );
}


double CompositeFunctionHodographicShaping::getComponentFunctionIntegralCurrentTime(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateIntegral( currentTime );
}

} // namespace shape_based_methods
} // namespace tudat
