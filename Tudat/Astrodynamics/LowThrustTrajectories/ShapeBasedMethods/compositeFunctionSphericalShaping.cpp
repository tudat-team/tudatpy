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


#include "compositeFunctionSphericalShaping.h"
#include <cmath>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{

void CompositeRadialFunctionSphericalShaping::resetCompositeFunctionCoefficients(
        const Eigen::VectorXd& compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.size() == static_cast< int >( compositeFunctionComponents_.size( ) ) )
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


double CompositeRadialFunctionSphericalShaping::evaluateCompositeFunction( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    return 1.0 / functionValue;
}


double CompositeRadialFunctionSphericalShaping::evaluateCompositeFunctionFirstDerivative( const double independentVariable )
{
    double sumComponentsFirstDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsFirstDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }
    return - sumComponentsFirstDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 );
}


double CompositeRadialFunctionSphericalShaping::evaluateCompositeFunctionSecondDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    double sumComponentsSecondDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsSecondDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }
    return - sumComponentsSecondDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 )
            + 2.0 * functionValue * std::pow( evaluateCompositeFunctionFirstDerivative( independentVariable ), 2.0 );
}

double CompositeRadialFunctionSphericalShaping::evaluateCompositeFunctionThirdDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }

    double sumComponentsFirstDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsFirstDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }

    double sumComponentsSecondDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsSecondDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }

    double sumComponentsThirdDerivativeValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        sumComponentsThirdDerivativeValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateThirdDerivative( independentVariable );
    }

    return - sumComponentsThirdDerivativeValue * std::pow( evaluateCompositeFunction( independentVariable ), 2.0 )
            - 2.0 * evaluateCompositeFunction( independentVariable ) * evaluateCompositeFunctionFirstDerivative( independentVariable )
            * sumComponentsSecondDerivativeValue
            + 2.0 * sumComponentsFirstDerivativeValue * std::pow( evaluateCompositeFunctionFirstDerivative( independentVariable ), 2.0 )
            + 4.0 * functionValue * evaluateCompositeFunctionFirstDerivative( independentVariable )
            * evaluateCompositeFunctionSecondDerivative( independentVariable );
}

double CompositeRadialFunctionSphericalShaping::getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction( currentTime );
}


double CompositeRadialFunctionSphericalShaping::getComponentFunctionFirstDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFirstDerivative( currentTime );
}

double CompositeRadialFunctionSphericalShaping::getComponentFunctionSecondDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateSecondDerivative( currentTime );
}

double CompositeRadialFunctionSphericalShaping::getComponentFunctionThirdDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateThirdDerivative( currentTime );
}


void CompositeElevationFunctionSphericalShaping::resetCompositeFunctionCoefficients(
        const Eigen::VectorXd& compositeFunctionCoefficients )
{
    // Check whether the size is correct.
    if( compositeFunctionCoefficients.size() == static_cast< int >( compositeFunctionComponents_.size( ) ) )
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


double CompositeElevationFunctionSphericalShaping::evaluateCompositeFunction( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFunction( independentVariable );
    }
    return functionValue;
}


double CompositeElevationFunctionSphericalShaping::evaluateCompositeFunctionFirstDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateFirstDerivative( independentVariable );
    }
    return functionValue;
}


double CompositeElevationFunctionSphericalShaping::evaluateCompositeFunctionSecondDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateSecondDerivative( independentVariable );
    }
    return functionValue;
}

double CompositeElevationFunctionSphericalShaping::evaluateCompositeFunctionThirdDerivative( const double independentVariable )
{
    double functionValue = 0.0;
    for( unsigned int i = 0; i < compositeFunctionComponents_.size(); i++ )
    {
        functionValue += compositeFunctionCoefficients_[i] * compositeFunctionComponents_[ i ]->evaluateThirdDerivative( independentVariable );
    }
    return functionValue;
}

double CompositeElevationFunctionSphericalShaping::getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFunction( currentTime );
}


double CompositeElevationFunctionSphericalShaping::getComponentFunctionFirstDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateFirstDerivative( currentTime );
}

double CompositeElevationFunctionSphericalShaping::getComponentFunctionSecondDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateSecondDerivative( currentTime );
}

double CompositeElevationFunctionSphericalShaping::getComponentFunctionThirdDerivative(
        const int componentIndex, const double currentTime )
{
    return compositeFunctionComponents_[ componentIndex ]->evaluateThirdDerivative( currentTime );
}


} // namespace shape_based_methods
} // namespace tudat
