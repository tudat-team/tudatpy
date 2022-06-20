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

#ifndef TUDAT_COMPOSITE_FUNCTION_SPHERICAL_SHAPING_H
#define TUDAT_COMPOSITE_FUNCTION_SPHERICAL_SHAPING_H

#include "tudat/astro/low_thrust/shape_based/baseFunctionsSphericalShaping.h"
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>

namespace tudat
{
namespace shape_based_methods
{

class CompositeFunctionSphericalShaping
{
public:

    //! Constructor.
    CompositeFunctionSphericalShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeFunctionSphericalShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    std::vector< std::shared_ptr< BaseFunctionSphericalShaping > > getCompositeComponents()
    {
        return compositeFunctionComponents_;
    }

    std::shared_ptr< BaseFunctionSphericalShaping > getCompositeComponent( int indexComponent )
    {
        return compositeFunctionComponents_[ indexComponent ];
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients( const Eigen::VectorXd& compositeFunctionCoefficients );



    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime );

    double getComponentFunctionFirstDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionSecondDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionThirdDerivative(
        const int componentIndex, const double currentTime );


protected:

private:

    std::vector< std::shared_ptr< BaseFunctionSphericalShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};

class CompositeRadialFunctionSphericalShaping : public CompositeFunctionSphericalShaping
{
public:

    //! Constructor.
    CompositeRadialFunctionSphericalShaping( const Eigen::VectorXd& compositeFunctionCoefficients ){

        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

        std::vector< std::shared_ptr< BaseFunctionSphericalShaping > > compositeFunctionComponents;
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( constantSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( linearSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( squaredSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( cosineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( powerCosineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( sineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( powerSineSphericalShaping ) );
        compositeFunctionComponents_ = compositeFunctionComponents;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeRadialFunctionSphericalShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients(
            const Eigen::VectorXd& compositeFunctionCoefficients );

    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime );

    double getComponentFunctionFirstDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionSecondDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionThirdDerivative(
        const int componentIndex, const double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionSphericalShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};


class CompositeElevationFunctionSphericalShaping : public CompositeFunctionSphericalShaping
{
public:

    //! Constructor.
    CompositeElevationFunctionSphericalShaping( const Eigen::VectorXd& compositeFunctionCoefficients ){

        compositeFunctionCoefficients_ = compositeFunctionCoefficients;

        std::vector< std::shared_ptr< BaseFunctionSphericalShaping > > compositeFunctionComponents;
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( cosineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( powerCosineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( sineSphericalShaping ) );
        compositeFunctionComponents.push_back( createBaseFunctionSphericalShaping( powerSineSphericalShaping ) );
        compositeFunctionComponents_ = compositeFunctionComponents;

    }


    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CompositeElevationFunctionSphericalShaping( ) { }

    //! Return coefficients.
    /*!
     * Returns coefficients of composite function.
     */
    Eigen::VectorXd getCompositeFunctionCoefficients()
    {
        return compositeFunctionCoefficients_;
    }

    //! Reset coefficients.
    /*!
     * Resets the coefficients of the composite function.
     */
    void resetCompositeFunctionCoefficients(
            const Eigen::VectorXd& compositeFunctionCoefficients );

    double evaluateCompositeFunction( const double independentVariable );

    double evaluateCompositeFunctionFirstDerivative( const double independentVariable );

    double evaluateCompositeFunctionSecondDerivative( const double independentVariable );

    double evaluateCompositeFunctionThirdDerivative( const double independentVariable );

    double getComponentFunctionCurrentValue(
        const int componentIndex, const double currentTime );

    double getComponentFunctionFirstDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionSecondDerivative(
        const int componentIndex, const double currentTime );

    double getComponentFunctionThirdDerivative(
        const int componentIndex, const double currentTime );


protected:

private:

    std::vector<  std::shared_ptr< BaseFunctionSphericalShaping > > compositeFunctionComponents_;

    //! Vector of coefficients.
    /*!
     * VectorXd containing coefficients corresponding to components of the composite function.
     */
    Eigen::VectorXd compositeFunctionCoefficients_;

};


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_COMPOSITE_FUNCTION_SPHERICAL_SHAPING_H
