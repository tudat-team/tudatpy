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

#ifndef TUDAT_BASE_FUNCTIONS_SPHERICAL_SHAPING_H
#define TUDAT_BASE_FUNCTIONS_SPHERICAL_SHAPING_H

#include <math.h>
#include <boost/make_shared.hpp>

namespace tudat
{
namespace shape_based_methods
{

enum baseFunctionSphericalShapingType
{
    constantSphericalShaping,
    linearSphericalShaping,
    squaredSphericalShaping,
    cosineSphericalShaping,
    powerCosineSphericalShaping,
    sineSphericalShaping,
    powerSineSphericalShaping
};

//! Base class for shaping functions to be used in the spherical shaping method.
class BaseFunctionSphericalShaping
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    BaseFunctionSphericalShaping( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~BaseFunctionSphericalShaping( ) { }

    virtual double evaluateFunction( double independentVariable ) = 0;

    virtual double evaluateFirstDerivative( double independentVariable ) = 0;

    virtual double evaluateSecondDerivative( double independentVariable ) = 0;

    virtual double evaluateThirdDerivative( double independentVariable ) = 0;

protected:

private:

};


//! Constant function.
class ConstantFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConstantFunctionSphericalShaping( ){}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConstantFunctionSphericalShaping( ) { }

    double evaluateFunction( double independentVariable )
    {
        return 1.0;
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }

protected:

private:

};


//! Linear function.
class LinearFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    LinearFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~LinearFunctionSphericalShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return independentVariable;
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 1.0;
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }


protected:

private:

};


//! x^2 function.
class SquaredFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SquaredFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SquaredFunctionSphericalShaping( ) { }

    double evaluateFunction( double independentVariable ){
        return std::pow( independentVariable, 2.0 );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return 2.0 * independentVariable;
    }

    double evaluateSecondDerivative( double independentVariable ){
        return 2.0;
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return 0.0;
    }


protected:

private:

};


//! cosine function.
class CosineFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    CosineFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CosineFunctionSphericalShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::cos( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return - std::sin( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - std::cos( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return std::sin( independentVariable );
    }

protected:

private:

};


//! Power cosine function.
class PowerCosineFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerCosineFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerCosineFunctionSphericalShaping( ) { }


    double evaluateFunction( double independentVariable ){
        return independentVariable * std::cos( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return std::cos( independentVariable ) - independentVariable * std::sin( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - 2.0 * std::sin( independentVariable ) - independentVariable * std::cos( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - 3.0 * std::cos( independentVariable ) + independentVariable * std::sin( independentVariable );
    }

protected:

private:

};


//! Sine function.
class SineFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SineFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SineFunctionSphericalShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return std::sin( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable )
    {
        return std::cos( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return - std::sin( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - std::cos( independentVariable );
    }

protected:

private:

};


//! Power times sine function.
class PowerSineFunctionSphericalShaping : public BaseFunctionSphericalShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerSineFunctionSphericalShaping( ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerSineFunctionSphericalShaping( ) { }


    double evaluateFunction( double independentVariable )
    {
        return independentVariable * std::sin( independentVariable );
    }

    double evaluateFirstDerivative( double independentVariable ){
        return std::sin( independentVariable ) + independentVariable * std::cos( independentVariable );
    }

    double evaluateSecondDerivative( double independentVariable )
    {
        return 2.0 * std::cos( independentVariable ) - independentVariable * std::sin( independentVariable );
    }

    double evaluateThirdDerivative( double independentVariable )
    {
        return - 3.0 * std::sin( independentVariable ) - independentVariable * std::cos( independentVariable );
    }

protected:

private:

};

std::shared_ptr< BaseFunctionSphericalShaping > createBaseFunctionSphericalShaping(
        const baseFunctionSphericalShapingType baseFunctionType );


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_BASE_FUNCTIONS_SPHERICAL_SHAPING_H
