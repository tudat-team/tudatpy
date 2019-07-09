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

#ifndef BASEFUNCTIONSHODOGRAPHICSHAPING_H
#define BASEFUNCTIONSHODOGRAPHICSHAPING_H

#include <math.h>
#include <boost/make_shared.hpp>

namespace tudat
{
namespace shape_based_methods
{

enum baseFunctionHodographicShapingType
{
    constant,
    sine,
    cosine,
    exponential,
    scaledExponential,
    exponentialSine,
    scaledExponentialSine,
    exponentialCosine,
    scaledExponentialCosine,
    power,
    scaledPower,
    powerCosine,
    scaledPowerCosine,
    powerSine,
    scaledPowerSine
};

//! Base class for shaping functions to be used in the hodographic shaping method.
class BaseFunctionHodographicShaping
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    BaseFunctionHodographicShaping(  ){ }//:
//    power_( 0.0 ),
//    frequency_( 0.0 ){ }

//    BaseFunctionHodographicShaping( double power ):
//    power_( power ),
//    frequency_( 0.0 ){}

//    BaseFunctionHodographicShaping( double power, double frequency ):
//    power_( power ),
//    frequency_( frequency ){}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~BaseFunctionHodographicShaping( ) { }

    virtual double evaluateFunction( double independentVariable ) = 0;

    virtual double evaluateDerivative( double independentVariable ) = 0;

    virtual double evaluateIntegral( double independentVariable ) = 0;

//    double getPowerBaseFunctionHodographicShaping(){
//        return power_;
//    }

//    double getFrequencyBaseFunctionHodographicShaping(){
//        return frequency_;
//    }

protected:

private:

//    double power_;
//    double frequency_;

};


//! Constant function.
class ConstantFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConstantFunction( ){}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConstantFunction( ) { }

    double evaluateFunction( double independentVariable )
    {
        return 1.0;
    }

    double evaluateDerivative( double independentVariable )
    {
        return 0.0;
    }

    double evaluateIntegral( double independentVariable )
    {
        return independentVariable;
    }

protected:

private:

};


//! Sine function.
class SineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    SineFunction( ) : frequency_( 0.0 ) { }

    //! Constructor taking the frequency of the sine function as input.
    SineFunction( double frequency )
    {
        frequency_ = frequency;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setFrequency( double frequency )
    {
        frequency_ = frequency;
    }

protected:

private:
    double frequency_;
};


//! Cosine function.
class CosineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    CosineFunction( ) : frequency_( 0.0 ) { }

    //! Constructor taking the frequency of the cosine function as input.
    CosineFunction( double frequency )
    {
        frequency_ = frequency;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CosineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setFrequency( double frequency )
    {
        frequency_ = frequency;
    }

protected:

private:
    double frequency_;
};


//! Exponential function.
class ExponentialFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialFunction( ) : exponent_( 0.0 ) { }

    //! Constructor taking the order of the exponential function as input.
    ExponentialFunction( double exponent )
    {
        exponent_ = exponent;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ExponentialFunction( ) { }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrder( double exponent )
    {
        exponent_ = exponent;
    }

protected:

private:
    double exponent_;

};

//! Scaled exponential function.
class ScaledExponentialFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ScaledExponentialFunction( ) : exponent_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Constructor taking the exponent and scaling factor as inputs.
    ScaledExponentialFunction( double exponent, double scaleFactor )
    {
        exponent_ = exponent;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledExponentialFunction( ) { }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrder( double exponent )
    {
        exponent_ = exponent;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponent_;
    double scaleFactor_;
};

//! Exponential times sine function.
class ExponentialSineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialSineFunction( ) : exponentExponentialFunction_( 0.0 ), frequencySineFunction_( 0.0 ) { }

    //! Constructor taking exponent of the exponential function and frequency of the sine function as inputs.
    ExponentialSineFunction( double exponentExponentialFunction, double frequencySineFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencySineFunction_ = frequencySineFunction;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ExponentialSineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderExponential( double exponentExponentialFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
    }

    void setFrequency( double frequencySineFunction )
    {
        frequencySineFunction_ = frequencySineFunction;
    }

protected:

private:
    double exponentExponentialFunction_;
    double frequencySineFunction_;
};

//! Scaled exponential times sine function.
class ScaledExponentialSineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ScaledExponentialSineFunction( ) : exponentExponentialFunction_( 0.0 ), frequencySineFunction_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Constructor taking exponent of the exponential function, frequency of the sine function, and scaling factor as inputs.
    ScaledExponentialSineFunction( double exponentExponentialFunction, double frequencySineFunction, double scaleFactor )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencySineFunction_ = frequencySineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledExponentialSineFunction( ) { }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderExponential( double exponentExponentialFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
    }

    void setFrequency( double frequencySineFunction )
    {
        frequencySineFunction_ = frequencySineFunction;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponentExponentialFunction_;
    double frequencySineFunction_;
    double scaleFactor_;
};

//! Exponential times cosine function.
class ExponentialCosineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ExponentialCosineFunction( ) : exponentExponentialFunction_( 0.0 ), frequencyCosineFunction_( 0.0 ) { }

    //! Constructor taking exponent of the exponential function and frequency of the cosine function as inputs.
    ExponentialCosineFunction( double exponentExponentialFunction, double frequencyCosineFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ExponentialCosineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderExponential( double exponentExponentialFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
    }

    void setFrequency( double frequencyCosineFunction )
    {
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

protected:

private:
    double exponentExponentialFunction_;
    double frequencyCosineFunction_;
};


//! Scaled exponential times cosine function.
class ScaledExponentialCosineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ScaledExponentialCosineFunction( ) : exponentExponentialFunction_( 0.0 ), frequencyCosineFunction_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Constructor taking order of the exponential function, frequency of the cosine function, and scaling factor as inputs.
    ScaledExponentialCosineFunction( double exponentExponentialFunction, double frequencyCosineFunction, double scaleFactor )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledExponentialCosineFunction( ) { }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderExponential( double exponentExponentialFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
    }

    void setFrequency( double frequencyCosineFunction )
    {
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponentExponentialFunction_;
    double frequencyCosineFunction_;
    double scaleFactor_;
};


//! Power function.
class PowerFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerFunction( ) : exponent_( 0.0 ) { }

    //! Constructor with the order of the power function as input.
    PowerFunction( double exponent )
    {
        exponent_ = exponent;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrder( double exponent )
    {
        exponent_ = exponent;
    }

protected:

private:
    double exponent_;

};


//! Scaled power function.
class ScaledPowerFunction : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the exponent and scaling factor as inputs.
    ScaledPowerFunction( ) : exponent_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledPowerFunction( ) { }

    ScaledPowerFunction( double exponent, double scaleFactor )
    {
        exponent_ = exponent;
        scaleFactor_ = scaleFactor;
    }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrder( double exponent )
    {
        exponent_ = exponent;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponent_;
    double scaleFactor_;
};


//! Power times sine function.
class PowerSineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerSineFunction( ) : exponentPowerFunction_( 0.0 ), frequencySineFunction_( 0.0 ) { }

    //! Constructor taking exponent of the power function and frequency of the sine function as inputs.
    PowerSineFunction( double exponentPowerFunction, double frequencySineFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencySineFunction_ = frequencySineFunction;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerSineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setExponentPowerFunction( double exponentPowerFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
    }

    void setFrequencySineFunction( double frequency )
    {
        frequencySineFunction_ = frequency;
    }

protected:

private:
    double exponentPowerFunction_;
    double frequencySineFunction_;
};


class ScaledPowerSineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ScaledPowerSineFunction( ) : exponentPowerFunction_( 0.0 ), frequencySineFunction_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledPowerSineFunction( ) { }

    ScaledPowerSineFunction( double exponentPowerFunction, double frequencySineFunction, double scaleFactor )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencySineFunction_ = frequencySineFunction;
        scaleFactor_ = scaleFactor;
    }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderPower( double exponentPowerFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
    }

    void setFrequency( double frequencySineFunction )
    {
        frequencySineFunction_ = frequencySineFunction;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponentPowerFunction_;
    double frequencySineFunction_;
    double scaleFactor_;
};


//! Power times cosine function.
class PowerCosineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    PowerCosineFunction( ) : exponentPowerFunction_( 0.0 ), frequencyCosineFunction_( 0.0 ) { }

    //! Constructor taking exponent of the power function and frequency of the cosine function as inputs.
    PowerCosineFunction( double exponentPowerFunction, double frequencyCosineFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~PowerCosineFunction( ) { }


    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderPower( double exponentPowerFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
    }

    void setFrequency( double frequencyCosineFunction )
    {
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

protected:

private:
    double exponentPowerFunction_;
    double frequencyCosineFunction_;
};


class ScaledPowerCosineFunction : public BaseFunctionHodographicShaping
{
public:
    //! Default constructor.
    /*!
     * Default constructor.
     */
    ScaledPowerCosineFunction( ) : exponentPowerFunction_( 0.0 ), frequencyCosineFunction_( 0.0 ), scaleFactor_( 0.0 ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ScaledPowerCosineFunction( ) { }

    ScaledPowerCosineFunction( double exponentPowerFunction, double frequencyCosineFunction, double scaleFactor )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
        scaleFactor_ = scaleFactor;
    }

    double evaluateFunction( double independentVariable );

    double evaluateDerivative( double independentVariable );

    double evaluateIntegral( double independentVariable );

    void setOrderPower( double exponentPowerFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
    }

    void setFrequency( double frequencyCosineFunction )
    {
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    void setScaleFactor( double scaleFactor )
    {
        scaleFactor_ = scaleFactor;
    }

protected:

private:
    double exponentPowerFunction_;
    double frequencyCosineFunction_;
    double scaleFactor_;
};


} // namespace shape_based_methods
} // namespace tudat

#endif // BASEFUNCTIONSHODOGRAPHICSHAPING_H
