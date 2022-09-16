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

#ifndef TUDAT_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H
#define TUDAT_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H

#include <vector>
#include <cmath>
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
    BaseFunctionHodographicShaping(  ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~BaseFunctionHodographicShaping( ) { }

    virtual double evaluateFunction( const double independentVariable ) = 0;

    virtual double evaluateDerivative( const double independentVariable ) = 0;

    virtual double evaluateIntegral( const double independentVariable ) = 0;

protected:

private:

};


typedef std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > HodographicBasisFunctionList;


//! Constant function.
class ConstantFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Default constructor.
    ConstantFunctionHodographicShaping( ){}

    //! Default destructor.
    ~ConstantFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable )
    {
        return 1.0;
    }

    double evaluateDerivative( const double independentVariable )
    {
        return 0.0;
    }

    double evaluateIntegral( const double independentVariable )
    {
        return independentVariable;
    }

protected:

private:

};


//! Sine function.
class SineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the frequency of the sine function as input.
    SineFunctionHodographicShaping( double frequency )
    {
        frequency_ = frequency;
    }

    //! Default destructor.
    ~SineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double frequency_;
};


//! Cosine function.
class CosineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the frequency of the cosine function as input.
    CosineFunctionHodographicShaping( double frequency )
    {
        frequency_ = frequency;
    }

    //! Default destructor.
    ~CosineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double frequency_;
};


//! Exponential function.
class ExponentialFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the order of the exponential function as input.
    ExponentialFunctionHodographicShaping( double exponent )
    {
        exponent_ = exponent;
    }

    //! Default destructor.
    ~ExponentialFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponent_;

};

//! Scaled exponential function.
class ScaledExponentialFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the exponent and scaling factor as inputs.
    ScaledExponentialFunctionHodographicShaping( double exponent, double scaleFactor )
    {
        exponent_ = exponent;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledExponentialFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponent_;
    double scaleFactor_;
};

//! Exponential times sine function.
class ExponentialSineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking exponent of the exponential function and frequency of the sine function as inputs.
    ExponentialSineFunctionHodographicShaping( double exponentExponentialFunction, double frequencySineFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencySineFunction_ = frequencySineFunction;
    }

    //! Default destructor.
    ~ExponentialSineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentExponentialFunction_;
    double frequencySineFunction_;
};

//! Scaled exponential times sine function.
class ScaledExponentialSineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking exponent of the exponential function, frequency of the sine function, and scaling factor as inputs.
    ScaledExponentialSineFunctionHodographicShaping( double exponentExponentialFunction, double frequencySineFunction, double scaleFactor )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencySineFunction_ = frequencySineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledExponentialSineFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentExponentialFunction_;
    double frequencySineFunction_;
    double scaleFactor_;
};

//! Exponential times cosine function.
class ExponentialCosineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking exponent of the exponential function and frequency of the cosine function as inputs.
    ExponentialCosineFunctionHodographicShaping( double exponentExponentialFunction, double frequencyCosineFunction )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    //! Default destructor.
    ~ExponentialCosineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentExponentialFunction_;
    double frequencyCosineFunction_;
};


//! Scaled exponential times cosine function.
class ScaledExponentialCosineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking order of the exponential function, frequency of the cosine function, and scaling factor as inputs.
    ScaledExponentialCosineFunctionHodographicShaping( double exponentExponentialFunction, double frequencyCosineFunction, double scaleFactor )
    {
        exponentExponentialFunction_ = exponentExponentialFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledExponentialCosineFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentExponentialFunction_;
    double frequencyCosineFunction_;
    double scaleFactor_;
};


//! Power function.
class PowerFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor with the order of the power function as input.
    PowerFunctionHodographicShaping( double exponent )
    {
        exponent_ = exponent;
    }

    //! Default destructor.
    ~PowerFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponent_;

};


//! Scaled power function.
class ScaledPowerFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking the exponent and scaling factor as inputs.
    ScaledPowerFunctionHodographicShaping( double exponent, double scaleFactor )
    {
        exponent_ = exponent;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledPowerFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponent_;
    double scaleFactor_;
};


//! Power times sine function.
class PowerSineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking exponent of the power function and frequency of the sine function as inputs.
    PowerSineFunctionHodographicShaping( double exponentPowerFunction, double frequencySineFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencySineFunction_ = frequencySineFunction;
    }

    //! Default destructor.
    ~PowerSineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentPowerFunction_;
    double frequencySineFunction_;
};


class ScaledPowerSineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    ScaledPowerSineFunctionHodographicShaping( double exponentPowerFunction, double frequencySineFunction, double scaleFactor )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencySineFunction_ = frequencySineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledPowerSineFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentPowerFunction_;
    double frequencySineFunction_;
    double scaleFactor_;
};


//! Power times cosine function.
class PowerCosineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    //! Constructor taking exponent of the power function and frequency of the cosine function as inputs.
    PowerCosineFunctionHodographicShaping( double exponentPowerFunction, double frequencyCosineFunction )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
    }

    //! Default destructor.
    ~PowerCosineFunctionHodographicShaping( ) { }


    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentPowerFunction_;
    double frequencyCosineFunction_;
};


class ScaledPowerCosineFunctionHodographicShaping : public BaseFunctionHodographicShaping
{
public:

    ScaledPowerCosineFunctionHodographicShaping( double exponentPowerFunction, double frequencyCosineFunction, double scaleFactor )
    {
        exponentPowerFunction_ = exponentPowerFunction;
        frequencyCosineFunction_ = frequencyCosineFunction;
        scaleFactor_ = scaleFactor;
    }

    //! Default destructor.
    ~ScaledPowerCosineFunctionHodographicShaping( ) { }

    double evaluateFunction( const double independentVariable );

    double evaluateDerivative( const double independentVariable );

    double evaluateIntegral( const double independentVariable );


protected:

private:
    double exponentPowerFunction_;
    double frequencyCosineFunction_;
    double scaleFactor_;
};

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographConstant( )
{
    return std::make_shared< ConstantFunctionHodographicShaping >( );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographSine(  double frequency )
{
    return std::make_shared< SineFunctionHodographicShaping >( frequency );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographCosine(  double frequency )
{
    return std::make_shared< CosineFunctionHodographicShaping >( frequency );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographExponential(  double exponent )
{
    return std::make_shared< ExponentialFunctionHodographicShaping >( exponent );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographScaledExponential(  double exponent, double scaleFactor )
{
    return std::make_shared< ScaledExponentialFunctionHodographicShaping >( exponent, scaleFactor );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographExponentialSine(  double exponentExponentialFunction, double frequencySineFunction )
{
    return std::make_shared< ExponentialSineFunctionHodographicShaping >( exponentExponentialFunction, frequencySineFunction );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographScaledExponentialSine(  double exponentExponentialFunction, double frequencySineFunction, double scaleFactor )
{
    return std::make_shared< ScaledExponentialSineFunctionHodographicShaping >( exponentExponentialFunction, frequencySineFunction, scaleFactor );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographExponentialCosine(  double exponentExponentialFunction, double frequencySineFunction )
{
    return std::make_shared< ExponentialCosineFunctionHodographicShaping >( exponentExponentialFunction, frequencySineFunction );
}

inline std::shared_ptr< BaseFunctionHodographicShaping >  hodographScaledExponentialCosine(  double exponentExponentialFunction, double frequencyCosineFunction, double scaleFactor )
{
    return std::make_shared< ScaledExponentialCosineFunctionHodographicShaping >( exponentExponentialFunction, frequencyCosineFunction, scaleFactor );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographPower( double exponent )
{
    return std::make_shared< PowerFunctionHodographicShaping >( exponent );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographScaledPower( double exponent, double scaleFactor )
{
    return std::make_shared< ScaledPowerFunctionHodographicShaping >( exponent, scaleFactor );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographPowerSine(
        double exponentPowerFunction, double frequencySineFunction )
{
    return std::make_shared< PowerSineFunctionHodographicShaping >( exponentPowerFunction, frequencySineFunction );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographScaledPowerSine(
        double exponentPowerFunction, double frequencySineFunction, double scaleFactor )
{
    return std::make_shared< ScaledPowerSineFunctionHodographicShaping >( exponentPowerFunction, frequencySineFunction, scaleFactor );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographPowerCosine(
        double exponentPowerFunction, double frequencyCosineFunction )
{
    return std::make_shared< PowerCosineFunctionHodographicShaping >( exponentPowerFunction, frequencyCosineFunction );
}

inline std::shared_ptr< BaseFunctionHodographicShaping > hodographScaledPowerCosine(
        double exponentPowerFunction, double frequencyCosineFunction, double scaleFactor )
{
    return std::make_shared< ScaledPowerCosineFunctionHodographicShaping >( exponentPowerFunction, frequencyCosineFunction, scaleFactor );
}

} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H
