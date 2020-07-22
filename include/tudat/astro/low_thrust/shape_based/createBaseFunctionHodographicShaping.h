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

#ifndef TUDAT_CREATE_BASE_FUNCTION_HODOGRAPHIC_SHAPING_H
#define TUDAT_CREATE_BASE_FUNCTION_HODOGRAPHIC_SHAPING_H

#include <cmath>
#include <boost/make_shared.hpp>
#include "tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h"

namespace tudat
{
namespace shape_based_methods
{

//! Settings class for hodographic shaping base functions.
class BaseFunctionHodographicShapingSettings
{
public:

    //! Default, empty constructor.
    BaseFunctionHodographicShapingSettings( ){ }

    //! Destructor.
    virtual ~BaseFunctionHodographicShapingSettings( ){ }

};


//! Settings class for hodographic shaping trigonometric function (either sine or cosine).
class TrigonometricFunctionHodographicShapingSettings: public BaseFunctionHodographicShapingSettings
{
public:

    //! Constructor to set frequency for a trigonometric function (either sine or cosine).
    /*!
     *  Constructor to set frequency for a trigonometric function (either sine or cosine).
     *  \param frequency Frequency of the trigonometric function.
     */
    TrigonometricFunctionHodographicShapingSettings( const double frequency ) :
        frequency_( frequency ){ }

    //! Frequency of the trigonometric function.
    double frequency_;
};


//! Settings class for hodographic shaping exponential function.
class ExponentialFunctionHodographicShapingSettings: public BaseFunctionHodographicShapingSettings
{
public:

    //! Constructor to set exponent, and possibly scaling factor for the exponential function.
    /*!
     *  Constructor to set exponent, and possibly scaling factor for the exponential function.
     *  \param exponent Exponent of the exponential function.
     *  \param scaleFactor Scaling factor for the exponential function.
     */
    ExponentialFunctionHodographicShapingSettings( const double exponent,
                                                   const double scaleFactor = 1.0 ) :
        exponent_( exponent ), scaleFactor_( scaleFactor ){ }

    //! Exponent of the exponential function.
    double exponent_;

    //! Scale factor, to be used for a scaled exponential function.
    double scaleFactor_;
};


//! Settings class for hodographic shaping exponential times trigonometric function.
class ExponentialTimesTrigonometricFunctionHodographicShapingSettings: public BaseFunctionHodographicShapingSettings
{
public:

    //! Constructor to set exponent, frequency, and possibly scaling factor for the exponential times sine or cosine function.
    /*!
     *  Constructor to set exponent, frequency, and possibly scaling factor for the exponential times sine or cosine function.
     *  \param exponent Exponent of the exponential function.
     *  \param frequency Frequency of the trigonometric function.
     *  \param scaleFactor Scaling factor for the exponential times trigonometric function.
     */
    ExponentialTimesTrigonometricFunctionHodographicShapingSettings( const double exponent,
                                                                     const double frequency,
                                                                     const double scaleFactor = 1.0 ) :
        exponent_( exponent ), frequency_( frequency ), scaleFactor_( scaleFactor ){ }

    //! Exponent of the exponential function.
    double exponent_;

    //! Frequency of the trigonometric function.
    double frequency_;

    //! Scale factor, to be used for a scaled exponential times trigonometric function.
    double scaleFactor_;
};


//! Settings class for hodographic shaping power function.
class PowerFunctionHodographicShapingSettings: public BaseFunctionHodographicShapingSettings
{
public:

    //! Constructor to set exponent, and possibly scaling factor for the power function.
    /*!
     *  Constructor to set exponent, and possibly scaling factor for the power function.
     *  \param exponent Exponent of the power function.
     *  \param scaleFactor Scaling factor for the power function.
     */
    PowerFunctionHodographicShapingSettings( const double exponent,
                                             const double scaleFactor = 1.0 ) :
        exponent_( exponent ), scaleFactor_( scaleFactor ){ }

    //! Exponent of the power function.
    double exponent_;

    //! Scale factor, to be used for a scaled power function.
    double scaleFactor_;
};


//! Settings class for hodographic shaping power times trigonometric function.
class PowerTimesTrigonometricFunctionHodographicShapingSettings: public BaseFunctionHodographicShapingSettings
{
public:

    //! Constructor to set exponent, frequency, and possibly scaling factor for the power times sine or cosine function.
    /*!
     *  Constructor to set exponent, frequency, and possibly scaling factor for the power times sine or cosine function.
     *  \param exponent Exponent of the power function.
     *  \param frequency Frequency of the trigonometric function.
     *  \param scaleFactor Scaling factor for the power times trigonometric function.
     */
    PowerTimesTrigonometricFunctionHodographicShapingSettings( const double exponent,
                                                               const double frequency,
                                                               const double scaleFactor = 1.0 ) :
        exponent_( exponent ), frequency_( frequency ), scaleFactor_( scaleFactor ){ }

    //! Exponent of the power function.
    double exponent_;

    //! Frequency of the trigonometric function.
    double frequency_;

    //! Scale factor, to be used for a scaled power times trigonometric function.
    double scaleFactor_;
};


//! Function to create a base function for hodographic shaping.
std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > createBaseFunctionHodographicShaping(
        const shape_based_methods::baseFunctionHodographicShapingType baseFunctionType,
        const std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > baseFunctionSettings );


} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_CREATE_BASE_FUNCTION_HODOGRAPHIC_SHAPING_H
