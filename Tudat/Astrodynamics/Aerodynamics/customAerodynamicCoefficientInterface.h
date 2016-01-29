/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      150416    D. Dirkx          File created.
 *
 */

#ifndef TUDAT_CUSTOM_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_CUSTOM_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace aerodynamics
{

//! Aerodynamic coefficient interface taking function pointers providing aerodynamics
//! coefficients as a function of independent variables (doubles).
/*!
 *  Aerodynamic coefficient interface taking function pointers providing aerodynamics
 *  coefficients as a function of independent variables (doubles). The origin of the coefficients
 *  or the nature of the independent variables is irrelevant for this class.
 *  A factory functios (createConstantCoefficientAerodynamicCoefficientInterface) is provided
 *  at the end of this file, which can be used to define constant coefficients.
 *  NOTE: Functionality of this class is tested in test_aerodynamic_coefficient_generator
 *  test suite.
 */
class CustomAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
{
public:

    //! Function to combined the force and moment coefficients from separate function pointers.
    /*!
     *  Function to combined the force and moment coefficients from separate function pointers.
     *  The output is the concatenated force and moment coefficient vector, evaluated
     *  at the current set of independent variables.
     *  \param forceCoefficientFunction Function returning the aerodynamic force coefficients as
     *  function of the set of independent variables.
     *  \param momentCoefficientFunction Function returning the aerodynamic force coefficients as
     *  function of the set of independent variables.
     *  \param independentVariables Current list of values of the independent variables upon
     *  which the coefficients depend.
     */
    basic_mathematics::Vector6d concatenateForceAndMomentCoefficients(
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >&
            forceCoefficientFunction,
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >&
            momentCoefficientFunction,
            const std::vector< double >& independentVariables )
    {
        return ( basic_mathematics::Vector6d( )<<forceCoefficientFunction( independentVariables ),
                 momentCoefficientFunction( independentVariables ) ).finished( );
    }

    //! Constructor.
    /*!
     *  Constructor.
     *  \param forceCoefficientFunction Function returning the aerodynamic force coefficients as
     *  function of the set of independent variables.
     *  \param momentCoefficientFunction Function returning the aerodynamic force coefficients as
     *  function of the set of independent variables.
     *  \param referenceLength Reference length with which aerodynamic moments
     *  (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
     *  is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. which aerodynamic moment is calculated.
     *  \param independentVariableNames Vector with identifiers for the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    CustomAerodynamicCoefficientInterface(
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >
            forceCoefficientFunction,
            const boost::function< Eigen::Vector3d( const std::vector< double >& ) >
            momentCoefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        AerodynamicCoefficientInterface( referenceLength, referenceArea, lateralReferenceLength,
                                         momentReferencePoint, independentVariableNames,
                                         areCoefficientsInAerodynamicFrame,
                                         areCoefficientsInNegativeAxisDirection )
    {
        coefficientFunction_ = boost::bind(
                    &CustomAerodynamicCoefficientInterface::concatenateForceAndMomentCoefficients,
                    this, forceCoefficientFunction, momentCoefficientFunction, _1 );
    }

    //! Constructor.
    /*!
     *  Constructor.

     *  \param coefficientFunction Function returning the concatenated aerodynamic force and moment
     *  coefficients as function of the set of independent variables.
     *  \param referenceLength Reference length with which aerodynamic moments
     *  (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
     *  is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. which aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers for the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    CustomAerodynamicCoefficientInterface(
            const boost::function< basic_mathematics::Vector6d( const std::vector< double >& ) >
            coefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        AerodynamicCoefficientInterface( referenceLength, referenceArea, lateralReferenceLength,
                                         momentReferencePoint, independentVariableNames,
                                         areCoefficientsInAerodynamicFrame,
                                         areCoefficientsInNegativeAxisDirection ),
        coefficientFunction_( coefficientFunction ){ }

    //! Compute the aerodynamic coefficients at current flight condition.
    /*!
     *  Compute the aerodynamic coefficients at current flight conditions (independent variables).
     *  Input is a set of independent variables (doubles) which represent the variables from which
     *  the coefficients are calculated. The physical nature of these variables depends on
     *  the coefficientFunction_ variables. The size of the independent variable vector must be
     *  numberOfIndependentVariables_
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     */
    virtual void updateCurrentCoefficients( const std::vector< double >& independentVariables )
    {
        // Check if the correct number of aerodynamic coefficients is provided.
        if( independentVariables.size( ) != numberOfIndependentVariables_ )
        {
            throw std::runtime_error(
                        "Error in CustomAerodynamicCoefficientInterface, number of "
                        "input variables is inconsistent " );
        }

        // Update current coefficients.
        basic_mathematics::Vector6d currentCoefficients = coefficientFunction_(
                    independentVariables );
        currentForceCoefficients_ = currentCoefficients.segment( 0, 3 );
        currentMomentCoefficients_ = currentCoefficients.segment( 3, 3 );
    }

private:

    //! Function returning the concatenated aerodynamic force and moment coefficients as function of
    //! the set of independent variables.
    boost::function< basic_mathematics::Vector6d( const std::vector< double >& ) >
    coefficientFunction_;


};

//! Function to create an aerodynamic coefficient interface containing constant coefficients.
/*!
 *  Function to create an aerodynamic coefficient interface containing constant coefficients,
 *  As a result, the generated coefficient interface depends on zero parameters.
 *  \param constantForceCoefficient Constant force coefficients.
 *  \param constantMomentCoefficient Constant moment coefficients.
 *  \param referenceLength Reference length with which aerodynamic moments
 *  (about x- and z- axes) are non-dimensionalized.
 *  \param referenceArea Reference area with which aerodynamic forces and moments are
 *  non-dimensionalized.
 *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
 *  is non-dimensionalized.
 *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
 *  coefficients are typically defined in negative direction.
 *  \return Aerodynamic coefficient interface with constant coefficients.
 */
boost::shared_ptr< AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = 0,
        const bool areCoefficientsInNegativeAxisDirection = 1 );

template< int NumberOfDimensions >
boost::shared_ptr< AerodynamicCoefficientInterface >
createTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::vector< std::vector< double > > independentVariables,
        const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients,
        const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
        independentVariableNames,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = 0,
        const bool areCoefficientsInNegativeAxisDirection = 1 )
{
    if( independentVariables.size( ) != NumberOfDimensions )
    {
        std::cerr<<"Error when creating tabulated aerodynamic coefficient interface, "
                <<"inconsistent variable vector dimensioning"<<std::endl;
    }

    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
        std::cerr<<"Error when creating tabulated aerodynamic coefficient interface, "
                <<"inconsistent variable name vector dimensioning"<<std::endl;

    }

    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > forceInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, forceCoefficients );

    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > momentInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, momentCoefficients );

    return  boost::make_shared< CustomAerodynamicCoefficientInterface >(
                boost::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             forceInterpolator, _1 ),
                boost::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             momentInterpolator, _1 ),
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}

}

}

#endif // CUSTOM_TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H
