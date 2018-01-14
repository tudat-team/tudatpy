/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CUSTOM_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_CUSTOM_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Basics/basicTypedefs.h"

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
 *  in the createFlightConditions file, which can be used to define constant coefficients.
 *  NOTE: Functionality of this class is tested in test_aerodynamic_coefficient_generator
 *  test suite.
 */
class CustomAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
{
public:

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
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz) (default true).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction (default true).
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
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true ):
        AerodynamicCoefficientInterface( referenceLength, referenceArea, lateralReferenceLength,
                                         momentReferencePoint, independentVariableNames,
                                         areCoefficientsInAerodynamicFrame,
                                         areCoefficientsInNegativeAxisDirection )
    {
        coefficientFunction_ = boost::bind(
                    &concatenateForceAndMomentCoefficients, forceCoefficientFunction, momentCoefficientFunction, _1 );
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
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz) (default true).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction (default true).
     */
    CustomAerodynamicCoefficientInterface(
            const boost::function< Eigen::Vector6d( const std::vector< double >& ) >
            coefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true ):
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
                        "Error in CustomAerodynamicCoefficientInterface, number of input variables is inconsistent " +
                        std::to_string( independentVariables.size( ) ) + ", " +
                        std::to_string( numberOfIndependentVariables_ ) );
        }

        // Update current coefficients.
        Eigen::Vector6d currentCoefficients = coefficientFunction_(
                    independentVariables );
        currentForceCoefficients_ = currentCoefficients.segment( 0, 3 );
        currentMomentCoefficients_ = currentCoefficients.segment( 3, 3 );
    }

    //! Function to reset the constant aerodynamic coefficients, only valid if coefficients are already constant
    /*!
     * Function to reset the constant aerodynamic coefficients, only valid if coefficients are already constant. Function
     * checks if the numberOfIndependentVariables_ is equal to zero, and throws an error if it is not.
     * \param constantCoefficients New force and moment coefficients (in that order) expressed in the same frame as existing
     * coefficients.
     */
    void resetConstantCoefficients(
            const Eigen::Vector6d& constantCoefficients )
    {
        if( numberOfIndependentVariables_ != 0 )
        {
            throw std::runtime_error( "Error when setting constant aerodynamic coefficients, numberOfIndependentVariables_ is not equal to 0 " );
        }
        coefficientFunction_ = boost::lambda::constant( constantCoefficients );
    }

private:

    //! Function returning the concatenated aerodynamic force and moment coefficients as function of
    //! the set of independent variables.
    boost::function< Eigen::Vector6d( const std::vector< double >& ) >
    coefficientFunction_;


};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_CUSTOM_AERODYNAMIC_COEFFICIENT_INTERFACE_H
