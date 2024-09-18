/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H
#define TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H

#include <vector>

#include <functional>


#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"



namespace tudat
{

namespace aerodynamics
{

//! Base class to hold an aerodynamic coefficient interface for a vehicle control surface
/*!
 *  Base class to hold an aerodynamic coefficient interface for a vehicle control surface.
 *  This interface can, for instance, be a database of coefficients or an aerodynamic analysis code  which generates
 *  coefficients.The aerodynamic coefficients are defined as a function of any number of independent
 *  variables, the physical meaning of which is stored in the coefficient interface. The control surface increment interface
 *  is always stored in an AerodynamicCoefficientInterface, and never used independently in propagation. As such, the
 *  aerodynamic reference quantities of this interface are not stored here, but are instead implicitly assumed to be equal
 *  to the AerodynamicCoefficientInterface hosting it.
 */
class ControlSurfaceIncrementAerodynamicInterface
{
public:

    //! Constructor
    /*!
     *  Constructor
     *  \param independentVariableNames List of independent variable names. Note that an error is shown if not exactly one
     *  of these is a control_surface_deflection_dependent.
     */
    ControlSurfaceIncrementAerodynamicInterface(
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        independentVariableNames_( independentVariableNames )
    {
        // Count number of control surface deflections.
        int numberOfControlSurfaceDeflections = 0;
        for( unsigned int i = 0; i < independentVariableNames_.size( ); i++ )
        {
            if( independentVariableNames_.at( i ) == control_surface_deflection_dependent )
            {
                numberOfControlSurfaceDeflections++;
            }
        }

        // Provide warning if number of control surface deflections is not equal to 1.
        if( numberOfControlSurfaceDeflections != 1 )
        {
            std::string errorMessage = "Warnining when making control surface deflection interface, " + std::to_string(
                        numberOfControlSurfaceDeflections ) + " deflection independent variables found, must be 1 per object ";
            throw std::runtime_error( errorMessage );
        }

        numberOfIndependentVariables_ = static_cast< unsigned int >( independentVariableNames_.size( ) );
    }

    //! Destructor
    virtual ~ControlSurfaceIncrementAerodynamicInterface( ){ }

    //! Compute the aerodynamic coefficient increments of the control surface.
    /*!
     *  Computes the current force and moment coefficients increments of the control surface, and is to be
     *  implemented in derived classes. Input is a set of independent variables
     *  (doubles) which represent the variables from which the coefficients are calculated
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     */
    virtual void updateCurrentCoefficients(
            std::vector< double >& independentVariables ) = 0;

    //! Function for returning current aerodynamic force coefficients
    /*!
     *  Function for returning current aerodynamic force coefficients.
     *  \return Force coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentForceCoefficients( )
    {
        return currentForceCoefficients_;
    }

    //! Function for returning current aerodynamic moment coefficients
    /*!
     *  Function for returning current aerodynamic moment coefficients.
     *  \return Moment coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentMomentCoefficients( )
    {
        return currentMomentCoefficients_;
    }

    //! Function for returning current aerodynamic force and moment coefficients
    /*!
     *  Function for returning current aerodynamic force and moment coefficients
     *  \return Force and moment coefficients at given independent variables
     */
    Eigen::Matrix< double, 6, 1 > getCurrentAerodynamicCoefficients(  )
    {
        Eigen::Matrix< double, 6, 1 > coefficients;
        coefficients.segment( 0, 3 ) = getCurrentForceCoefficients( );
        coefficients.segment( 3, 3 ) = getCurrentMomentCoefficients( );
        return coefficients;
    }

    //! Function to return the identifiers of the physical meaning of each independent variable.
    /*!
     *  Function to return the identifiers of the physical meaning of each independent variable
     *  of the aerodynamic coefficient interface.
     *  \return A vector with the identifiers of the physical meaning of each independent variable.
     */
    std::vector< AerodynamicCoefficientsIndependentVariables > getIndependentVariableNames( )
    {
        return independentVariableNames_;
    }

    //! Function to return a single identifier of the physical meaning of one independent variable.
    /*!
     *  Function to return a single identifier of the physical meaning of one of the independent
     *  independent variable of the coefficient interface. The index of the variable is defined
     *  by the input variable.
     *  \param index Index of list of identfiers to return
     *  \return The identifiers of the physical meaning of the independent variable at the position
     *  of the input variable.
     */
    AerodynamicCoefficientsIndependentVariables getIndependentVariableName(
            const unsigned int index )
    {
        if( index >= numberOfIndependentVariables_ )
        {
            throw std::runtime_error(
                        std::string( "Error when retrieving aerodynamic coefficient interface " ) +
                        ( " variable name, requested variable index " ) +
                        std::to_string( index ) +
                        ", but only " + std::to_string(
                            numberOfIndependentVariables_ ) + " variables available." );
        }

        return independentVariableNames_.at( index );
    }

    //! Function to return the number of independent variables upon which the coeficients depend.
    /*!
     *  Function to return the number of independent variables upon which the coeficients depend.
     *  The size of the vector used as input for updateCurrentCoefficients should always have the
     *  size returned by this variable.
     *  \return Number of independent variables upon which the coeficients depend
     */
    unsigned int getNumberOfIndependentVariables( )
    {
        return numberOfIndependentVariables_;
    }
protected:


    //! Vector with identifiers for the physical meaning of each independent variable of the aerodynamic coefficients.
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames_;

    //! Number of independent variables upon which the force and moment coefficient increments depend.
    /*!
     *  Number of independent variables upon which the force and moment coefficient increments depend, i.e.
     *  the length of the vectors that should be used as input to updateCurrentCoefficients
     */
    unsigned int numberOfIndependentVariables_;

    //! The current force coefficient increments.
    /*!
     * The force coefficients increments at the current flight condition.
     */
    Eigen::Vector3d currentForceCoefficients_;

    //! The current moment coefficient increments.
    /*!
     * The moment coefficient increments at the current flight condition.
     */
    Eigen::Vector3d currentMomentCoefficients_;
};


//! Control surface aerodynamic coefficient interface taking function pointers providing aerodynamics
//! coefficient increments as a function of independent variables (doubles).
/*!
 *  Control surface aerodynamic coefficient interface taking function pointers providing aerodynamic
 *  coefficient increments as a function of independent variables (doubles). The origin of the coefficients
 *  or the nature of the independent variables is irrelevant for this class. A typical use of this class is to
 *  interface interpolated tabulated control surface increments.
 */
class CustomControlSurfaceIncrementAerodynamicInterface: public ControlSurfaceIncrementAerodynamicInterface
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param coefficientFunction Function returning the concatenated aerodynamic force and moment
     *  coefficient increments as function of the set of independent variables.
     *  \param independentVariableNames Vector with identifiers for the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     */
    CustomControlSurfaceIncrementAerodynamicInterface(
            const std::function< Eigen::Vector6d( const std::vector< double >& ) > coefficientFunction,
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        ControlSurfaceIncrementAerodynamicInterface( independentVariableNames ),
        coefficientFunction_( coefficientFunction ){ }

    CustomControlSurfaceIncrementAerodynamicInterface(
            const std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction,
            const std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction,
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        ControlSurfaceIncrementAerodynamicInterface( independentVariableNames )
    {
        coefficientFunction_ = std::bind(
                    &concatenateForceAndMomentCoefficients, forceCoefficientFunction, momentCoefficientFunction,
                    std::placeholders::_1 );
    }


    //! Destructor
    ~CustomControlSurfaceIncrementAerodynamicInterface( ){ }

    //! Compute the aerodynamic coefficient increments of the control surface.
    /*!
     *  Computes the current force and moment coefficients increments of the control surface.
     *  Input is a set of independent variables (doubles) which represent the variables from which the coefficients are
     *  calculated
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     */
    void updateCurrentCoefficients( std::vector< double >& independentVariables )
    {
        // Check if the correct number of aerodynamic coefficients is provided.
        if( independentVariables.size( ) != numberOfIndependentVariables_ )
        {
            throw std::runtime_error(
                        "Error in CustomControlSurfaceIncrementAerodynamicInterface, number of "
                        "input variables is inconsistent " );
        }

        // Update current coefficients.
        Eigen::Vector6d currentCoefficients = coefficientFunction_(
                    independentVariables );
        currentForceCoefficients_ = currentCoefficients.segment( 0, 3 );
        currentMomentCoefficients_ = currentCoefficients.segment( 3, 3 );
    }

protected:

    //! Function returning the concatenated aerodynamic force and moment coefficient increments as function of the set of
    //! independent variables.
    std::function< Eigen::Vector6d( const std::vector< double >& ) > coefficientFunction_;
};

}

}

#endif // TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H
