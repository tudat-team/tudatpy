
#ifndef TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H
#define TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H

#include <vector>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
#include <Tudat/Astrodynamics/Aerodynamics/aerodynamics.h>

namespace tudat
{

namespace aerodynamics
{

class ControlSurfaceIncrementAerodynamicInterface
{
public:

    ControlSurfaceIncrementAerodynamicInterface(
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        independentVariableNames_( independentVariableNames )
    {
        int numberOfControlSurfaceDeflections = 0;
        for( unsigned int i = 0; i < independentVariableNames_.size( ); i++ )
        {
            if( independentVariableNames_.at( i ) == control_surface_deflection_dependent )
            {
                numberOfControlSurfaceDeflections++;
            }
        }

        if( numberOfControlSurfaceDeflections != 1 )
        {
            std::string errorMessage = "Error when making control surface deflection interface, " + boost::lexical_cast< std::string >(
                        numberOfControlSurfaceDeflections ) + " deflection independent variables found, must be 1 per object ";
            throw std::runtime_error( errorMessage );
        }

        numberOfIndependentVariables_ = independentVariableNames_.size( );
    }

    virtual ~ControlSurfaceIncrementAerodynamicInterface( ){ }

    virtual void updateCurrentCoefficients(
            std::vector< double >& independentVariables ) = 0;

    //! Function for turning aerodynamic force coefficients
    /*!
     *  Function for returning aerodynamic force coefficients.
     *  \return Force coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentForceCoefficients( )
    {
        return currentForceCoefficients_;
    }

    //! Function for returning aerodynamic moment coefficients
    /*!
     *  Function for returning aerodynamic moment coefficients.
     *  \return Moment coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentMomentCoefficients( )
    {
        return currentMomentCoefficients_;
    }

    //! Function for returning aerodynamic force and moment coefficients
    /*!
     *  Function for returning aerodynamic force and moment coefficients
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
                        boost::lexical_cast< std::string >( index ) +
                        ", but only " + boost::lexical_cast< std::string >(
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
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames_;

    unsigned int numberOfIndependentVariables_;

    Eigen::Vector3d currentForceCoefficients_;

    Eigen::Vector3d currentMomentCoefficients_;
};

class CustomControlSurfaceIncrementAerodynamicInterface: public ControlSurfaceIncrementAerodynamicInterface
{
public:
    CustomControlSurfaceIncrementAerodynamicInterface(
            const boost::function< basic_mathematics::Vector6d( const std::vector< double >& ) > coefficientFunction,
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        ControlSurfaceIncrementAerodynamicInterface( independentVariableNames ),
        coefficientFunction_( coefficientFunction )
    {

    }

    ~CustomControlSurfaceIncrementAerodynamicInterface( ){ }

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
        basic_mathematics::Vector6d currentCoefficients = coefficientFunction_(
                    independentVariables );
        currentForceCoefficients_ = currentCoefficients.segment( 0, 3 );
        currentMomentCoefficients_ = currentCoefficients.segment( 3, 3 );
    }
protected:

    boost::function< basic_mathematics::Vector6d( const std::vector< double >& ) > coefficientFunction_;
};

}

}

#endif // TUDAT_CONTROLSURFACEAERODYNAMICCOEFFICIENTINTERFACE_H
