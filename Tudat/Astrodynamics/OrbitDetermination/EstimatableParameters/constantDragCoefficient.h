/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTANTDRAGCOEFFICIENT_H
#define TUDAT_CONSTANTDRAGCOEFFICIENT_H

#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a constant body drag coefficient
class ConstantDragCoefficient: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param coefficientInterface Object that contains the aerodynamic coefficients. Constructor checks whether it is
     * consistent with this class, e.g. if the existing aerodynamic coefficients are constant
     * \param associatedBody Body for which the drag coefficient is considered.
     */
    ConstantDragCoefficient(
            const std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface,
            const std::string& associatedBody ):
        EstimatableParameter< double >( constant_drag_coefficient, associatedBody ),
        coefficientInterface_( coefficientInterface )
    {
        if( coefficientInterface->getNumberOfIndependentVariables( ) != 0 )
        {
            throw std::runtime_error( "Error when making ConstantDragCoefficient, coefficient interface is inconsistent" );
        }
    }

    //! Destructor
    ~ConstantDragCoefficient( ) { }

    //! Function to get the current value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to get the current value of the constant drag coefficient that is to be estimated.
     * \return Current value of the constant drag coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return coefficientInterface_->getConstantCoefficients( )( 0 );
    }

    //! Function to reset the value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to reset the value of the constant drag coefficient that is to be estimated.
     * \param parameterValue New value of the constant drag coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        Eigen::Vector6d currentCoefficientSet =
                coefficientInterface_->getCurrentAerodynamicCoefficients( );
        currentCoefficientSet( 0 ) = parameterValue;
        coefficientInterface_->resetConstantCoefficients( currentCoefficientSet );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Object that contains the aerodynamic coefficients
    std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_CONSTANTDRAGCOEFFICIENT_H
