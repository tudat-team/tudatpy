/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_A1MARSDENPARAMETER_H
#define TUDAT_A1MARSDENPARAMETER_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/basic_astro/marsdenAcceleration.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a constant body drag coefficient
class A1MarsdenParameter : public EstimatableParameter< double >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param coefficientInterface
     * \param associatedBody Body for which the A1 coefficient is considered.
     */
    A1MarsdenParameter( const std::shared_ptr< basic_astrodynamics::MarsdenAcceleration > marsdenAcceleration,
                        const std::string& associatedBody,
                        const std::string& centralBody ):
        EstimatableParameter< double >( A1_parameter_marsden_acceleration, associatedBody, centralBody ),
        marsdenAcceleration_( marsdenAcceleration )
    {
        if( centralBody == "" )
        {
            throw std::runtime_error( "Error when creating A1 Marsden parameter; central body is empty." );
        }
    }

    //! Destructor
    ~A1MarsdenParameter( ) {}

    //! Function to get the current value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to get the current value of the constant drag coefficient that is to be estimated.
     * \return Current value of the constant drag coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return marsdenAcceleration_->getA1( );
    }

    //! Function to reset the value of the constant A1 coefficient that is to be estimated.
    /*!
     * \param parameterValue New value of the constant A1 coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        marsdenAcceleration_->resetA1( parameterValue );
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
    std::shared_ptr< basic_astrodynamics::MarsdenAcceleration > marsdenAcceleration_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_A1MARSDENPARAMETER_H
