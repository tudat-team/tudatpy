/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_YARKOVSKYPARAMETER_H
#define TUDAT_YARKOVSKYPARAMETER_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a constant body drag coefficient
class YarkovskyParameter: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param coefficientInterface Object that contains the aerodynamic coefficients. Constructor checks whether it is
     * consistent with this class, e.g. if the existing aerodynamic coefficients are constant
     * \param associatedBody Body for which the drag coefficient is considered.
     */
    YarkovskyParameter(
            const std::shared_ptr< electromagnetism::YarkovskyAcceleration > yarkovskyAcceleration,
            const std::string& associatedBody,
            const std::string& centralBody ):
        EstimatableParameter< double >( yarkovsky_parameter, associatedBody, centralBody ),
        yarkovskyAcceleration_( yarkovskyAcceleration )
    {
        if( centralBody == "" )
        {
            throw std::runtime_error( "Error when creating Yarkovsky parameter; central body is empty." );
        }
    }

    //! Destructor
    ~YarkovskyParameter( ) { }

    //! Function to get the current value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to get the current value of the constant drag coefficient that is to be estimated.
     * \return Current value of the constant drag coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return yarkovskyAcceleration_->getYarkovskyParameter( );
    }

    //! Function to reset the value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to reset the value of the constant drag coefficient that is to be estimated.
     * \param parameterValue New value of the constant drag coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        yarkovskyAcceleration_->setYarkovskyParameter( parameterValue );
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
    std::shared_ptr< electromagnetism::YarkovskyAcceleration > yarkovskyAcceleration_;
};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_YARKOVSKYPARAMETER_H
