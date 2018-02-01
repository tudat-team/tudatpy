/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITATIONALPARAMETER_H
#define TUDAT_GRAVITATIONALPARAMETER_H

#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a gravitational parameter
class GravitationalParameter: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param gravityFieldModel Gravity field object containing the gravitational parameter to be estimated.
     * \param associatedBody Name of body containing the gravityFieldModel object
     */
    GravitationalParameter(
            const boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel, const std::string& associatedBody ):
        EstimatableParameter< double >( gravitational_parameter, associatedBody ),
        gravityFieldModel_( gravityFieldModel ){ }

    //! Destructor
    ~GravitationalParameter( ) { }

    //! Function to get the current value of the gravitational parameter that is to be estimated.
    /*!
     * Function to get the current value of the gravitational parameter that is to be estimated.
     * \return Current value of the gravitational parameter that is to be estimated.
     */
    double getParameterValue( )
    {
        return gravityFieldModel_->getGravitationalParameter( );
    }

    //! Function to reset the value of the gravitational parameter that is to be estimated.
    /*!
     * Function to reset the value of the gravitational parameter that is to be estimated.
     * \param parameterValue New value of the gravitational parameter that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        gravityFieldModel_->resetGravitationalParameter( parameterValue );
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

    //! Gravity field object containing the gravitational parameter to be estimated.
    boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_GRAVITATIONALPARAMETER_H
