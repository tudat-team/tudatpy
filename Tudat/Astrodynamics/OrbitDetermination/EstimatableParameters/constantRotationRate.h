/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTANTROTATIONRATE_H
#define TUDAT_CONSTANTROTATIONRATE_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's constant rotation rate parameter.
/*!
 * Interface class for estimation of a body's constant rotation rate parameter. Interfaces the estimation with the rotation
 * rate parameter of a SimpleRotationalEphemeris object
 */
class RotationRate: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel SimpleRotationalEphemeris object of which rotation rate parameter is a property
     *  \param associatedBody Name of body of which parameter is a property.
     */
    RotationRate( const std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel,
                  const std::string& associatedBody ):
        EstimatableParameter< double >( constant_rotation_rate, associatedBody ),
        rotationModel_( rotationModel ){ }

    //! Destructor
    ~RotationRate( ) { }

    //! Get value of rotation rate.
    /*!
     *  Get value of rotation rate.
     *  \return Value of rotation rate
     */
    double getParameterValue( )
    {
        return rotationModel_->getRotationRate( );
    }

    //! Reset value of rotation rate.
    /*!
     *  Reset value of rotation rate.
     *  \param parameterValue New value of rotation rate
     */
    void setParameterValue( double parameterValue )
    {
        rotationModel_->resetRotationRate( parameterValue );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, 1 for this parameter
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! SimpleRotationalEphemeris object of which rotation rate parameter is a property
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_CONSTANTROTATIONRATE_H
