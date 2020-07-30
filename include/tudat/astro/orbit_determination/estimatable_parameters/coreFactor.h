/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_COREFACTOR_H
#define TUDAT_COREFACTOR_H


#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's core factor (for a full planetary rotational model).
/*!
 *  Interface class for estimation of a body's core factor (for a full planetary rotational model).
 *  Interfaces the estimation with the core factor of a PlanetaryRotationModel object
 */
class CoreFactor: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel PlanetaryRotationModel object of which the core factor is a property.
     *  \param associatedBody Name of body of which parameter is a property.
     */
    CoreFactor(
            const std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< double >( core_factor, associatedBody ),
        rotationModel_( rotationModel ) { }

    //! Destructor
    ~CoreFactor( ) { }


    //! Get value of the core factor of the body whose rotational ephemeris is described by a full planetary rotational model.
    /*!
     *  Get value of of the core factor of the body whose rotational ephemeris is described by a full planetary rotational model.
     *  \return Core factor
     */
    double getParameterValue( )
    {
        return rotationModel_->getPlanetaryOrientationAngleCalculator( )->getCorefactor( );
    }


    //! Reset value of the core factor of the body whose rotational ephemeris is described by a full planetary rotational model.
    /*!
     *  Reset value of the core factor of the body whose rotational ephemeris is described by a full planetary rotational model.
     *  \param parameterValue New value of the core factor.
     */
    void setParameterValue( const double parameterValue )
    {
        rotationModel_->getPlanetaryOrientationAngleCalculator( )->resetCoreFactor( parameterValue );
    }


    //! Function to retrieve the size of the parameter (always 1)
    /*!
     *  Function to retrieve the size of the parameter (always 1)
     *  \return Size of parameter.
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! PlanetaryRotationModel object of which the core factor is a property
    std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_COREFACTOR_H
