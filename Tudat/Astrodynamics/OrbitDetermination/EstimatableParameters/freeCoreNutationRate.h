/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FREECORENUTATIONRATE_H
#define TUDAT_FREECORENUTATIONRATE_H


#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/fullPlanetaryRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's free core nutation rate (for a full planetary rotational model).
/*!
 *  Interface class for estimation of a body's free core nutation rate (for a full planetary rotational model).
 *  Interfaces the estimation with the free core nutation rate of a PlanetaryRotationModel object
 */
class FreeCoreNutationRate: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel PlanetaryRotationModel object of which the free core nutation rate is a property.
     *  \param associatedBody Name of body of which parameter is a property.
     */
    FreeCoreNutationRate(
            const std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< double >( free_core_nutation_rate, associatedBody ),
        rotationModel_( rotationModel ) { }

    //! Destructor
    ~FreeCoreNutationRate( ) { }


    //! Get value of the free core nutation rate of the body whose rotational ephemeris is described by a full planetary rotational model.
    /*!
     *  Get value of of the free core nutation rate of the body whose rotational ephemeris is described by a full planetary rotational model.
     *  \return Free core nutation rate
     */
    double getParameterValue( )
    {
        return rotationModel_->getPlanetaryOrientationAngleCalculator()->getFreeCoreNutationRate();
    }


    //! Reset value of the free core nutation rate of the body whose rotational ephemeris is described by a full planetary rotational model.
    /*!
     *  Reset value of the free core nutation rate of the body whose rotational ephemeris is described by a full planetary rotational model.
     *  \param parameterValue New value of the free core nutation rate.
     */
    void setParameterValue( const double parameterValue )
    {
        rotationModel_->getPlanetaryOrientationAngleCalculator()->resetFreeCoreNutationRate( parameterValue );
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

    //! PlanetaryRotationModel object of which the free core nutation rate is a property
    std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_FREECORENUTATIONRATE_H
