/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_IAUROTATIONMODELPARAMETERS_H
#define TUDAT_IAUROTATIONMODELPARAMETERS_H


#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/ephemerides/iauRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

class NominalRotationPoleParameter: public EstimatableParameter< Eigen::VectorXd >
{
public:

    NominalRotationPoleParameter(
            const std::shared_ptr< ephemerides::IauRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< Eigen::VectorXd >( nominal_rotation_pole_position, associatedBody ),
        rotationModel_( rotationModel ){ }

    //! Destructor
    ~NominalRotationPoleParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        return rotationModel_->getNominalPole( );
    }

    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        if( parameterValue.rows( ) != 2 )
        {
            throw std::runtime_error( "Error when resetting nominal pole position, new parameter size is incompatible" );
        }
        rotationModel_->setNominalPole( parameterValue );

    }

    int getParameterSize( )
    {
        return 2;
    }

protected:

private:

    std::shared_ptr< ephemerides::IauRotationModel > rotationModel_;
};


class RotationPoleRateParameter: public EstimatableParameter< Eigen::VectorXd >
{
public:

    RotationPoleRateParameter(
            const std::shared_ptr< ephemerides::IauRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< Eigen::VectorXd >( rotation_pole_position_rate, associatedBody ),
        rotationModel_( rotationModel ){ }

    //! Destructor
    ~RotationPoleRateParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        return rotationModel_->getPolePrecession( );
    }

    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        if( parameterValue.rows( ) != 2 )
        {
            throw std::runtime_error( "Error when resetting nominal pole position, new parameter size is incompatible" );
        }
        rotationModel_->setPolePrecession( parameterValue );

    }

    int getParameterSize( )
    {
        return 2;
    }

protected:

private:

    std::shared_ptr< ephemerides::IauRotationModel > rotationModel_;
};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_IAUROTATIONMODELPARAMETERS_H
