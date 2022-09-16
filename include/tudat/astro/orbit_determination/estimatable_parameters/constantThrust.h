/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTANTTHRUST_H
#define TUDAT_CONSTANTTHRUST_H

#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{


class ConstantThrustMagnitudeParameter: public EstimatableParameter< double >
{

public:


    ConstantThrustMagnitudeParameter(
            const std::shared_ptr< propulsion::ConstantThrustMagnitudeWrapper > thrustWrapper,
            const std::string& associatedBody,
            const std::string& engineId ):
        EstimatableParameter< double  >( constant_thrust_magnitude_parameter, associatedBody, engineId ),
        thrustWrapper_( thrustWrapper ) { }

    ~ConstantThrustMagnitudeParameter( ) { }


    double getParameterValue( )
    {
        return thrustWrapper_->getCurrentThrustForceMagnitude( );
    }

    void setParameterValue( const double parameterValue )
    {
        thrustWrapper_->resetConstantThrustForceMagnitude( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    const std::shared_ptr< propulsion::ConstantThrustMagnitudeWrapper > thrustWrapper_;
};

template< typename SpecificImpulseSource >
class ConstantSpecificImpulseParameter: public EstimatableParameter< double >
{

public:


    ConstantSpecificImpulseParameter(
            const std::shared_ptr< SpecificImpulseSource > thrustWrapper,
            const std::string& associatedBody,
            const std::string& engineId ):
        EstimatableParameter< double  >( constant_specific_impulse, associatedBody, engineId ),
        thrustWrapper_( thrustWrapper ) { }

    ~ConstantSpecificImpulseParameter( ) { }


    double getParameterValue( )
    {
        return thrustWrapper_->getCurrentSpecificImpulse( );
    }

    void setParameterValue( const double parameterValue )
    {
        thrustWrapper_->resetConstantSpecificImpulse( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    const std::shared_ptr< SpecificImpulseSource > thrustWrapper_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_CONSTANTTHRUST_H
