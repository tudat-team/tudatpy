/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LONGITUDELIBRATIONAMPLITUDE_H
#define TUDAT_LONGITUDELIBRATIONAMPLITUDE_H

#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{


class ScaledLongitudeLibrationAmplitude: public EstimatableParameter< double >
{

public:


    ScaledLongitudeLibrationAmplitude(
            const std::shared_ptr< ephemerides::DirectLongitudeLibrationCalculator > librationCalculator,
            const std::string& associatedBody ):
        EstimatableParameter< double  >( scaled_longitude_libration_amplitude, associatedBody ),
        librationCalculator_( librationCalculator ) { }

    ~ScaledLongitudeLibrationAmplitude( ) { }


    double getParameterValue( )
    {
        return librationCalculator_->getScaledLibrationAmplitude( );
    }

    void setParameterValue( const double parameterValue )
    {
        librationCalculator_->setScaledLibrationAmplitude( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    std::shared_ptr< ephemerides::DirectLongitudeLibrationCalculator > librationCalculator_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_LONGITUDELIBRATIONAMPLITUDE_H
