/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MCD_ATMOSPHERE_MODEL_H
#define TUDAT_MCD_ATMOSPHERE_MODEL_H

#include <iostream>
#include <vector>
#include <map>
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/interface/mcd/marsClimateDatabase.h"

namespace tudat
{
namespace aerodynamics
{

//! MCD Atmosphere Model class
/*!
 * Class for Mars Climate Database atmosphere model.
 * This class interfaces with the MCD Fortran routines to provide
 * atmospheric properties (density, temperature, pressure, winds) for Mars.
 *
 * ALTITUDE CONVENTION:
 * --------------------
 * Input altitude is "height above local surface" (matching Tudat convention).
 * Internally uses MCD's zkey=3 mode, which allows MCD to handle the conversion
 * to radial distance using its own areoid and topography models.
 *
 * When highResolutionMode=1, MCD uses MOLA topography for accurate surface height.
 * When highResolutionMode=0, MCD uses GCM resolution topography.
 *
 * PARAMETERS:
 * -----------
 * See constructor documentation for detailed parameter descriptions.
 *
 * THREAD SAFETY:
 * --------------
 * This class is NOT thread-safe due to cached internal state.
 * Each thread should use its own instance.
 */
class McdAtmosphereModel : public AtmosphereModel
{
public:

    McdAtmosphereModel( const std::shared_ptr< mcd_interface::MarsClimateDatabase > marsClimateDatabase ) :
        marsClimateDatabase_( marsClimateDatabase ) 
    {
        requiredExtVar_ = { 
            mcd_interface::ExtVar::ratio_of_specific_heats, 
            mcd_interface::ExtVar::reduced_molecular_gas_constant
        };
    }

    //! Destructor
    virtual ~McdAtmosphereModel( ) {}

    virtual double getDensity( );

    virtual double getPressure( );

    virtual double getTemperature( );

    virtual double getSpeedOfSound( );

    double getZonalWind( ) const
    {
        return marsClimateDatabase_->getZonalWind( );
    }

    double getMeridionalWind( ) const
    {
        return marsClimateDatabase_->getMeridionalWind( );
    }

protected:

    std::shared_ptr< mcd_interface::MarsClimateDatabase > marsClimateDatabase_;

    std::vector< mcd_interface::ExtVar > requiredExtVar_;
    
};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_MCD_ATMOSPHERE_MODEL_H