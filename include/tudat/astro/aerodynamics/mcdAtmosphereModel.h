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

        marsClimateDatabase_->addExtraVariableKeys( requiredExtVar_ );
        marsClimateDatabase_->setZkey( 3 );
        
    }

    //! Destructor
    virtual ~McdAtmosphereModel( ) {}

    double getDensity( [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double ) override;

    double getPressure( [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double ) override;

    double getTemperature( [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double ) override;

    double getSpeedOfSound( [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double, [[maybe_unused]] double ) override;

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