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
#include "tudat/interface/mcd/marsClimateDatabaseClimateModel.h"

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

    McdAtmosphereModel( const std::shared_ptr< mcd_interface::MarsClimateDatabaseClimateModel > marsClimateDatabaseClimateModel ) :
        marsClimateDatabaseClimateModel_( marsClimateDatabaseClimateModel ) 
    {
        requiredExtVar_ = { 
            mcd_interface::ExtVar::ratio_of_specific_heats, 
            mcd_interface::ExtVar::reduced_molecular_gas_constant
        };

        setRequiresClimateModel( );

        marsClimateDatabaseClimateModel_->addExtraVariableKeys( requiredExtVar_ );
        marsClimateDatabaseClimateModel_->setZkey( 3 );
        
    }

    //! Destructor
    virtual ~McdAtmosphereModel( ) {}

    double getDensity( double, double, double, double ) override;

    double getPressure( double, double, double, double ) override;

    double getTemperature( double, double, double, double ) override;

    double getSpeedOfSound( double, double, double, double ) override;

    double getZonalWind( double, double, double, double ) const;

    double getMeridionalWind( double, double, double, double ) const;

protected:

    std::shared_ptr< mcd_interface::MarsClimateDatabaseClimateModel > marsClimateDatabaseClimateModel_;

    std::vector< mcd_interface::ExtVar > requiredExtVar_;
    
};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_MCD_ATMOSPHERE_MODEL_H