/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATECLIMATEMODEL_H
#define TUDAT_CREATECLIMATEMODEL_H

#include <string>
#include <map>
#include <memory>

#include "tudat/astro/basic_astro/climateModel.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

enum ClimateModelTypes {

#if TUDAT_BUILD_WITH_MCD_INTERFACE

    mars_climate_database,

#endif

};

class ClimateModelSettings
{
public:
    ClimateModelSettings( ClimateModelTypes climateModelType ): climateModelType_( climateModelType ) {}

    virtual ~ClimateModelSettings( ) = default;

    ClimateModelTypes climateModelType_;
};

#if TUDAT_BUILD_WITH_MCD_INTERFACE
class MarsClimateDatabaseClimateModelSettings : public ClimateModelSettings
{
public:
    explicit MarsClimateDatabaseClimateModelSettings( const std::string& mcdDataPath = "",
                                                      const int dustScenario = 1,
                                                      const int perturbationKey = 0,
                                                      const double perturbationSeed = 0.0,
                                                      const double gravityWaveLength = 0.0,
                                                      const int highResolutionMode = 0,
                                                      const int maximumCacheSize = 1000 ):
        ClimateModelSettings( mars_climate_database ), mcdDataPath_( mcdDataPath ), dustScenario_( dustScenario ),
        perturbationKey_( perturbationKey ), perturbationSeed_( perturbationSeed ), gravityWaveLength_( gravityWaveLength ),
        highResolutionMode_( highResolutionMode ), maximumCacheSize_( maximumCacheSize )
    {}

    //! Path to MCD data files
    std::string mcdDataPath_;

    //! Dust and solar EUV scenario (1-8 or 24-35)
    int dustScenario_;

    //! Perturbation type
    int perturbationKey_;

    //! Perturbation seed
    double perturbationSeed_;

    //! Gravity wave wavelength
    double gravityWaveLength_;

    //! High resolution topography flag (0 or 1)
    int highResolutionMode_;

    //! Maximum number of MCD query results retained in the cache
    int maximumCacheSize_;
};

inline std::shared_ptr< ClimateModelSettings > marsClimateDatabaseClimateModelSettings( const std::string& mcdDataPath = "",
                                                                                        const int dustScenario = 1,
                                                                                        const int perturbationKey = 0,
                                                                                        const double perturbationSeed = 0.0,
                                                                                        const double gravityWaveLength = 0.0,
                                                                                        const int highResolutionMode = 0,
                                                                                        const int maximumCacheSize = 1000 )
{
    return std::make_shared< MarsClimateDatabaseClimateModelSettings >(
            mcdDataPath, dustScenario, perturbationKey, perturbationSeed, gravityWaveLength, highResolutionMode, maximumCacheSize );
}

#endif

std::shared_ptr< environment::ClimateModel > createClimateModel( std::shared_ptr< ClimateModelSettings > ClimateModelSettings,
                                                                 std::shared_ptr< simulation_setup::Body > body );

}  // namespace simulation_setup

}  // namespace tudat

#endif
