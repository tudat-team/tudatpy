/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEBODYDEFORMATIONMODEL_H
#define TUDAT_CREATEBODYDEFORMATIONMODEL_H

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/ground_stations/oceanPoleTideDeformation.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

class SystemOfBodies;

std::shared_ptr< basic_astrodynamics::Iers2010EarthDeformation > createDefaultEarthIers2010DeformationModel(
        const std::shared_ptr< ephemerides::Ephemeris > earthEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > lunarEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > solarEphemeris,
        const std::shared_ptr< ephemerides::RotationalEphemeris > earthRotation,
        const std::function< double( ) > gravitionalParametersOfEarth,
        const std::function< double( ) > gravitionalParametersOfMoon,
        const std::function< double( ) > gravitionalParametersOfSun,
        const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction = {} );

class BodyDeformationSettings
{
public:
    BodyDeformationSettings( const gravitation::BodyDeformationTypes bodyDeformationType ): bodyDeformationType_( bodyDeformationType ) {}

    virtual ~BodyDeformationSettings( ) {}

    gravitation::BodyDeformationTypes getBodyDeformationType( )
    {
        return bodyDeformationType_;
    }

protected:
    gravitation::BodyDeformationTypes bodyDeformationType_;
};

class BasicSolidBodyDeformationSettings : public BodyDeformationSettings
{
public:
    BasicSolidBodyDeformationSettings( const std::vector< std::string > deformingBodies,
                                       const std::map< int, std::pair< double, double > > displacementLoveNumbers,
                                       const double bodyReferenceRadius = TUDAT_NAN ):
        BodyDeformationSettings( gravitation::basic_solid_body ), deformingBodies_( deformingBodies ),
        displacementLoveNumbers_( displacementLoveNumbers ), bodyReferenceRadius_( bodyReferenceRadius )
    {}

    std::vector< std::string > getDeformingBodies( )
    {
        return deformingBodies_;
    }
    std::map< int, std::pair< double, double > > getDisplacementLoveNumbers( )
    {
        return displacementLoveNumbers_;
    }
    double getBodyReferenceRadius( )
    {
        return bodyReferenceRadius_;
    }

protected:
    std::vector< std::string > deformingBodies_;
    std::map< int, std::pair< double, double > > displacementLoveNumbers_;
    double bodyReferenceRadius_;
};

//! Preloaded or plain-text ocean pole-load grid; loading is performed once when the model is constructed.
enum class AtmosphericTideInputFormat { grid, stationCoefficients };
class AtmosphericTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    AtmosphericTideBodyDeformationSettings( const input_output::AtmosphericTideGrid& grid,
                                            const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), grid_( grid ), settings_( settings ),
        format_( AtmosphericTideInputFormat::grid ), preloaded_( true )
    {}
    AtmosphericTideBodyDeformationSettings( const input_output::AtmosphericTideStationData& stations,
                                            const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), stations_( stations ), settings_( settings ),
        format_( AtmosphericTideInputFormat::stationCoefficients ), preloaded_( true )
    {}
    AtmosphericTideBodyDeformationSettings( const std::string& path,
                                            const AtmosphericTideInputFormat format,
                                            const input_output::AtmosphericTideMetadata& metadata,
                                            const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), file_( path ), metadata_( metadata ), settings_( settings ),
        format_( format ), preloaded_( false )
    {}
    const input_output::AtmosphericTideGrid& getGrid( ) const
    {
        return grid_;
    }
    const input_output::AtmosphericTideStationData& getStations( ) const
    {
        return stations_;
    }
    const input_output::AtmosphericTideMetadata& getMetadata( ) const
    {
        return metadata_;
    }
    const basic_astrodynamics::AtmosphericTideEvaluationSettings& getEvaluationSettings( ) const
    {
        return settings_;
    }
    const std::string& getFile( ) const
    {
        return file_;
    }
    AtmosphericTideInputFormat getInputFormat( ) const
    {
        return format_;
    }
    bool hasPreloadedData( ) const
    {
        return preloaded_;
    }

private:
    input_output::AtmosphericTideGrid grid_;
    input_output::AtmosphericTideStationData stations_;
    std::string file_;
    input_output::AtmosphericTideMetadata metadata_;
    basic_astrodynamics::AtmosphericTideEvaluationSettings settings_;
    AtmosphericTideInputFormat format_;
    bool preloaded_;
};
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGrid(
        const input_output::AtmosphericTideGrid& grid,
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( grid, settings );
}
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStations(
        const input_output::AtmosphericTideStationData& stations,
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( stations, settings );
}
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGridFile(
        const std::string& path,
        const input_output::AtmosphericTideMetadata& metadata = {},
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( path, AtmosphericTideInputFormat::grid, metadata, settings );
}
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStationFile(
        const std::string& path,
        const input_output::AtmosphericTideMetadata& metadata = {},
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >(
            path, AtmosphericTideInputFormat::stationCoefficients, metadata, settings );
}

class OceanPoleTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    OceanPoleTideBodyDeformationSettings( const input_output::OceanPoleTideGrid& grid,
                                          const basic_astrodynamics::OceanPoleTideLatitudePolicy policy ):
        BodyDeformationSettings( gravitation::ocean_pole_tide ), grid_( grid ), policy_( policy ), preloaded_( true )
    {}
    OceanPoleTideBodyDeformationSettings( const std::string& gridFile, const basic_astrodynamics::OceanPoleTideLatitudePolicy policy ):
        BodyDeformationSettings( gravitation::ocean_pole_tide ), gridFile_( gridFile ), policy_( policy ), preloaded_( false )
    {}
    const input_output::OceanPoleTideGrid& getGrid( ) const
    {
        return grid_;
    }
    const std::string& getGridFile( ) const
    {
        return gridFile_;
    }
    basic_astrodynamics::OceanPoleTideLatitudePolicy getLatitudePolicy( ) const
    {
        return policy_;
    }
    bool hasPreloadedGrid( ) const
    {
        return preloaded_;
    }

private:
    input_output::OceanPoleTideGrid grid_;
    std::string gridFile_;
    basic_astrodynamics::OceanPoleTideLatitudePolicy policy_;
    bool preloaded_;
};

inline std::shared_ptr< OceanPoleTideBodyDeformationSettings > oceanPoleTideBodyShapeDeformation(
        const std::string& gridFile,
        const basic_astrodynamics::OceanPoleTideLatitudePolicy policy =
                basic_astrodynamics::OceanPoleTideLatitudePolicy::requireGridCoverage )
{
    return std::make_shared< OceanPoleTideBodyDeformationSettings >( gridFile, policy );
}
inline std::shared_ptr< OceanPoleTideBodyDeformationSettings > oceanPoleTideBodyShapeDeformationFromData(
        const input_output::OceanPoleTideGrid& grid,
        const basic_astrodynamics::OceanPoleTideLatitudePolicy policy =
                basic_astrodynamics::OceanPoleTideLatitudePolicy::requireGridCoverage )
{
    return std::make_shared< OceanPoleTideBodyDeformationSettings >( grid, policy );
}

class OceanTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    OceanTideBodyDeformationSettings(
            const std::vector< std::string > blqFiles,
            const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp ):
        BodyDeformationSettings( gravitation::ocean_tide ), blqFiles_( blqFiles ), synthesis_( synthesis )
    {}

    OceanTideBodyDeformationSettings(
            const input_output::OceanLoadingData& data,
            const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp ):
        BodyDeformationSettings( gravitation::ocean_tide ), preloadedData_( data ), usePreloadedData_( true ), synthesis_( synthesis )
    {
        input_output::validateOceanLoadingData( preloadedData_ );
    }

    std::vector< std::string > getBlqFiles( )
    {
        return blqFiles_;
    }
    const input_output::OceanLoadingData& getOceanLoadingData( ) const
    {
        return preloadedData_;
    }
    bool usesPreloadedData( ) const
    {
        return usePreloadedData_;
    }
    basic_astrodynamics::OceanTideSynthesis getSynthesis( ) const
    {
        return synthesis_;
    }

protected:
    std::vector< std::string > blqFiles_;
    input_output::OceanLoadingData preloadedData_;
    bool usePreloadedData_ = false;
    basic_astrodynamics::OceanTideSynthesis synthesis_;
};

inline std::shared_ptr< BasicSolidBodyDeformationSettings > basicTidalBodyShapeDeformation(
        const std::vector< std::string > deformingBodies,
        const std::map< int, std::pair< double, double > > displacementLoveNumbers,
        const double bodyReferenceRadius = TUDAT_NAN )
{
    return std::make_shared< BasicSolidBodyDeformationSettings >( deformingBodies, displacementLoveNumbers, bodyReferenceRadius );
}

inline std::shared_ptr< BasicSolidBodyDeformationSettings > degreeTwoBasicTidalBodyShapeDeformation(
        const std::vector< std::string > deformingBodies,
        const double loveNumber,
        const double shidaNumber,
        const double bodyReferenceRadius = TUDAT_NAN )
{
    std::map< int, std::pair< double, double > > displacementLoveNumbers;
    displacementLoveNumbers[ 2 ] = std::make_pair( loveNumber, shidaNumber );
    return std::make_shared< BasicSolidBodyDeformationSettings >( deformingBodies, displacementLoveNumbers, bodyReferenceRadius );
}

inline std::shared_ptr< BodyDeformationSettings > iers2010TidalBodyShapeDeformation( )
{
    return std::make_shared< BodyDeformationSettings >( gravitation::iers_2010 );
}

inline std::shared_ptr< BodyDeformationSettings > oceanTideBodyShapeDeformation(
        const std::vector< std::string > blqFiles,
        const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp )
{
    return std::make_shared< OceanTideBodyDeformationSettings >( blqFiles, synthesis );
}

inline std::shared_ptr< BodyDeformationSettings > oceanTideBodyShapeDeformationFromData(
        const input_output::OceanLoadingData& data,
        const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp )
{
    return std::make_shared< OceanTideBodyDeformationSettings >( data, synthesis );
}

inline std::shared_ptr< BodyDeformationSettings > poleTideBodyShapeDeformation( )
{
    return std::make_shared< BodyDeformationSettings >( gravitation::pole_tide );
}

std::shared_ptr< basic_astrodynamics::BodyDeformationModel > createBodyDeformationModel(
        const std::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const SystemOfBodies& bodyMap );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_CREATEBODYDEFORMATIONMODEL_H
