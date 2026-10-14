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
#include <variant>
#include <vector>

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/ground_stations/stationLoadingDeformation.h"
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

class StationLoadingBodyDeformationSettings : public BodyDeformationSettings
{
public:
    StationLoadingBodyDeformationSettings( const input_output::StationLoadingData& data,
                                           const basic_astrodynamics::StationLoadingEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::station_loading ), data_( data ), settings_( settings )
    {
        input_output::validateStationLoadingData( data_ );
    }
    const input_output::StationLoadingData& getData( ) const
    {
        return data_;
    }
    const basic_astrodynamics::StationLoadingEvaluationSettings& getEvaluationSettings( ) const
    {
        return settings_;
    }

private:
    input_output::StationLoadingData data_;
    basic_astrodynamics::StationLoadingEvaluationSettings settings_;
};

inline std::shared_ptr< StationLoadingBodyDeformationSettings > stationLoadingBodyShapeDeformation(
        const input_output::StationLoadingData& data,
        const basic_astrodynamics::StationLoadingEvaluationSettings& settings )
{
    return std::make_shared< StationLoadingBodyDeformationSettings >( data, settings );
}

//! The two physical atmospheric coefficient representations; file input is loaded at settings creation.
enum class AtmosphericTideInputFormat { grid, stationCoefficients };

//! Atmospheric loading settings containing grid or station-specific provider coefficients.
class AtmosphericTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    //! Copy grid data; geometry and origin policies are resolved when the model is created.
    AtmosphericTideBodyDeformationSettings( const input_output::AtmosphericTideGrid& grid,
                                            const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), data_( grid ), settings_( settings )
    {}

    //! Copy station-specific coefficients and their source metadata.
    AtmosphericTideBodyDeformationSettings( const input_output::AtmosphericTideStationData& stations,
                                            const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), data_( stations ), settings_( settings )
    {}

    //! Return the active grid representation; callers must select the matching input format.
    const input_output::AtmosphericTideGrid& getGrid( ) const
    {
        return std::get< input_output::AtmosphericTideGrid >( data_ );
    }

    //! Return the active station representation; callers must select the matching input format.
    const input_output::AtmosphericTideStationData& getStations( ) const
    {
        return std::get< input_output::AtmosphericTideStationData >( data_ );
    }

    //! Return the requested scientific geometry, origin and endpoint policies.
    const basic_astrodynamics::AtmosphericTideEvaluationSettings& getEvaluationSettings( ) const
    {
        return settings_;
    }

    //! Return whether these settings contain grid or station-specific coefficients.
    AtmosphericTideInputFormat getInputFormat( ) const
    {
        return std::holds_alternative< input_output::AtmosphericTideGrid >( data_ ) ? AtmosphericTideInputFormat::grid
                                                                                    : AtmosphericTideInputFormat::stationCoefficients;
    }

private:
    //! One provider data object, including its own metadata and file provenance.
    std::variant< input_output::AtmosphericTideGrid, input_output::AtmosphericTideStationData > data_;

    //! Policies resolved against the active provider data during environment creation.
    basic_astrodynamics::AtmosphericTideEvaluationSettings settings_;
};

//! Construct settings from a copied atmospheric grid.
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGrid(
        const input_output::AtmosphericTideGrid& grid,
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( grid, settings );
}

//! Construct settings from copied atmospheric station coefficients.
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStations(
        const input_output::AtmosphericTideStationData& stations,
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( stations, settings );
}

//! Read and validate the grid now; subsequent model creation performs no atmospheric file access.
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGridFile(
        const std::string& path,
        const input_output::AtmosphericTideMetadata& metadata = {},
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return atmosphericTideBodyShapeDeformationFromGrid( input_output::readAtmosphericTideGrid( path, metadata ), settings );
}

//! Read and validate station coefficients now; subsequent model creation uses the same preloaded settings.
inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStationFile(
        const std::string& path,
        const input_output::AtmosphericTideMetadata& metadata = {},
        const basic_astrodynamics::AtmosphericTideEvaluationSettings& settings = {} )
{
    return atmosphericTideBodyShapeDeformationFromStations( input_output::readAtmosphericTideStationData( path, metadata ), settings );
}

class OceanPoleTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    //! Copy the provider grid, preserving its declared reference origin.
    OceanPoleTideBodyDeformationSettings( const input_output::OceanPoleTideGrid& grid,
                                          const basic_astrodynamics::OceanPoleTideLatitudePolicy policy ):
        BodyDeformationSettings( gravitation::ocean_pole_tide ), grid_( grid ), policy_( policy )
    {}

    //! Return the single preloaded provider grid.
    const input_output::OceanPoleTideGrid& getGrid( ) const
    {
        return grid_;
    }

    //! Return the requested treatment of latitudes outside provider coverage.
    basic_astrodynamics::OceanPoleTideLatitudePolicy getLatitudePolicy( ) const
    {
        return policy_;
    }

private:
    //! Provider coefficients, coverage and provenance.
    input_output::OceanPoleTideGrid grid_;

    //! Strict coverage or explicitly requested latitude clamping.
    basic_astrodynamics::OceanPoleTideLatitudePolicy policy_;
};

//! Read the grid at settings creation, retaining the existing filename convenience interface.
inline std::shared_ptr< OceanPoleTideBodyDeformationSettings > oceanPoleTideBodyShapeDeformation(
        const std::string& gridFile,
        const basic_astrodynamics::OceanPoleTideLatitudePolicy policy =
                basic_astrodynamics::OceanPoleTideLatitudePolicy::requireGridCoverage )
{
    return std::make_shared< OceanPoleTideBodyDeformationSettings >( input_output::readOceanPoleTideGrid( gridFile ), policy );
}

//! Construct the same ocean pole-tide settings from preloaded grid data.
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
    //! Read and combine BLQ inputs now, retaining provenance and rejecting duplicate station identifiers.
    OceanTideBodyDeformationSettings(
            const std::vector< std::string > blqFiles,
            const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp );

    //! Copy BLQ amplitudes in metres and phase lags in degrees; empty data permits later model append.
    OceanTideBodyDeformationSettings(
            const input_output::OceanLoadingData& data,
            const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp ):
        BodyDeformationSettings( gravitation::ocean_tide ), data_( data ), synthesis_( synthesis )
    {
        input_output::validateOceanLoadingData( data_, true );
    }

    //! Return BLQ station coefficients and their provider provenance.
    const input_output::OceanLoadingData& getOceanLoadingData( ) const
    {
        return data_;
    }

    //! Return the selected conventional or explicit legacy synthesis method.
    basic_astrodynamics::OceanTideSynthesis getSynthesis( ) const
    {
        return synthesis_;
    }

    //! Return the BLQ files from which the data were read (empty for preloaded data).
    std::vector< std::string > getBlqFiles( )
    {
        return blqFiles_;
    }

protected:
    //! BLQ files from which the data were read, if any.
    std::vector< std::string > blqFiles_;

    //! BLQ coefficients, including provider headers and per-station source provenance.
    input_output::OceanLoadingData data_;

    //! Synthesis convention applied when the environment constructs the model.
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

//! Read one or more BLQ files at settings creation; later model creation uses only normalized station data.
inline std::shared_ptr< BodyDeformationSettings > oceanTideBodyShapeDeformation(
        const std::vector< std::string > blqFiles,
        const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp )
{
    return std::make_shared< OceanTideBodyDeformationSettings >( blqFiles, synthesis );
}

//! Construct ocean-loading settings from preloaded BLQ coefficients and source metadata.
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
