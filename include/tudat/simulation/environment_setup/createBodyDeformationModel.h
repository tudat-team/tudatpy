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

//! Settings for a body deformation model computed from a station loading series (for instance atmospheric loading).
class StationLoadingBodyDeformationSettings : public BodyDeformationSettings
{
public:
    StationLoadingBodyDeformationSettings( const input_output::StationLoadingData& data,
                                           const std::map< std::string, std::string >& stationBindings,
                                           const double maximumGapSeconds ):
        BodyDeformationSettings( gravitation::station_loading ), data_( data ), stationBindings_( stationBindings ),
        maximumGapSeconds_( maximumGapSeconds )
    {}

    //! Station loading series, the site identifier per ground station name, and the maximum interpolation gap (see
    //! basic_astrodynamics::StationLoadingDeformation).
    input_output::StationLoadingData data_;
    std::map< std::string, std::string > stationBindings_;
    double maximumGapSeconds_;
};

inline std::shared_ptr< StationLoadingBodyDeformationSettings > stationLoadingBodyShapeDeformation(
        const input_output::StationLoadingData& data,
        const std::map< std::string, std::string >& stationBindings = {},
        const double maximumGapSeconds = 86400.0 )
{
    return std::make_shared< StationLoadingBodyDeformationSettings >( data, stationBindings, maximumGapSeconds );
}

//! Settings for S1/S2 atmospheric tide loading from a grid or from station coefficients (see
//! basic_astrodynamics::AtmosphericTideDeformation).
class AtmosphericTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    AtmosphericTideBodyDeformationSettings( const input_output::AtmosphericTideGrid& grid, const bool convertToCentreOfMass ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), grid_( grid ), convertToCentreOfMass_( convertToCentreOfMass )
    {}

    AtmosphericTideBodyDeformationSettings( const std::map< std::string, input_output::AtmosphericTideCoefficients >& stationCoefficients,
                                            const bool convertToCentreOfMass ):
        BodyDeformationSettings( gravitation::atmospheric_tide ), stationCoefficients_( stationCoefficients ),
        convertToCentreOfMass_( convertToCentreOfMass )
    {}

    //! Grid of coefficients (empty if station coefficients are used), coefficients per station name, and whether the
    //! translation to the centre of mass of the whole Earth is added.
    input_output::AtmosphericTideGrid grid_;
    std::map< std::string, input_output::AtmosphericTideCoefficients > stationCoefficients_;
    bool convertToCentreOfMass_;
};

inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGrid(
        const input_output::AtmosphericTideGrid& grid,
        const bool convertToCentreOfMass = false )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( grid, convertToCentreOfMass );
}

inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStations(
        const std::map< std::string, input_output::AtmosphericTideCoefficients >& stationCoefficients,
        const bool convertToCentreOfMass = false )
{
    return std::make_shared< AtmosphericTideBodyDeformationSettings >( stationCoefficients, convertToCentreOfMass );
}

inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromGridFile(
        const std::string& path,
        const bool convertToCentreOfMass = false )
{
    return atmosphericTideBodyShapeDeformationFromGrid( input_output::readAtmosphericTideGrid( path ), convertToCentreOfMass );
}

inline std::shared_ptr< AtmosphericTideBodyDeformationSettings > atmosphericTideBodyShapeDeformationFromStationFile(
        const std::string& path,
        const bool convertToCentreOfMass = false )
{
    return atmosphericTideBodyShapeDeformationFromStations( input_output::readAtmosphericTideStationData( path ), convertToCentreOfMass );
}

//! Settings for the ocean pole tide loading deformation model (IERS 2010 Conventions, Section 7.1.5).
class OceanPoleTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    //! Copy the provider grid, preserving its declared reference origin.
    OceanPoleTideBodyDeformationSettings( const input_output::OceanPoleTideGrid& grid ):
        BodyDeformationSettings( gravitation::ocean_pole_tide ), grid_( grid )
    {}

    //! Return the single preloaded provider grid.
    const input_output::OceanPoleTideGrid& getGrid( ) const
    {
        return grid_;
    }

private:
    //! Provider coefficients, coverage and provenance.
    input_output::OceanPoleTideGrid grid_;
};

//! Read the grid at settings creation, retaining the existing filename convenience interface.
inline std::shared_ptr< OceanPoleTideBodyDeformationSettings > oceanPoleTideBodyShapeDeformation( const std::string& gridFile )
{
    return std::make_shared< OceanPoleTideBodyDeformationSettings >( input_output::readOceanPoleTideGrid( gridFile ) );
}

//! Construct the same ocean pole-tide settings from preloaded grid data.
inline std::shared_ptr< OceanPoleTideBodyDeformationSettings > oceanPoleTideBodyShapeDeformationFromData(
        const input_output::OceanPoleTideGrid& grid )
{
    return std::make_shared< OceanPoleTideBodyDeformationSettings >( grid );
}

//! Settings for the ocean tide loading deformation model, from ocean loading coefficients per station (BLQ format).
class OceanTideBodyDeformationSettings : public BodyDeformationSettings
{
public:
    OceanTideBodyDeformationSettings(
            const std::vector< std::string > blqFiles,
            const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp ):
        BodyDeformationSettings( gravitation::ocean_tide ), blqFiles_( blqFiles ), synthesis_( synthesis )
    {}

    //! Return the selected conventional or explicit legacy synthesis method.
    basic_astrodynamics::OceanTideSynthesis getSynthesis( ) const
    {
        return synthesis_;
    }

    //! Return the BLQ files with the ocean loading coefficients.
    std::vector< std::string > getBlqFiles( )
    {
        return blqFiles_;
    }

protected:
    //! BLQ files with the ocean loading coefficients.
    std::vector< std::string > blqFiles_;

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

inline std::shared_ptr< BodyDeformationSettings > oceanTideBodyShapeDeformation(
        const std::vector< std::string > blqFiles,
        const basic_astrodynamics::OceanTideSynthesis synthesis = basic_astrodynamics::OceanTideSynthesis::iersHardisp )
{
    return std::make_shared< OceanTideBodyDeformationSettings >( blqFiles, synthesis );
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
