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

#include <memory>

#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/iers2010SolidTidalBodyDeformation.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"


namespace tudat
{

namespace simulation_setup
{


std::shared_ptr< basic_astrodynamics::Iers2010EarthDeformation > createDefaultEarthIers2010DeformationModel(
        const std::shared_ptr< ephemerides::Ephemeris > earthEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > lunarEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > solarEphemeris,
        const std::shared_ptr< ephemerides::RotationalEphemeris > earthRotation,
        const std::function< double( ) > gravitionalParametersOfEarth,
        const std::function< double( ) > gravitionalParametersOfMoon,
        const std::function< double( ) > gravitionalParametersOfSun,
        const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction =
        [](const double time){return sofa_interface::calculateDoodsonFundamentalArguments( time ); } );


class BodyDeformationSettings
{
public:
    BodyDeformationSettings(
            const gravitation::BodyDeformationTypes bodyDeformationType ):bodyDeformationType_( bodyDeformationType ){ }

    virtual ~BodyDeformationSettings( ){ }

    gravitation::BodyDeformationTypes getBodyDeformationType( ){ return bodyDeformationType_; }

protected:
    gravitation::BodyDeformationTypes bodyDeformationType_;
};

class BasicSolidBodyDeformationSettings: public BodyDeformationSettings
{
public:
    BasicSolidBodyDeformationSettings( const std::vector< std::string > deformingBodies,
                                       const std::map< int, std::pair< double, double > > displacementLoveNumbers,
                                       const double bodyReferenceRadius = TUDAT_NAN ):
        BodyDeformationSettings( gravitation::basic_solid_body ),
        deformingBodies_( deformingBodies ),
        displacementLoveNumbers_( displacementLoveNumbers ),
        bodyReferenceRadius_( bodyReferenceRadius )
    { }

    std::vector< std::string > getDeformingBodies( ){ return deformingBodies_;}
    std::map< int, std::pair< double, double > > getDisplacementLoveNumbers( ){ return displacementLoveNumbers_; }
    double getBodyReferenceRadius( ){ return bodyReferenceRadius_; }

protected:
    std::vector< std::string > deformingBodies_;
    std::map< int, std::pair< double, double > > displacementLoveNumbers_;
    double bodyReferenceRadius_;
};

class OceanTideBodyDeformationSettings: public BodyDeformationSettings
{
public:
    OceanTideBodyDeformationSettings( const std::vector<std::string> blqFiles ):
        BodyDeformationSettings( gravitation::ocean_tide ),
        blqFiles_( blqFiles ){ }

    std::vector<std::string> getBlqFiles( ){ return blqFiles_; }

protected:

    std::vector<std::string> blqFiles_;
};

inline std::shared_ptr< BasicSolidBodyDeformationSettings > basicTidalBodyShapeDeformation (
        const std::vector< std::string > deformingBodies,
        const std::map< int, std::pair< double, double > > displacementLoveNumbers,
        const double bodyReferenceRadius = TUDAT_NAN )
{
    return std::make_shared< BasicSolidBodyDeformationSettings >(
                deformingBodies, displacementLoveNumbers, bodyReferenceRadius );
}

inline std::shared_ptr< BasicSolidBodyDeformationSettings > degreeTwoBasicTidalBodyShapeDeformation (
        const std::vector< std::string > deformingBodies,
        const double loveNumber,
        const double shidaNumber,
        const double bodyReferenceRadius = TUDAT_NAN )
{
    std::map< int, std::pair< double, double > > displacementLoveNumbers;
    displacementLoveNumbers[ 2 ] = std::make_pair( loveNumber, shidaNumber );
    return std::make_shared< BasicSolidBodyDeformationSettings >(
                deformingBodies, displacementLoveNumbers, bodyReferenceRadius );
}

inline std::shared_ptr< BodyDeformationSettings > iers2010TidalBodyShapeDeformation( )
{
    return std::make_shared< BodyDeformationSettings >( gravitation::iers_2010 );
}

inline std::shared_ptr< BodyDeformationSettings > oceanTideBodyShapeDeformation( const std::vector<std::string> blqFiles )
{
    return std::make_shared< OceanTideBodyDeformationSettings >( blqFiles );
}

inline std::shared_ptr< BodyDeformationSettings > poleTideBodyShapeDeformation( )
{
    return std::make_shared< BodyDeformationSettings >( gravitation::pole_tide );
}

std::shared_ptr< basic_astrodynamics::BodyDeformationModel > createBodyDeformationModel(
        const std::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const SystemOfBodies& bodyMap );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEBODYDEFORMATIONMODEL_H
