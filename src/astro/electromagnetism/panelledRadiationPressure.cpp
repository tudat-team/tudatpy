//RP-OLD
#include <iostream>

#include "tudat/astro/electromagnetism/panelledRadiationPressure.h"

namespace tudat
{

namespace electromagnetism
{

//! Function to calculate radiation pressure force on a single partially reflecting panel
Eigen::Vector3d computeSinglePanelNormalizedRadiationPressureForce(
        const Eigen::Vector3d& normalizedVectorToSource, const Eigen::Vector3d& panelSurfaceNormal,
        const double panelArea, const double panelEmissivitty, const double panelDiffuseReflectionCoefficient )
{
    // Calculate cosine of the angle between panel surface normal and vector from accelerated to radiating body
    double cosineOfPanelInclination = normalizedVectorToSource.dot( panelSurfaceNormal );

    // Initialize force to zero
    Eigen::Vector3d panelRadiationPressureForce = Eigen::Vector3d::Zero( );

    // If cosineOfPanelInclination is larger than zero (i.e inclination is smaller than 90 degrees), calculated acceleration force (zero otherwise)
    if( cosineOfPanelInclination > 0.0 )
    {
        // Evaluate Eq. (3.72) of Montenbruck & Gill (2000)
        panelRadiationPressureForce = -cosineOfPanelInclination * panelArea * (
                    ( 1.0 - panelEmissivitty ) * normalizedVectorToSource + 2.0 * (
                        panelEmissivitty * cosineOfPanelInclination + panelDiffuseReflectionCoefficient / 3.0 ) * panelSurfaceNormal );
    }

    return panelRadiationPressureForce;
}

//! Constructor for setting up the acceleration model, with input RadiationPressureInterface (and massFunction)
PanelledRadiationPressureAcceleration::PanelledRadiationPressureAcceleration(
        const std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface,
        const std::function< double( ) > massFunction ):
    massFunction_( massFunction )
{
    sourcePositionFunction_ =
            std::bind( &RadiationPressureInterface::getCurrentSolarVector, radiationPressureInterface );
    acceleratedBodyPositionFunction_ = [ ]( ){ return Eigen::Vector3d::Zero( ); };

    radiationPressureFunction_ = std::bind( &RadiationPressureInterface::getCurrentRadiationPressure, radiationPressureInterface );

    numberOfPanels_ = radiationPressureInterface->getNumberOfPanels( );

    for( int i = 0; i < radiationPressureInterface->getNumberOfPanels( ); i++ )
    {
        panelEmissivittyFunctions_.push_back(
                    std::bind( &PanelledRadiationPressureInterface::getEmissivity, radiationPressureInterface, i ) );
        panelDiffuseReflectionCoefficientFunctions_.push_back(
                    std::bind( &PanelledRadiationPressureInterface::getDiffuseReflectionCoefficient, radiationPressureInterface, i ) );
        panelSurfaceNormalFunctions_.push_back(
                    std::bind( &PanelledRadiationPressureInterface::getCurrentSurfaceNormal, radiationPressureInterface, i ) );
        panelAreaFunctions_.push_back(
                    std::bind( &PanelledRadiationPressureInterface::getArea, radiationPressureInterface, i ) );
    }

    currentPanelAccelerations_.resize( radiationPressureInterface->getNumberOfPanels( ) );
}

}

}
