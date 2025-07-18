/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
#ifndef TUDAT_GASSURFACEINTERACTIONMODEL_H
#define TUDAT_GASSURFACEINTERACTIONMODEL_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"

namespace tudat
{
namespace aerodynamics
{

enum GasSurfaceInteractionModelType
{
    constantCoefficients,
    newton,
    storch,
    sentman,
    cook
};

class GasSurfaceInteractionModel
{
public:
GasSurfaceInteractionModel( const GasSurfaceInteractionModelType modelType,
                            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                            const double referenceArea,
                            const int maximumNumberOfPixels, 
                            const bool onlyDrag ):
                            modelType_( modelType), allPanels_( allPanels ),
                            referenceArea_( referenceArea ),
                            maximumNumberOfPixels_( maximumNumberOfPixels ),
                            onlyDrag_( onlyDrag )                     
{ 
    totalNumberOfPanels_ = allPanels_.size( );
    surfacePanelCosines_.resize( totalNumberOfPanels_ );
    if ( maximumNumberOfPixels_ > 2 )
    {
        aerodynamicSelfShadowing_ = std::make_shared< tudat::system_models::SelfShadowing >( allPanels_, maximumNumberOfPixels_ );
    }
    else
    {
        aerodynamicSelfShadowing_ = nullptr;
    }
    unityIlluminationFraction_ = std::vector< double >( totalNumberOfPanels_, 1.0);
}

virtual Eigen::Vector3d computeAerodynamicCoefficients( ) = 0;

virtual ~GasSurfaceInteractionModel( ) = default;

void updateMembers( );

std::vector< double >& getSurfacePanelCosines( )
{
    return surfacePanelCosines_;
}

std::vector< double >& getIlluminatedPanelFractions( )
{
    return illuminatedPanelFractions_;
}

double getReferenceArea( ) const
{
    return referenceArea_;
}

void setIncomingDirection( const Eigen::Vector3d incomingDirection )
{
    incomingDirection_ = incomingDirection;
}

void setFreeStreamTemperature( const double freeStreamTemperature )
{
    freeStreamTemperature_ = freeStreamTemperature;
}

void setAirSpeed( const double airSpeed )
{
    airSpeed_ = airSpeed;
}

void setSpecifiGasConstant( const double specificGasConstant )
{
    specificGasConstant_ = specificGasConstant;
}

GasSurfaceInteractionModelType getGasSurfaceInteractionModelType( ) const
{
    return modelType_;
}

void setConstantAerodynamicCoefficients( const Eigen::Vector3d constantAerodynamicCoefficients )
{
    constantAerodynamicCoefficients_ = constantAerodynamicCoefficients;
}

protected:
GasSurfaceInteractionModelType modelType_;

const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels_;

std::shared_ptr< tudat::system_models::SelfShadowing > aerodynamicSelfShadowing_;

double referenceArea_;

int maximumNumberOfPixels_;

bool onlyDrag_;

Eigen::Vector3d incomingDirection_;

std::vector< double > unityIlluminationFraction_;

std::vector< double > illuminatedPanelFractions_;

std::vector< double > surfacePanelCosines_;

int totalNumberOfPanels_;

double speedRatio_;

double airSpeed_;

double freeStreamTemperature_;

double incidentTemperature_;

double specificGasConstant_;

// for constant coefficients but variable cross-section
Eigen::Vector3d constantAerodynamicCoefficients_;

};

class ConstantInteractionModel : public GasSurfaceInteractionModel
{
public:
ConstantInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                  const double referenceArea,
                                  const int maximumNumberOfPixels, 
                                  const bool onlyDrag ):
                                  GasSurfaceInteractionModel( constantCoefficients, allPanels, referenceArea, maximumNumberOfPixels, onlyDrag )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class NewtonGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
NewtonGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                  const double referenceArea,
                                  const int maximumNumberOfPixels, 
                                  const bool onlyDrag ):
                                  GasSurfaceInteractionModel( newton, allPanels, referenceArea, maximumNumberOfPixels, onlyDrag )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class StorchGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
StorchGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                  const double referenceArea,
                                  const int maximumNumberOfPixels,
                                  const bool onlyDrag ):
                                  GasSurfaceInteractionModel( storch, allPanels, referenceArea, maximumNumberOfPixels, onlyDrag )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class SentmanGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
SentmanGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                   const double referenceArea,
                                   const int maximumNumberOfPixels,
                                   const bool onlyDrag ):
                                   GasSurfaceInteractionModel( sentman, allPanels, referenceArea, maximumNumberOfPixels, onlyDrag )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

class CookGasSurfaceInteractionModel : public GasSurfaceInteractionModel
{
public:
CookGasSurfaceInteractionModel( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                const double referenceArea,
                                const int maximumNumberOfPixels,
                                const bool onlyDrag ):
                                GasSurfaceInteractionModel( cook, allPanels, referenceArea, maximumNumberOfPixels, onlyDrag )
{ }

Eigen::Vector3d computeAerodynamicCoefficients( );

};

inline std::vector< AerodynamicCoefficientsIndependentVariables > createIndependentVariablesNamesForGasSurfaceInteractionModel( 
    const GasSurfaceInteractionModelType gasSurfaceInteractionModelType )
{
    switch( gasSurfaceInteractionModelType )
    {
        case constantCoefficients: {
            return {
                angle_of_attack_dependent,
                angle_of_sideslip_dependent
            };
        }
        case newton: {
            return {
                angle_of_attack_dependent,
                angle_of_sideslip_dependent
            };
        }
        case storch: {
            return {
                angle_of_attack_dependent,
                angle_of_sideslip_dependent
            };
        }
        case sentman: {
            return {
                angle_of_attack_dependent,
                angle_of_sideslip_dependent,
                temperature_dependent,
                velocity_dependent
            };
        }
        case cook: {
            return {
                angle_of_attack_dependent,
                angle_of_sideslip_dependent,
                temperature_dependent
            };
        }
        default:
            throw std::runtime_error( "Error, unknown gas surface interaction model " );
    }
}

inline std::shared_ptr< GasSurfaceInteractionModel > createGasSurfaceInteractionModel( const GasSurfaceInteractionModelType gasSurfaceInteractionModelType,
                                                                                const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
                                                                                const double referenceArea,
                                                                                const int maximumNumberOfPixels,
                                                                                const bool onlyDrag )
{
    switch( gasSurfaceInteractionModelType )
    {
        case constantCoefficients: {
            return std::make_shared< ConstantInteractionModel >( allPanels, referenceArea, maximumNumberOfPixels, onlyDrag );
        }
        case newton: {
            return std::make_shared< NewtonGasSurfaceInteractionModel >( allPanels, referenceArea, maximumNumberOfPixels, onlyDrag );
        }
        case storch: {
            // checking material properties
            for ( auto panel: allPanels )
            {
                if ( panel->getNormalAccomodationCoefficient( ) == -1 )
                {
                    throw std::runtime_error( "Error, normal accomodation coefficient for panel type " + panel->getPanelTypeId( ) +
                            " not defined but is required for Storch model" );
                }
                if ( panel->getTangentialAccomodationCoefficient( ) == -1 )
                {
                    throw std::runtime_error( "Error, tangential accomodation coefficient for panel type " + panel->getPanelTypeId( ) +
                            " not defined but is required for Storch model" );
                }
                if ( panel->getNormalVelocityAtWallRatio( ) == -1 )
                {
                    throw std::runtime_error( "Error, normal velocity ratio for panel type " + panel->getPanelTypeId( ) +
                            " not defined but is required for Storch model" );
                }
            }
            return std::make_shared< StorchGasSurfaceInteractionModel >( allPanels, referenceArea, maximumNumberOfPixels, onlyDrag );
        }
        case sentman: {
            // checking material properties
            for ( auto panel: allPanels )
            {
                if ( panel->getEnergyAccomodationCoefficient( ) == -1 )
                {
                    throw std::runtime_error( "Error, energy accomodation coefficient for panel type " + panel->getPanelTypeId( ) +
                            " not defined but is required for Sentman model" );
                }
            }
            return std::make_shared< SentmanGasSurfaceInteractionModel >( allPanels, referenceArea, maximumNumberOfPixels, onlyDrag );
        }
        case cook: {
            // checking material properties
            for ( auto panel: allPanels )
            {
                if ( panel->getEnergyAccomodationCoefficient( ) == -1 )
                {
                    throw std::runtime_error( "Error, energy accomodation coefficient for panel type " + panel->getPanelTypeId( ) +
                            " not defined but is required for Cook model" );
                }
            }
            return std::make_shared< CookGasSurfaceInteractionModel >( allPanels, referenceArea, maximumNumberOfPixels, onlyDrag );
        }
        default:
            throw std::runtime_error( "Error, unknown gas surface interaction model " );
    }
}


}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_GASSURFACEINTERACTIONMODEL_H