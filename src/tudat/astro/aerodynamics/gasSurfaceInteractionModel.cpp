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
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <Eigen/Core>
#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/aerodynamics/gasSurfaceInteractionModel.h"

namespace tudat
{
namespace aerodynamics
{

void GasSurfaceInteractionModel::updateMembers( )
{
    if( maximumNumberOfPixels_ == 0 )
    {
        // SSH off
        illuminatedPanelFractions_ = unityIlluminationFraction_;
    }
    else
    {
        // SSH on
        aerodynamicSelfShadowing_->reset( );
        aerodynamicSelfShadowing_->updateIlluminatedPanelFractions( incomingDirection_ );
        illuminatedPanelFractions_ = aerodynamicSelfShadowing_->getIlluminatedPanelFractions( );
    }
}

Eigen::Vector3d ConstantInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );

    double cosineDelta;
    double actualCrossSectionArea = 0;
    Eigen::Vector3d panelNormal;

    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        actualCrossSectionArea += allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
    }
    referenceArea_ = actualCrossSectionArea;
    return constantAerodynamicCoefficients_;
}

Eigen::Vector3d NewtonGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );

    double cosineDelta;
    double Cp, panelArea;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // Cp
        Cp = 2 * cosineDelta * cosineDelta;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame += -Cp * panelNormal * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d StorchGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    double cosineDelta, sineDelta;
    double Cp, Ct, panelArea;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;
    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt( std::max( 0.0, 1 - cosineDelta * cosineDelta ) );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // Cp
        Cp = 2 * cosineDelta *
                ( allPanels_[ i ]->getNormalAccomodationCoefficient( ) * allPanels_[ i ]->getNormalVelocityAtWallRatio( ) +
                  ( 2 - allPanels_[ i ]->getNormalAccomodationCoefficient( ) ) * cosineDelta );
        // Ct
        Ct = 2 * allPanels_[ i ]->getTangentialAccomodationCoefficient( ) * sineDelta * cosineDelta;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame +=
                ( -Cp * panelNormal - Ct * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d StorchGasSurfaceInteractionModel::computeAerodynamicCoefficientsPartial( const PanelMaterialPropertyType propertyType,
                                                                                         const std::string& panelGroupId )
{
    // Replays the loop of computeAerodynamicCoefficients( ), accumulating the derivative of the per-panel Cp, Ct w.r.t. the
    // selected material property for the panels belonging to panelGroupId. The geometry (cosines, illuminated fractions,
    // reference area) is reused from the most recent forward evaluation.
    double cosineDelta, sineDelta;
    double dCp, dCt, panelArea;
    Eigen::Vector3d coefficientPartialBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;
    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        if( allPanels_[ i ]->getPanelTypeId( ) != panelGroupId )
        {
            continue;
        }
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        if( cosineDelta == 0 )
        {
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt( std::max( 0.0, 1 - cosineDelta * cosineDelta ) );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        dCp = 0.0;
        dCt = 0.0;
        switch( propertyType )
        {
            case normal_accommodation_property:
                // Cp = 2 cosD ( sigma_n V_w + (2 - sigma_n) cosD ) -> dCp/dsigma_n = 2 cosD ( V_w - cosD )
                dCp = 2 * cosineDelta * ( allPanels_[ i ]->getNormalVelocityAtWallRatio( ) - cosineDelta );
                break;
            case normal_velocity_ratio_property:
                // dCp/dV_w = 2 cosD sigma_n
                dCp = 2 * cosineDelta * allPanels_[ i ]->getNormalAccomodationCoefficient( );
                break;
            case tangential_accommodation_property:
                // Ct = 2 sigma_t sinD cosD -> dCt/dsigma_t = 2 sinD cosD
                dCt = 2 * sineDelta * cosineDelta;
                break;
            default:
                // Energy accommodation does not enter the Storch model: partial is exactly zero.
                break;
        }
        coefficientPartialBodyFrame +=
                ( -dCp * panelNormal - dCt * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    coefficientPartialBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        coefficientPartialBodyFrame = incomingDirection_ * ( coefficientPartialBodyFrame.dot( incomingDirection_ ) );
    }
    return coefficientPartialBodyFrame;
}

Eigen::Vector3d SentmanGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    speedRatio_ = airSpeed_ / std::sqrt( 2 * specificGasConstant_ * freeStreamTemperature_ );
    double cosineDelta, sineDelta;
    double Cp, Ct, panelArea;
    double erf, exp;
    double sqrtPi = std::sqrt( mathematical_constants::PI );
    incidentTemperature_ = 2.0 / 3.0 * speedRatio_ * speedRatio_ * freeStreamTemperature_;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        erf = std::erf( speedRatio_ * cosineDelta );
        exp = std::exp( -speedRatio_ * speedRatio_ * cosineDelta * cosineDelta );
        sineDelta = std::sqrt( std::max( 0.0, 1 - cosineDelta * cosineDelta ) );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // Cp
        Cp = ( cosineDelta * cosineDelta ) * ( 1 + erf ) + cosineDelta / ( speedRatio_ * sqrtPi ) * exp +
                0.5 *
                        std::sqrt( 2.0 / 3.0 *
                                   ( 1 +
                                     ( allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * allPanels_[ i ]->getPanelTemperature( ) ) /
                                             ( incidentTemperature_ - 1 ) ) ) *
                        ( sqrtPi * cosineDelta * ( 1 + erf ) + 1.0 / speedRatio_ * exp );
        // Ct
        Ct = sineDelta * cosineDelta * ( 1 + erf ) + sineDelta / ( speedRatio_ * sqrtPi ) * exp;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame +=
                ( -Cp * panelNormal - Ct * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d SentmanGasSurfaceInteractionModel::computeAerodynamicCoefficientsPartial( const PanelMaterialPropertyType propertyType,
                                                                                          const std::string& panelGroupId )
{
    // Only the energy accommodation coefficient enters the Sentman model; the partial w.r.t. the other properties is zero.
    if( propertyType != energy_accommodation_property )
    {
        return Eigen::Vector3d::Zero( );
    }

    speedRatio_ = airSpeed_ / std::sqrt( 2 * specificGasConstant_ * freeStreamTemperature_ );
    incidentTemperature_ = 2.0 / 3.0 * speedRatio_ * speedRatio_ * freeStreamTemperature_;
    double sqrtPi = std::sqrt( mathematical_constants::PI );
    double cosineDelta;
    double dCp, dCt, panelArea;
    double erf, exp;
    Eigen::Vector3d coefficientPartialBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;
    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        if( allPanels_[ i ]->getPanelTypeId( ) != panelGroupId )
        {
            continue;
        }
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        if( cosineDelta == 0 )
        {
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        erf = std::erf( speedRatio_ * cosineDelta );
        exp = std::exp( -speedRatio_ * speedRatio_ * cosineDelta * cosineDelta );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // Cp = ... + 0.5 sqrt(B) ( sqrtPi cosD (1+erf) + (1/S) exp ), with B = 2/3 ( 1 + alpha_e Tp / (Ti - 1) ).
        // d sqrt(B)/d alpha_e = ( dB/d alpha_e ) / ( 2 sqrt(B) ), dB/d alpha_e = 2/3 Tp / (Ti - 1).
        double panelTemperature = allPanels_[ i ]->getPanelTemperature( );
        double bTerm = 2.0 / 3.0 *
                ( 1 + ( allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * panelTemperature ) / ( incidentTemperature_ - 1 ) );
        double bTermDerivative = 2.0 / 3.0 * panelTemperature / ( incidentTemperature_ - 1 );
        dCp = 0.5 * bTermDerivative / ( 2.0 * std::sqrt( bTerm ) ) * ( sqrtPi * cosineDelta * ( 1 + erf ) + 1.0 / speedRatio_ * exp );
        // Ct is independent of the energy accommodation coefficient.
        dCt = 0.0;
        coefficientPartialBodyFrame +=
                ( -dCp * panelNormal - dCt * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    coefficientPartialBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        coefficientPartialBodyFrame = incomingDirection_ * ( coefficientPartialBodyFrame.dot( incomingDirection_ ) );
    }
    return coefficientPartialBodyFrame;
}

Eigen::Vector3d CookGasSurfaceInteractionModel::computeAerodynamicCoefficients( )
{
    updateMembers( );
    double cosineDelta, sineDelta;
    double Cd, Cl, Cp, Ct, panelArea;
    double sqrt;
    Eigen::Vector3d currentForceCoefficientsBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;

    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        surfacePanelCosines_[ i ] = cosineDelta;
        if( cosineDelta == 0 )
        {
            illuminatedPanelFractions_[ i ] = 0.0;
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt( std::max( 0.0, 1 - cosineDelta * cosineDelta ) );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        sqrt = std::sqrt( 1 +
                          ( allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * allPanels_[ i ]->getPanelTemperature( ) ) /
                                  ( freeStreamTemperature_ - 1 ) );
        // Cd
        Cd = 2 * cosineDelta * ( 1 + 2.0 / 3.0 * cosineDelta * sqrt );
        // Cl
        Cl = 4.0 / 3.0 * sineDelta * cosineDelta * sqrt;
        // convert cd, cd to cp, ct
        Cp = cosineDelta * Cd + sineDelta * Cl;
        Ct = sineDelta * Cd - cosineDelta * Cl;
        // add panel contribution to the force coefficient vector (pre-multiplied by the surface
        // of the panel, which will be removed later)
        currentForceCoefficientsBodyFrame +=
                ( -Cp * panelNormal - Ct * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    // divide by total reference area to obtain true aerodynamic coefficients
    currentForceCoefficientsBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        currentForceCoefficientsBodyFrame = incomingDirection_ * ( currentForceCoefficientsBodyFrame.dot( incomingDirection_ ) );
    }
    return currentForceCoefficientsBodyFrame;
}

Eigen::Vector3d CookGasSurfaceInteractionModel::computeAerodynamicCoefficientsPartial( const PanelMaterialPropertyType propertyType,
                                                                                       const std::string& panelGroupId )
{
    // Only the energy accommodation coefficient enters the Cook model; the partial w.r.t. the other properties is zero.
    if( propertyType != energy_accommodation_property )
    {
        return Eigen::Vector3d::Zero( );
    }

    double cosineDelta, sineDelta;
    double dCd, dCl, dCp, dCt, panelArea;
    Eigen::Vector3d coefficientPartialBodyFrame = Eigen::Vector3d::Zero( );
    Eigen::Vector3d panelNormal;
    for( int i = 0; i < totalNumberOfPanels_; i++ )
    {
        if( allPanels_[ i ]->getPanelTypeId( ) != panelGroupId )
        {
            continue;
        }
        panelNormal = allPanels_[ i ]->getBodyFixedSurfaceNormal( )( );
        cosineDelta = panelNormal.dot( -incomingDirection_ );
        cosineDelta = cosineDelta > 0.0 ? cosineDelta : 0.0;
        if( cosineDelta == 0 )
        {
            continue;
        }
        if( illuminatedPanelFractions_[ i ] == 0.0 )
        {
            continue;
        }
        sineDelta = std::sqrt( std::max( 0.0, 1 - cosineDelta * cosineDelta ) );
        panelArea = allPanels_[ i ]->getPanelArea( ) * illuminatedPanelFractions_[ i ];
        // sqrtTerm = sqrt( 1 + alpha_e Tp / (T_inf - 1) ); d sqrtTerm/d alpha_e = ( Tp / (T_inf - 1) ) / ( 2 sqrtTerm ).
        double panelTemperature = allPanels_[ i ]->getPanelTemperature( );
        double sqrtTerm = std::sqrt(
                1 + ( allPanels_[ i ]->getEnergyAccomodationCoefficient( ) * panelTemperature ) / ( freeStreamTemperature_ - 1 ) );
        double sqrtTermDerivative = ( panelTemperature / ( freeStreamTemperature_ - 1 ) ) / ( 2.0 * sqrtTerm );
        // Cd = 2 cosD ( 1 + 2/3 cosD sqrtTerm ), Cl = 4/3 sinD cosD sqrtTerm.
        dCd = 4.0 / 3.0 * cosineDelta * cosineDelta * sqrtTermDerivative;
        dCl = 4.0 / 3.0 * sineDelta * cosineDelta * sqrtTermDerivative;
        // Cp = cosD Cd + sinD Cl, Ct = sinD Cd - cosD Cl.
        dCp = cosineDelta * dCd + sineDelta * dCl;
        dCt = sineDelta * dCd - cosineDelta * dCl;
        coefficientPartialBodyFrame +=
                ( -dCp * panelNormal - dCt * ( incomingDirection_.cross( panelNormal ) ).cross( panelNormal ) ) * panelArea;
    }
    coefficientPartialBodyFrame /= referenceArea_;
    if( onlyDrag_ )
    {
        coefficientPartialBodyFrame = incomingDirection_ * ( coefficientPartialBodyFrame.dot( incomingDirection_ ) );
    }
    return coefficientPartialBodyFrame;
}

}  // namespace aerodynamics

}  // namespace tudat
