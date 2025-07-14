/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#ifndef TUDAT_RADIATIONPRESSURETARGETMODEL_H
#define TUDAT_RADIATIONPRESSURETARGETMODEL_H

#include <map>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/selfShadowing.h"

namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling a target that is accelerated by radiation pressure.
 *
 * Parameters to all functions are in local target-fixed coordinates. Proper translation and rotation is ensured by the
 * RadiationPressureAcceleration class.
 */
class RadiationPressureTargetModel
{
public:
    explicit RadiationPressureTargetModel( const std::map< std::string, std::vector< std::string > >& sourceToTargetOccultingBodies = { } ):
        sourceToTargetOccultingBodies_( sourceToTargetOccultingBodies ), computeTorques_( false ), centerOfMassFunction_( nullptr )
    { }

    virtual ~RadiationPressureTargetModel( ) = default;

    virtual void enableTorqueComputation( const std::function< Eigen::Vector3d( ) > centerOfMassFunction ) = 0;

    void updateMembers( double currentTime );

    /*!
     * Calculate radiation pressure force from incident radiation using geometrical/optical target properties.
     *
     * @param sourceIrradiance Incident irradiance magnitude [W/m²]
     * @param sourceToTargetDirectionLocalFrame Direction of incoming radiation
     * @return Radiation pressure force vector in local (i.e. target-fixed) coordinates [N]
     */
    virtual void updateRadiationPressureForcing( const double sourceIrradiance,
                                                 const Eigen::Vector3d& sourceToTargetDirection,
                                                 const bool resetForces,
                                                 const std::string sourceName = "" ) = 0;

    std::map< std::string, std::vector< std::string > > getSourceToTargetOccultingBodies( ) const
    {
        return sourceToTargetOccultingBodies_;
    }

    virtual bool forceFunctionRequiresLocalFrameInputs( ) = 0;

    Eigen::Vector3d getCurrentRadiationPressureForce( const std::string& sourceName = "" )
    {
        return currentRadiationPressureForce_.at( sourceName );
    }

    Eigen::Vector3d getCurrentRadiationPressureTorque( const std::string& sourceName = "" )
    {
        return currentRadiationPressureTorque_.at( sourceName );
    }

    Eigen::Vector3d updateAndGetRadiationPressureForce( const double sourceIrradiance,
                                                        const Eigen::Vector3d& sourceToTargetDirection,
                                                        const bool resetForces,
                                                        const std::string sourceName = "" )
    {
        updateRadiationPressureForcing( sourceIrradiance, sourceToTargetDirection, resetForces, sourceName );
        return currentRadiationPressureForce_.at( sourceName );
    }

    Eigen::Vector3d updateAndGetRadiationPressureTorque( const double sourceIrradiance,
                                                         const Eigen::Vector3d& sourceToTargetDirection,
                                                         const bool resetForces,
                                                         const std::string sourceName = "" )
    {
        updateRadiationPressureForcing( sourceIrradiance, sourceToTargetDirection, resetForces, sourceName );
        return currentRadiationPressureTorque_.at( sourceName );
    }

    virtual void resetDerivedComputations( const std::string sourceName ) { }

    void resetComputations( const std::string& sourceName )
    {
        currentRadiationPressureForce_[ sourceName ] = Eigen::Vector3d::Zero( );
        currentRadiationPressureTorque_[ sourceName ] = Eigen::Vector3d::Zero( );

        resetDerivedComputations( sourceName );
    }

    virtual void saveLocalComputations( const std::string sourceName, const bool saveCosines ) { }

protected:
    virtual void updateMembers_( const double currentTime ) {};

    double currentTime_{ TUDAT_NAN };
    // Only needed to transfer occultation settings from body setup to acceleration setup
    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies_;

    // Source-specific variables
    std::map< std::string, Eigen::Vector3d > currentRadiationPressureForce_;
    std::map< std::string, Eigen::Vector3d > currentRadiationPressureTorque_;
    bool computeTorques_;

    std::function< Eigen::Vector3d( ) > centerOfMassFunction_;
};

/*!
 * Class modeling a target as sphere ("cannonball"). The sphere has an isotropic radiation pressure coefficient and
 * the same cross-sectional area from any direction.
 */
class CannonballRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    /*!
     * Constructor.
     *
     * @param area Cross-sectional area (i.e. projected area of the sphere)
     * @param coefficient Radiation pressure coefficient (between 1 [pure absorption] and 2 [pure specular reflection])
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    CannonballRadiationPressureTargetModel(
            double area,
            double coefficient,
            const std::map< std::string, std::vector< std::string > >& sourceToTargetOccultingBodies = { } ):
        RadiationPressureTargetModel( sourceToTargetOccultingBodies ), area_( area ), coefficientFunction_( nullptr ),
        currentCoefficient_( coefficient )
    { }

    CannonballRadiationPressureTargetModel(
            double area,
            std::function< double( const double ) > coefficientFunction,
            const std::map< std::string, std::vector< std::string > >& sourceToTargetOccultingBodies = { } ):
        RadiationPressureTargetModel( sourceToTargetOccultingBodies ), area_( area ), coefficientFunction_( coefficientFunction ),
        currentCoefficient_( TUDAT_NAN )
    { }

    void enableTorqueComputation( const std::function< Eigen::Vector3d( ) > centerOfMassFunction ) override
    {
        if( centerOfMassFunction == nullptr )
        {
            throw std::runtime_error( "Error when enable radiation pressure torque in panelled model; center of mass function is null." );
        }

        if( computeTorques_ )
        {
            throw std::runtime_error( "Error when enable radiation pressure torque in panelled model; torque is already enabled." );
        }

        computeTorques_ = true;
        centerOfMassFunction_ = centerOfMassFunction;
    }

    void updateRadiationPressureForcing( double sourceIrradiance,
                                         const Eigen::Vector3d& sourceToTargetDirection,
                                         const bool resetForces,
                                         const std::string sourceName = "" ) override;

    double getArea( ) const
    {
        return area_;
    }

    double getCoefficient( ) const
    {
        return currentCoefficient_;
    }

    void resetCoefficient( const double coefficient )
    {
        currentCoefficient_ = coefficient;
    }

    std::function< double( const double ) > getCoefficientFunction( )
    {
        return coefficientFunction_;
    }

    void resetCoefficientFunction( const std::function< double( const double ) > coefficientFunction )
    {
        coefficientFunction_ = coefficientFunction;
    }

    bool forceFunctionRequiresLocalFrameInputs( ) override
    {
        return false;
    }

private:
    virtual void updateMembers_( const double currentTime ) override
    {
        if( coefficientFunction_ != nullptr )
        {
            currentCoefficient_ = coefficientFunction_( currentTime );
        }
    };

    double area_;

    std::function< double( const double ) > coefficientFunction_;

    double currentCoefficient_;
};

/*!
 * Class modeling a target as collection of panels, e.g., representing the box body and solar panels.
 */
class PaneledRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    /*!
     * Constructor.
     *
     * @param panels Panels comprising this paneled target
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    explicit PaneledRadiationPressureTargetModel(
            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& bodyFixedPanels,
            const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels,
            const std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >& segmentFixedPanels =
                    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >( ),
            const std::map< std::string, std::function< Eigen::Quaterniond( ) > >& segmentFixedToBodyFixedRotations =
                    std::map< std::string, std::function< Eigen::Quaterniond( ) > >( ),
            const std::map< std::string, std::vector< std::string > >& sourceToTargetOccultingBodies = { },
            const std::map< std::string, int > maximumNumberOfPixelsPerSource = { },
            bool panelGeometryDefined = false ):
        RadiationPressureTargetModel( sourceToTargetOccultingBodies ), bodyFixedPanels_( bodyFixedPanels ),
        segmentFixedPanels_( segmentFixedPanels ), segmentFixedToBodyFixedRotations_( segmentFixedToBodyFixedRotations ),
        maximumNumberOfPixelsPerSource_( maximumNumberOfPixelsPerSource ), allPanels_( allPanels ),
        panelGeometryDefined_( panelGeometryDefined )
    {
        totalNumberOfPanels_ = allPanels.size( );
        panelForces_.resize( totalNumberOfPanels_ );
        surfacePanelCosines_.resize( totalNumberOfPanels_ );
        surfaceNormals_.resize( totalNumberOfPanels_ );
        illuminatedPanelFractions_.resize( totalNumberOfPanels_ );

        // create panelTypeIdList for dependent variable save
        panelTypeIdList_ = std::vector< std::string >( totalNumberOfPanels_ );
        for( int i = 0; i < totalNumberOfPanels_; i++ )
        {
            panelTypeIdList_[ i ] = allPanels_[ i ]->getPanelTypeId( );
        }
        // create map of self-shadowing objects per source
        for( auto it: maximumNumberOfPixelsPerSource_ )
        {
            if( it.second < 0 || it.second == 1 )
            {
                throw std::runtime_error( "Error, invalid maximum number of pixels for source " + it.first + ", value should be > 2." );
            }
            if( it.second > 1 && !panelGeometryDefined_ )
            {
                throw std::runtime_error( "Error, maximum number of pixels given for source " + it.first +
                                          ", however, no panel geometry is defined." );
            }
            if( it.second == 0 )
            {
                continue;
            }
            selfShadowingPerSource_[ it.first ] =
                    std::make_shared< system_models::SelfShadowing >( system_models::SelfShadowing( allPanels_, it.second ) );
        }

        unityIlluminationFraction_ = std::vector< double >( totalNumberOfPanels_, 1.0 );
    }

    void enableTorqueComputation( const std::function< Eigen::Vector3d( ) > centerOfMassFunction ) override
    {
        if( centerOfMassFunction == nullptr )
        {
            throw std::runtime_error( "Error when enable radiation pressure torque in panelled model; center of mass function is null." );
        }

        if( computeTorques_ )
        {
            throw std::runtime_error( "Error when enable radiation pressure torque in panelled model; torque is already enabled." );
        }

        computeTorques_ = true;
        centerOfMassFunction_ = centerOfMassFunction;
        panelTorques_.resize( totalNumberOfPanels_ );
        panelCentroidMomentArms_.resize( totalNumberOfPanels_ );
    }

    void updateRadiationPressureForcing( double sourceIrradiance,
                                         const Eigen::Vector3d& sourceToTargetDirectionLocalFrame,
                                         const bool resetForces,
                                         const std::string sourceName = "" ) override;

    Eigen::Vector3d evaluateRadiationPressureForcePartialWrtDiffuseReflectivity( double sourceIrradiance,
                                                                                 const Eigen::Vector3d& sourceToTargetDirectionLocalFrame );

    Eigen::Vector3d evaluateRadiationPressureForcePartialWrtSpecularReflectivity(
            double sourceIrradiance,
            const Eigen::Vector3d& sourceToTargetDirectionLocalFrame );

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& getBodyFixedPanels( )
    {
        return bodyFixedPanels_;
    }

    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >& getSegmentFixedPanels( )
    {
        return segmentFixedPanels_;
    }

    bool forceFunctionRequiresLocalFrameInputs( ) override
    {
        return true;
    }

    std::vector< Eigen::Vector3d >& getSurfaceNormals( )
    {
        return surfaceNormals_;
    }

    std::vector< double >& getSurfacePanelCosines( const std::string& sourceName )
    {
        if( surfacePanelCosinesPerSource_.count( sourceName ) == 0 )
        {
            throw std::runtime_error( "Error wen getting panelled radiation pressure target surface cosines from body " + sourceName +
                                      ", no such source is saved" );
        }
        return surfacePanelCosinesPerSource_[ sourceName ];
    }

    std::vector< Eigen::Vector3d >& getPanelForces( const std::string& sourceName )
    {
        if( panelForcesPerSource_.count( sourceName ) == 0 )
        {
            throw std::runtime_error( "Error wen getting panelled radiation pressure panel force from body " + sourceName +
                                      ", no such source is saved" );
        }
        return panelForcesPerSource_[ sourceName ];
    }

    const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& getAllPanels( ) const
    {
        return allPanels_;
    }

    void saveLocalComputations( const std::string sourceName, const bool saveCosines ) override;

    std::vector< double > getIlluminatedPanelFractions( const std::string& sourceName )
    {
        if( illuminatedPanelFractionsPerSource_.count( sourceName ) == 0 )
        {
            throw std::runtime_error( "Error when getting panelled radiation pressure illuminated panel fractions from body " + sourceName +
                                      ", no such source is saved" );
        }
        return illuminatedPanelFractionsPerSource_[ sourceName ];
    }

    std::vector< std::string > getPanelTypeIdList( ) const
    {
        return panelTypeIdList_;
    }

    std::map< std::string, std::shared_ptr< system_models::SelfShadowing > > getSelfShadowingPerSources( ) const
    {
        return selfShadowingPerSource_;
    }

    int getTotalNumberOfPanels( ) const
    {
        return totalNumberOfPanels_;
    }

private:
    void updateMembers_( double currentTime ) override;

    void resetDerivedComputations( const std::string sourceName ) override
    {
        for( unsigned int i = 0; i < panelForces_.size( ); i++ )
        {
            panelForces_.at( i ).setZero( );
        }
        panelForcesPerSource_[ sourceName ] = panelForces_;

        for( unsigned int i = 0; i < panelTorques_.size( ); i++ )
        {
            panelTorques_.at( i ).setZero( );
        }
        panelTorquesPerSource_[ sourceName ] = panelTorques_;

        for( unsigned int i = 0; i < surfacePanelCosines_.size( ); i++ )
        {
            surfacePanelCosines_.at( i ) = 0.0;
        }
        surfacePanelCosinesPerSource_[ sourceName ] = surfacePanelCosines_;

        for( unsigned int i = 0; i < illuminatedPanelFractions_.size( ); i++ )
        {
            illuminatedPanelFractions_.at( i ) = 0.0;
        }
        illuminatedPanelFractionsPerSource_[ sourceName ] = illuminatedPanelFractions_;
        for( auto it: selfShadowingPerSource_ )
        {
            it.second->reset( );
        }
    }

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > bodyFixedPanels_;

    std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > > segmentFixedPanels_;

    std::map< std::string, std::function< Eigen::Quaterniond( ) > > segmentFixedToBodyFixedRotations_;

    int totalNumberOfPanels_;

    std::vector< Eigen::Vector3d > surfaceNormals_;

    std::vector< Eigen::Vector3d > panelCentroidMomentArms_;

    // Source-specific variables
    std::vector< double > surfacePanelCosines_;
    std::vector< Eigen::Vector3d > panelForces_;
    std::vector< Eigen::Vector3d > panelTorques_;
    std::vector< double > illuminatedPanelFractions_;

    std::map< std::string, std::vector< double > > surfacePanelCosinesPerSource_;
    std::map< std::string, std::vector< Eigen::Vector3d > > panelForcesPerSource_;
    std::map< std::string, std::vector< Eigen::Vector3d > > panelTorquesPerSource_;
    std::map< std::string, std::vector< double > > illuminatedPanelFractionsPerSource_;

    // ssh variables
    std::map< std::string, int > maximumNumberOfPixelsPerSource_;
    std::map< std::string, std::shared_ptr< system_models::SelfShadowing > > selfShadowingPerSource_;
    std::vector< std::string > panelTypeIdList_;
    std::vector< double > unityIlluminationFraction_;
    const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& allPanels_;
    bool panelGeometryDefined_;
};

}  // namespace electromagnetism
}  // namespace tudat

#endif  // TUDAT_RADIATIONPRESSURETARGETMODEL_H
