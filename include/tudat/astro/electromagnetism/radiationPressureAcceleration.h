/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSUREACCELERATION_H
#define TUDAT_RADIATIONPRESSUREACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/electromagnetism/occultationModel.h"


namespace tudat
{
namespace electromagnetism
{


/*!
 * Class modeling radiation pressure acceleration. Radiation pressure accelerates a target due to electromagnetic
 * radiation from a source.
 */
class RadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:

    virtual ~RadiationPressureAcceleration( ) { }

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime) override;

    virtual std::shared_ptr<RadiationSourceModel> getSourceModel() const = 0;

    std::shared_ptr<RadiationPressureTargetModel> getTargetModel() const
    {
        return targetModel_;
    }

    std::shared_ptr<OccultationModel> getSourceToTargetOccultationModel() const
    {
        return sourceToTargetOccultationModel_;
    }

    double getReceivedIrradiance() const
    {
        return receivedIrradiance;
    }

    std::function<double()>  getTargetMassFunction()
    {
        return targetMassFunction_;
    }

    std::function<Eigen::Vector3d()> getSourcePositionFunction( )
    {
        return sourcePositionFunction_;
    }

    std::function<Eigen::Vector3d()> getTargetPositionFunction( )
    {
        return targetPositionFunction_;
    }

    std::function<Eigen::Quaterniond()> getTargetRotationFromLocalToGlobalFrameFunction( )
    {
        return targetRotationFromLocalToGlobalFrameFunction_;
    }


    double getCurrentRadiationPressure( )
    {
        return getReceivedIrradiance( ) / physical_constants::SPEED_OF_LIGHT;
    }

    void enableScaling( )
    {
        isScalingModelSet_ = true;
    }

    double getSourceDirectionScaling( )
    {
        return sourceDirectionScaling_;
    }
    
    void setSourceDirectionScaling( const double sourceDirectionScaling )
    {
        enableScaling( );
        sourceDirectionScaling_ = sourceDirectionScaling;
    }
    
    double getPerpendicularSourceDirectionScaling( )
    {
        return perpendicularSourceDirectionScaling_;
    }

    void setPerpendicularSourceDirectionScaling( const double perpendicularSourceDirectionScaling )
    {
        enableScaling( );
        perpendicularSourceDirectionScaling_ = perpendicularSourceDirectionScaling;
    }

    Eigen::Vector3d getTargetCenterPositionInSourceFrame( )
    {
        return targetCenterPositionInSourceFrame_;
    }

    Eigen::Vector3d getCurrentUnscaledAcceleration( )
    {
        return currentUnscaledAcceleration_;
    }



protected:
    RadiationPressureAcceleration(const std::function<Eigen::Vector3d()>& sourcePositionFunction,
                                  const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
                                  const std::function<Eigen::Vector3d()>& targetPositionFunction,
                                  const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
                                  const std::function<double()>& targetMassFunction,
                                  const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel,
                                  const std::string& sourceName ) :
            sourcePositionFunction_(sourcePositionFunction),
            targetModel_(targetModel),
            targetPositionFunction_(targetPositionFunction),
            targetRotationFromLocalToGlobalFrameFunction_(targetRotationFromLocalToGlobalFrameFunction),
            targetMassFunction_(targetMassFunction),
            sourceToTargetOccultationModel_(sourceToTargetOccultationModel),
            currentUnscaledAcceleration_( Eigen::Vector3d::Constant( TUDAT_NAN ) ),
            receivedIrradiance(TUDAT_NAN),
            isScalingModelSet_( false ),
            sourceDirectionScaling_( 1.0 ),
            perpendicularSourceDirectionScaling_( 1.0 ),
            sourceName_( sourceName )
            {}

    virtual void computeAcceleration( ) = 0;

    virtual void scaleRadiationPressureAcceleration( )
    {
        if( !isScalingModelSet_ )
        {
            currentAcceleration_ = currentUnscaledAcceleration_;
        }
        else
        {
            Eigen::Vector3d targetUnitVector = targetCenterPositionInSourceFrame_.normalized( );
            Eigen::Vector3d toTargetComponent =  currentUnscaledAcceleration_ - targetUnitVector.dot( currentUnscaledAcceleration_ ) * targetUnitVector;
            currentAcceleration_ = sourceDirectionScaling_ * toTargetComponent + perpendicularSourceDirectionScaling_ * ( currentUnscaledAcceleration_ - toTargetComponent );
        }
    }



    std::function<Eigen::Vector3d()> sourcePositionFunction_;
    std::shared_ptr<RadiationPressureTargetModel> targetModel_;
    std::function<Eigen::Vector3d()> targetPositionFunction_;
    std::function<Eigen::Quaterniond()> targetRotationFromLocalToGlobalFrameFunction_;
    std::function<double()> targetMassFunction_;
    std::shared_ptr<OccultationModel> sourceToTargetOccultationModel_;


    Eigen::Vector3d currentUnscaledAcceleration_;


    Eigen::Vector3d sourceCenterPositionInGlobalFrame_;
    Eigen::Vector3d targetCenterPositionInGlobalFrame_;
    Eigen::Vector3d targetCenterPositionInSourceFrame_;

    // For dependent variable
    double receivedIrradiance;
    double currentRadiationPressure_;

    bool isScalingModelSet_;
    double sourceDirectionScaling_;
    double perpendicularSourceDirectionScaling_;

    std::string sourceName_;

};

/*!
 * Class modeling radiation pressure acceleration from an isotropic point source.
 *
 * Assumptions:
 *  - Radiation from the source is evaluated at the target center (i.e. the target extent is neglected,
 *    which is acceptable, e.g., when evaluating solar radiation at Earth since rays are
 *    essentially parallel and the target extent vanishes compared to the 1 AU distance).
 */
class IsotropicPointSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with an isotropic point source.
    *
    * @param sourceModel Radiation source model of the source body
    * @param sourceBodyShapeModel Body shape model of the source body
    * @param sourcePositionFunction Position function of the source body
    * @param targetModel Radiation pressure target model of the target body
    * @param targetPositionFunction Position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the target body
    * @param targetMassFunction Mass function of the target body
    * @param sourceToTargetOccultationModel Occultation model between source and target
    */
    IsotropicPointSourceRadiationPressureAcceleration(
            const std::shared_ptr<IsotropicPointRadiationSourceModel>& sourceModel,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel) :
            RadiationPressureAcceleration(sourcePositionFunction, targetModel,
                                          targetPositionFunction, targetRotationFromLocalToGlobalFrameFunction,
                                          targetMassFunction, sourceToTargetOccultationModel,
                                          sourceModel->getSourceName( ) ),
                                          sourceModel_(sourceModel),
                                          sourceBodyShapeModel_(sourceBodyShapeModel),
            sourceToTargetReceivedFraction(TUDAT_NAN) {}

    ~IsotropicPointSourceRadiationPressureAcceleration( ){ }

    std::shared_ptr<RadiationSourceModel> getSourceModel() const override
    {
        return sourceModel_;
    }

    double getSourceToTargetReceivedFraction() const
    {
        return sourceToTargetReceivedFraction;
    }

    Eigen::Vector3d getSourcePosition( )
    {
        return sourceCenterPositionInGlobalFrame_;
    }

    Eigen::Vector3d getTargetPosition( )
    {
        return targetCenterPositionInGlobalFrame_;
    }

    Eigen::Vector3d getTargetPositionWrtSource( )
    {
        return targetCenterPositionInSourceFrame_;
    }

    double getCurrentTargetMass( )
    {
        return currentTargetMass_;
    }

    Eigen::Quaterniond getTargetRotationFromLocalToGlobalFrame( )
    {
        return targetRotationFromLocalToGlobalFrame_;
    }

    Eigen::Quaterniond getTargetRotationFromGlobalToLocalFrame( )
    {
        return targetRotationFromGlobalToLocalFrame_;
    }

private:
    void computeAcceleration() override;

    std::shared_ptr<IsotropicPointRadiationSourceModel> sourceModel_;

    std::shared_ptr<basic_astrodynamics::BodyShapeModel> sourceBodyShapeModel_;

    double currentTargetMass_;

    // For dependent variable
    double sourceToTargetReceivedFraction;

    Eigen::Quaterniond targetRotationFromLocalToGlobalFrame_;
    Eigen::Quaterniond targetRotationFromGlobalToLocalFrame_;
};

/*!
 * Class modeling radiation pressure acceleration from a paneled source. An original source illuminates the
 * source, which emits the radiation that accelerates the target. Only isotropic point sources are supported as
 * original sources for now. For occultation, panel-center-to-target-center visibility is used for each panel instead
 * of a shadow function.
 *
 * Assumptions:
 *  - Radiation from the original source is evaluated at the source center (i.e. the source extent is
 *    neglected, which is acceptable, e.g., when evaluating solar radiation at Earth for albedo). Occultation using
 *    the shadow function is done for each panel at the panel center, however.
 *  - Radiation from each panel/subsource of the paneled source is evaluated at the target center (i.e. the target
 *    extent is neglected). This also applies for occultation calculations.
 */
 // If other original source types are to be supported in the future, add the rotational state update to
 // createEnvironmentUpdater.cpp under case radiation_pressure
class PaneledSourceRadiationPressureAcceleration : public RadiationPressureAcceleration
{
public:
    /*!
    * Creates a radiation pressure acceleration model with a paneled source.
    *
    * @param sourceModel Radiation source model of the source body
    * @param sourcePositionFunction Position function of the source body
    * @param sourceRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the source body
    * @param targetModel Radiation pressure target model of the target body
    * @param targetPositionFunction Position function of the target body
    * @param targetRotationFromLocalToGlobalFrameFunction Local-to-global rotation function of the target body
    * @param targetMassFunction Mass function of the target body
    * @param sourceToTargetOccultationModel Occultation model between source and target
    */
    PaneledSourceRadiationPressureAcceleration(
            const std::shared_ptr<PaneledRadiationSourceModel>& sourceModel,
            const std::function<Eigen::Vector3d()>& sourcePositionFunction,
            const std::function<Eigen::Quaterniond()>& sourceRotationFromLocalToGlobalFrameFunction,
            const std::shared_ptr<RadiationPressureTargetModel>& targetModel,
            const std::function<Eigen::Vector3d()>& targetPositionFunction,
            const std::function<Eigen::Quaterniond()>& targetRotationFromLocalToGlobalFrameFunction,
            const std::function<double()>& targetMassFunction,
            const std::shared_ptr<OccultationModel>& sourceToTargetOccultationModel) :
            RadiationPressureAcceleration(
                    sourcePositionFunction,
                    targetModel,
                    targetPositionFunction,
                    targetRotationFromLocalToGlobalFrameFunction,
                    targetMassFunction,
                    sourceToTargetOccultationModel,
                    sourceModel->getSourceName( ) ),
            sourceModel_(sourceModel),
            sourceRotationFromLocalToGlobalFrameFunction_(sourceRotationFromLocalToGlobalFrameFunction),
            visibleAndEmittingSourcePanelCount(-1) {}

    ~PaneledSourceRadiationPressureAcceleration( ) { }

    std::shared_ptr<RadiationSourceModel> getSourceModel() const override
    {
        return sourceModel_;
    }

    std::shared_ptr<PaneledRadiationSourceModel> getPaneledSourceModel( )
    {
        return sourceModel_;
    }


    unsigned int getVisibleAndEmittingSourcePanelCount() const
    {
        return visibleAndEmittingSourcePanelCount;
    }

    /*
     * Visible area of a paneled source, calculated from sum of areas of visible panels (including non-emitting panels).
     */
    double getVisibleSourceArea() const
    {
        return sourceModel_->getVisibleArea();
    }


    void enableSavePanellingGeometry( )
    {
        savePanellingGeometry_ = true;
        savedPanelIrradiances_.resize( sourceModel_->getNumberOfPanels( ) );
    }

    void enableSavePanellingIrradiance( )
    {
        savePanellingIrradiance_ = true;
        savedPanelGeometries_.resize( sourceModel_->getNumberOfPanels( ) );
    }

    std::vector< double >& getSavedPanelIrradiances( )
    {
        return savedPanelIrradiances_;
    }

    std::vector< Eigen::Vector7d >& getSavedPanelGeometries( )
    {
        return savedPanelGeometries_;
    }

private:
    void computeAcceleration() override;

    std::shared_ptr<PaneledRadiationSourceModel> sourceModel_;
    std::function<Eigen::Quaterniond()> sourceRotationFromLocalToGlobalFrameFunction_;

    // For dependent variable
    unsigned int visibleAndEmittingSourcePanelCount;

    bool savePanellingGeometry_ = false;

    bool savePanellingIrradiance_ = false;

    std::vector< double > savedPanelIrradiances_;

    std::vector< Eigen::Vector7d > savedPanelGeometries_;


};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSUREACCELERATION_H
