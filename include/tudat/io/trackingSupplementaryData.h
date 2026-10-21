/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_SUPPLEMENTARY_DATA_H
#define TUDAT_TRACKING_SUPPLEMENTARY_DATA_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"

namespace tudat
{

namespace data
{

class FrequencySupplementaryData
{
public:
    FrequencySupplementaryData( const std::string& frequencySupplementaryDataKind ):
        frequencySupplementaryDataKind_( frequencySupplementaryDataKind )
    {}

    virtual ~FrequencySupplementaryData( ) = default;

    std::string getFrequencySupplementaryDataKind( ) const
    {
        return frequencySupplementaryDataKind_;
    }

private:
    std::string frequencySupplementaryDataKind_;
};

class RampedFrequencySupplementaryData : public FrequencySupplementaryData
{
public:
    struct FrequencyRamp {
        FrequencyRamp( ) = default;

        FrequencyRamp( const Time startTime, const Time endTime, const double startFrequency, const double frequencyRate ):
            startTime_( startTime ), endTime_( endTime ), startFrequency_( startFrequency ), frequencyRate_( frequencyRate )
        {}

        Time startTime_ = 0.0;
        Time endTime_ = 0.0;
        double startFrequency_ = 0.0;
        double frequencyRate_ = 0.0;
    };

    RampedFrequencySupplementaryData( ): FrequencySupplementaryData( "ramped_frequency" ) {}

    RampedFrequencySupplementaryData( const std::vector< FrequencyRamp >& frequencyRamps ):
        FrequencySupplementaryData( "ramped_frequency" ), frequencyRamps_( frequencyRamps )
    {}

    void addFrequencyRamp( const Time startTime, const Time endTime, const double startFrequency, const double frequencyRate )
    {
        frequencyRamps_.emplace_back( startTime, endTime, startFrequency, frequencyRate );
    }

    const std::vector< FrequencyRamp >& getFrequencyRamps( ) const
    {
        return frequencyRamps_;
    }

private:
    std::vector< FrequencyRamp > frequencyRamps_;
};

class PiecewiseConstantFrequencySupplementaryData : public FrequencySupplementaryData
{
public:
    PiecewiseConstantFrequencySupplementaryData( ): FrequencySupplementaryData( "piecewise_constant_frequency" ) {}

    PiecewiseConstantFrequencySupplementaryData( const std::map< double, double >& frequencyHistory ):
        FrequencySupplementaryData( "piecewise_constant_frequency" ), frequencyHistory_( frequencyHistory )
    {}

    void setFrequency( const double time, const double frequency )
    {
        frequencyHistory_[ time ] = frequency;
    }

    const std::map< double, double >& getFrequencyHistory( ) const
    {
        return frequencyHistory_;
    }

private:
    std::map< double, double > frequencyHistory_;
};

class InstrumentSupplementaryData
{
public:
    explicit InstrumentSupplementaryData( const std::string& instrumentSupplementaryDataKind ):
        instrumentSupplementaryDataKind_( instrumentSupplementaryDataKind )
    {}

    virtual ~InstrumentSupplementaryData( ) = default;

    std::string getInstrumentSupplementaryDataKind( ) const
    {
        return instrumentSupplementaryDataKind_;
    }

private:
    std::string instrumentSupplementaryDataKind_;
};

//! Calibration of one picture, for a camera whose calibration changes between pictures (for instance by pixel binning).
struct CameraPictureCalibration {
    double focalLength_ = 0.0;
    Eigen::Vector2d principalPoint_ = Eigen::Vector2d::Zero( );
    Eigen::Vector4d fieldOfViewBounds_ = Eigen::Vector4d::Zero( );
    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Zero( );

    bool operator==( const CameraPictureCalibration& other ) const
    {
        return focalLength_ == other.focalLength_ && principalPoint_ == other.principalPoint_ &&
                fieldOfViewBounds_ == other.fieldOfViewBounds_ && kMatrix_ == other.kMatrix_;
    }
};

class CameraInstrumentSupplementaryData : public InstrumentSupplementaryData
{
public:
    CameraInstrumentSupplementaryData( ): InstrumentSupplementaryData( "camera_settings" ) {}

    CameraInstrumentSupplementaryData( const std::string& cameraId,
                                       const double focalLength,
                                       const Eigen::Vector2d& principalPoint,
                                       const Eigen::Vector4d& fieldOfViewBounds,
                                       const Eigen::Matrix< double, 2, 3 >& kMatrix,
                                       const Eigen::Matrix< double, 6, 1 >& distortionCoefficients,
                                       const Eigen::Vector3d& mountingOffsets,
                                       const std::map< double, Eigen::Quaterniond >& rotationFromInertialToCameraFrameHistory =
                                               std::map< double, Eigen::Quaterniond >( ) ):
        InstrumentSupplementaryData( "camera_settings" ), cameraId_( cameraId ), focalLength_( focalLength ),
        principalPoint_( principalPoint ), fieldOfViewBounds_( fieldOfViewBounds ), kMatrix_( kMatrix ),
        distortionCoefficients_( distortionCoefficients ), mountingOffsets_( mountingOffsets ),
        rotationFromInertialToCameraFrameHistory_( rotationFromInertialToCameraFrameHistory )
    {}

    const std::string& getCameraId( ) const
    {
        return cameraId_;
    }

    double getFocalLength( ) const
    {
        return focalLength_;
    }

    const Eigen::Vector2d& getPrincipalPoint( ) const
    {
        return principalPoint_;
    }

    const Eigen::Vector4d& getFieldOfViewBounds( ) const
    {
        return fieldOfViewBounds_;
    }

    const Eigen::Matrix< double, 2, 3 >& getKMatrix( ) const
    {
        return kMatrix_;
    }

    const Eigen::Matrix< double, 6, 1 >& getDistortionCoefficients( ) const
    {
        return distortionCoefficients_;
    }

    const Eigen::Vector3d& getMountingOffsets( ) const
    {
        return mountingOffsets_;
    }

    const std::map< double, Eigen::Quaterniond >& getRotationFromInertialToCameraFrameHistory( ) const
    {
        return rotationFromInertialToCameraFrameHistory_;
    }

    //! Body whose body-fixed frame the pointing history starts from.
    /*!
     * When empty, the pointing history holds rotations from the inertial frame to the camera frame. When set, it holds
     * rotations from the body-fixed frame of this body to the camera frame. The body's rotation model then completes
     * the rotation from the inertial frame when the camera is created.
     */
    const std::string& getPointingReferenceBodyName( ) const
    {
        return pointingReferenceBodyName_;
    }

    void setPointingReferenceBodyName( const std::string& pointingReferenceBodyName )
    {
        pointingReferenceBodyName_ = pointingReferenceBodyName;
    }

    //! Calibration of each picture, per picture epoch of the pointing history. When empty, all pictures use the calibration of
    //! the camera.
    const std::map< double, CameraPictureCalibration >& getPictureCalibrationHistory( ) const
    {
        return pictureCalibrationHistory_;
    }

    void setPictureCalibrationHistory( const std::map< double, CameraPictureCalibration >& pictureCalibrationHistory )
    {
        pictureCalibrationHistory_ = pictureCalibrationHistory;
    }

    //! Identifier of each picture, per picture epoch of the pointing history. Empty when not provided.
    const std::map< double, std::string >& getPictureIdHistory( ) const
    {
        return pictureIdHistory_;
    }

    void setPictureIdHistory( const std::map< double, std::string >& pictureIdHistory )
    {
        pictureIdHistory_ = pictureIdHistory;
    }

    //! A-priori uncertainty of the camera pointing about the three camera axes [rad], per picture epoch of the pointing history.
    //! Empty when not provided.
    const std::map< double, Eigen::Vector3d >& getPointingUncertaintyHistory( ) const
    {
        return pointingUncertaintyHistory_;
    }

    void setPointingUncertaintyHistory( const std::map< double, Eigen::Vector3d >& pointingUncertaintyHistory )
    {
        pointingUncertaintyHistory_ = pointingUncertaintyHistory;
    }

    //! Merge compatible camera histories; conflicting calibrations, attitudes, picture identifiers or pointing uncertainties are
    //! errors.
    void merge( const CameraInstrumentSupplementaryData& other )
    {
        if( cameraId_ != other.cameraId_ || focalLength_ != other.focalLength_ || !principalPoint_.isApprox( other.principalPoint_ ) ||
            !fieldOfViewBounds_.isApprox( other.fieldOfViewBounds_ ) || !kMatrix_.isApprox( other.kMatrix_ ) ||
            !distortionCoefficients_.isApprox( other.distortionCoefficients_ ) || !mountingOffsets_.isApprox( other.mountingOffsets_ ) ||
            pointingReferenceBodyName_ != other.pointingReferenceBodyName_ )
        {
            throw std::runtime_error( "Inconsistent calibration for PSF camera " + cameraId_ );
        }
        auto mergedHistory = rotationFromInertialToCameraFrameHistory_;
        for( const auto& entry : other.rotationFromInertialToCameraFrameHistory_ )
        {
            const auto inserted = mergedHistory.emplace( entry );
            if( !inserted.second && !inserted.first->second.toRotationMatrix( ).isApprox( entry.second.toRotationMatrix( ) ) )
            {
                throw std::runtime_error( "Conflicting pointing epochs for PSF camera " + cameraId_ );
            }
        }
        auto mergedPictureCalibrationHistory = pictureCalibrationHistory_;
        for( const auto& entry : other.pictureCalibrationHistory_ )
        {
            const auto inserted = mergedPictureCalibrationHistory.emplace( entry );
            if( !inserted.second && !( inserted.first->second == entry.second ) )
            {
                throw std::runtime_error( "Conflicting picture calibrations for PSF camera " + cameraId_ );
            }
        }
        auto mergedPictureIdHistory = pictureIdHistory_;
        for( const auto& entry : other.pictureIdHistory_ )
        {
            const auto inserted = mergedPictureIdHistory.emplace( entry );
            if( !inserted.second && inserted.first->second != entry.second )
            {
                throw std::runtime_error( "Conflicting picture identifiers for PSF camera " + cameraId_ );
            }
        }
        auto mergedUncertaintyHistory = pointingUncertaintyHistory_;
        for( const auto& entry : other.pointingUncertaintyHistory_ )
        {
            const auto inserted = mergedUncertaintyHistory.emplace( entry );
            if( !inserted.second && inserted.first->second != entry.second )
            {
                throw std::runtime_error( "Conflicting pointing uncertainties for PSF camera " + cameraId_ );
            }
        }
        rotationFromInertialToCameraFrameHistory_.swap( mergedHistory );
        pictureCalibrationHistory_.swap( mergedPictureCalibrationHistory );
        pictureIdHistory_.swap( mergedPictureIdHistory );
        pointingUncertaintyHistory_.swap( mergedUncertaintyHistory );
    }

private:
    std::string cameraId_;

    double focalLength_ = 0.0;

    Eigen::Vector2d principalPoint_ = Eigen::Vector2d::Zero( );

    Eigen::Vector4d fieldOfViewBounds_ = Eigen::Vector4d::Zero( );

    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Zero( );

    Eigen::Matrix< double, 6, 1 > distortionCoefficients_ = Eigen::Matrix< double, 6, 1 >::Zero( );

    Eigen::Vector3d mountingOffsets_ = Eigen::Vector3d::Zero( );

    std::map< double, Eigen::Quaterniond > rotationFromInertialToCameraFrameHistory_;

    std::string pointingReferenceBodyName_;

    std::map< double, CameraPictureCalibration > pictureCalibrationHistory_;

    std::map< double, std::string > pictureIdHistory_;

    std::map< double, Eigen::Vector3d > pointingUncertaintyHistory_;
};

class TranslationalStateSupplementaryData
{
public:
    TranslationalStateSupplementaryData( ) = default;

    //! Constructor for supplementary body-state histories loaded with tracking data.
    /*!
     * The epoch keys in the state history are tagged with a time scale. The default
     * is TDB to preserve the existing ephemeris convention for callers that do not
     * explicitly provide a time scale. The default frame orientation is J2000,
     * matching the tracking-data readers that currently create this object.
     */
    TranslationalStateSupplementaryData( const std::map< double, Eigen::Vector6d >& stateHistory,
                                         const std::string& frameOrigin,
                                         const bool isVelocityDefined,
                                         const std::string& timeScale = "TDB",
                                         const std::string& frameOrientation = "J2000" ):
        stateHistory_( stateHistory ), frameOrigin_( frameOrigin ), frameOrientation_( frameOrientation ),
        isVelocityDefined_( isVelocityDefined ), timeScale_( timeScale )
    {}

    const std::map< double, Eigen::Vector6d >& getStateHistory( ) const
    {
        return stateHistory_;
    }

    const std::string& getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

    const std::string& getFrameOrientation( ) const
    {
        return frameOrientation_;
    }

    bool isVelocityDefined( ) const
    {
        return isVelocityDefined_;
    }

    //! Time scale used by the state-history epoch keys.
    const std::string& getTimeScale( ) const
    {
        return timeScale_;
    }

private:
    std::map< double, Eigen::Vector6d > stateHistory_;

    std::string frameOrigin_;

    std::string frameOrientation_ = "J2000";

    bool isVelocityDefined_ = false;

    std::string timeScale_ = "TDB";
};

class RotationalStateSupplementaryData
{
public:
    const std::map< double, Eigen::Vector7d >& getRotationalStateHistory( ) const
    {
        return rotationalStateHistory_;
    }

    const std::string& getBaseFrameOrientation( ) const
    {
        return baseFrameOrientation_;
    }

    bool isAngularVelocityDefined( ) const
    {
        return isAngularVelocityDefined_;
    }

private:
    std::map< double, Eigen::Vector7d > rotationalStateHistory_;

    std::string baseFrameOrientation_;

    bool isAngularVelocityDefined_ = false;
};

//! Position of a reference point that is fixed to a body, such as a ground station or a surface landmark.
class ReferencePointSupplementaryData
{
public:
    explicit ReferencePointSupplementaryData( const Eigen::Vector3d& bodyFixedPosition ): bodyFixedPosition_( bodyFixedPosition ) {}

    //! Cartesian position of the reference point in the body-fixed frame [m].
    const Eigen::Vector3d& getBodyFixedPosition( ) const
    {
        return bodyFixedPosition_;
    }

private:
    Eigen::Vector3d bodyFixedPosition_;
};

class TrackingSupplementaryData
{
public:
    TrackingSupplementaryData( ) = default;

    TrackingSupplementaryData( const std::string& bodyName, const std::string& referencePointName ):
        bodyName_( bodyName ), referencePointName_( referencePointName ), isPassiveRadarReflector_( false )
    {}

    void setTranslationalStateSupplementaryData( const TranslationalStateSupplementaryData& translationalStateSupplementaryData )
    {
        translationalStateSupplementaryData_ = translationalStateSupplementaryData;
    }

    const TranslationalStateSupplementaryData& getTranslationalStateSupplementaryData( ) const
    {
        return translationalStateSupplementaryData_;
    }

    void setRotationalStateSupplementaryData( const RotationalStateSupplementaryData& rotationalStateSupplementaryData )
    {
        rotationalStateSupplementaryData_ = rotationalStateSupplementaryData;
    }

    const RotationalStateSupplementaryData& getRotationalStateSupplementaryData( ) const
    {
        return rotationalStateSupplementaryData_;
    }

    void setFrequencySupplementaryData( const std::vector< std::shared_ptr< FrequencySupplementaryData > >& frequencySupplementaryData )
    {
        frequencySupplementaryData_ = frequencySupplementaryData;
    }

    const std::vector< std::shared_ptr< FrequencySupplementaryData > >& getFrequencySupplementaryData( ) const
    {
        return frequencySupplementaryData_;
    }

    void setInstrumentSupplementaryData( const std::vector< std::shared_ptr< InstrumentSupplementaryData > >& instrumentSupplementaryData )
    {
        instrumentSupplementaryData_ = instrumentSupplementaryData;
    }

    const std::vector< std::shared_ptr< InstrumentSupplementaryData > >& getInstrumentSupplementaryData( ) const
    {
        return instrumentSupplementaryData_;
    }

    void setBodyName( const std::string& bodyName )
    {
        bodyName_ = bodyName;
    }

    const std::string& getBodyName( ) const
    {
        return bodyName_;
    }

    void setReferencePointName( const std::string& referencePointName )
    {
        referencePointName_ = referencePointName;
    }

    const std::string& getReferencePointName( ) const
    {
        return referencePointName_;
    }

    //! Body-fixed position of the reference point. Null when the reference point position is not provided.
    void setReferencePointSupplementaryData( const std::shared_ptr< ReferencePointSupplementaryData >& referencePointSupplementaryData )
    {
        referencePointSupplementaryData_ = referencePointSupplementaryData;
    }

    const std::shared_ptr< ReferencePointSupplementaryData >& getReferencePointSupplementaryData( ) const
    {
        return referencePointSupplementaryData_;
    }

    bool isPassiveRadarReflector( ) const
    {
        return isPassiveRadarReflector_;
    }

    void setIsPassiveRadarReflector( const bool isPassiveRadarReflector )
    {
        isPassiveRadarReflector_ = isPassiveRadarReflector;
    }

private:
    TranslationalStateSupplementaryData translationalStateSupplementaryData_;

    RotationalStateSupplementaryData rotationalStateSupplementaryData_;

    std::vector< std::shared_ptr< FrequencySupplementaryData > > frequencySupplementaryData_;

    std::vector< std::shared_ptr< InstrumentSupplementaryData > > instrumentSupplementaryData_;

    std::shared_ptr< ReferencePointSupplementaryData > referencePointSupplementaryData_;

    std::string bodyName_;

    std::string referencePointName_;

    bool isPassiveRadarReflector_ = false;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_SUPPLEMENTARY_DATA_H
