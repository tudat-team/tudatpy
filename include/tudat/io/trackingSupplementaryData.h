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
#include <vector>
#include "tudat/basics/basicTypedefs.h"

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

        FrequencyRamp( const double startTime, const double endTime, const double startFrequency, const double frequencyRate ):
            startTime_( startTime ), endTime_( endTime ), startFrequency_( startFrequency ), frequencyRate_( frequencyRate )
        {}

        double startTime_ = 0.0;
        double endTime_ = 0.0;
        double startFrequency_ = 0.0;
        double frequencyRate_ = 0.0;
    };

    RampedFrequencySupplementaryData( ): FrequencySupplementaryData( "ramped_frequency" ) {}

    RampedFrequencySupplementaryData( const std::vector< FrequencyRamp >& frequencyRamps ):
        FrequencySupplementaryData( "ramped_frequency" ), frequencyRamps_( frequencyRamps )
    {}

    void addFrequencyRamp( const double startTime, const double endTime, const double startFrequency, const double frequencyRate )
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

private:
    std::string cameraId_;

    double focalLength_ = 0.0;

    Eigen::Vector2d principalPoint_ = Eigen::Vector2d::Zero( );

    Eigen::Vector4d fieldOfViewBounds_ = Eigen::Vector4d::Zero( );

    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Zero( );

    Eigen::Matrix< double, 6, 1 > distortionCoefficients_ = Eigen::Matrix< double, 6, 1 >::Zero( );

    Eigen::Vector3d mountingOffsets_ = Eigen::Vector3d::Zero( );

    std::map< double, Eigen::Quaterniond > rotationFromInertialToCameraFrameHistory_;
};

class TranslationalStateSupplementaryData
{
public:
    TranslationalStateSupplementaryData( ) = default;

    TranslationalStateSupplementaryData( const std::map< double, Eigen::Vector6d >& stateHistory,
                                         const std::string& frameOrigin,
                                         const bool isVelocityDefined ):
        stateHistory_( stateHistory ), frameOrigin_( frameOrigin ), isVelocityDefined_( isVelocityDefined )
    {}

    const std::map< double, Eigen::Vector6d >& getStateHistory( ) const
    {
        return stateHistory_;
    }

    const std::string& getFrameOrigin( ) const
    {
        return frameOrigin_;
    }

    bool isVelocityDefined( ) const
    {
        return isVelocityDefined_;
    }

private:
    std::map< double, Eigen::Vector6d > stateHistory_;

    std::string frameOrigin_;

    bool isVelocityDefined_ = false;
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

class TrackingSupplementaryData
{
public:
    TrackingSupplementaryData( ) = default;

    TrackingSupplementaryData( const std::string& bodyName, const std::string& referencePointName ):
        bodyName_( bodyName ), referencePointName_( referencePointName )
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

private:
    TranslationalStateSupplementaryData translationalStateSupplementaryData_;

    RotationalStateSupplementaryData rotationalStateSupplementaryData_;

    std::vector< std::shared_ptr< FrequencySupplementaryData > > frequencySupplementaryData_;

    std::vector< std::shared_ptr< InstrumentSupplementaryData > > instrumentSupplementaryData_;

    std::string bodyName_;

    std::string referencePointName_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_SUPPLEMENTARY_DATA_H
