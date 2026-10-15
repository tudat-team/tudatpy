/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READ_SUM_LMK_FILES_H
#define TUDAT_READ_SUM_LMK_FILES_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/timeType.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace input_output
{

namespace sum_lmk
{

//! ---------------------------------------------------------------------------------------------------
//! SPC (stereophotoclinometry) SUM and LMK file contents
//! ---------------------------------------------------------------------------------------------------
//!
//! SPC is R. Gaskell's shape-and-landmark software, hosted by PSI. There is no formal published record
//! specification for the ASCII SUM/LMK files; the field meanings below are assembled from:
//!
//!   [WIKI]  The SPC Wiki, PSI: https://web.psi.edu/spc_wiki/  In particular the MAPFILES page, which
//!           documents the binary maplet header and so defines the maplet local frame shared with LMK.
//!   [GIANT] NASA Goddard's GIANT optical-navigation library, whose Summary and Landmark classes read
//!           these same files and document each record:
//!           https://aliounis.github.io/giant_documentation/utilities/giant.utilities.stereophotoclinometry.html
//!   [MEAS]  Measured in this repository against the archived 67P dataset in tests/data/sum_lmk, and
//!           locked in by the unit tests named alongside each field.
//!
//! Every field is tagged with how well its meaning is established. Please keep these tags honest: a
//! wrong assumption here shows up as a silent few-pixel bias, not as a failure.
//!
//!   [CERTAIN]    Stated by [WIKI] or [GIANT], and where it affects geometry also confirmed by [MEAS].
//!   [UNVERIFIED] Plausible reading, consistent with the data, but not stated by any source we have.
//!   [UNKNOWN]    Meaning not established. Parsed and carried, never interpreted. Do not use these
//!                without first establishing what they are.
//!
//! Units: the SPC files express all lengths in KILOMETRES. Where this reader converts to metres it is
//! noted per field, because the conversion is NOT applied uniformly - some fields are deliberately left
//! in kilometres, and mixing them silently costs a factor of 1000.
//!
//! Records present in the files but NOT represented in these structs:
//!   SUM "LIMB FITS"    - parsed into limbFitObservations_, but never converted to observations.
//!   LMK "MAP OVERLAPS" - skipped entirely. [WIKI] describes it as "the relative locations of
//!                        overlapping maplets", i.e. SPC's internal maplet-to-maplet ties. Of no use
//!                        without the maplets themselves, which live in the binary .MAP files this
//!                        codebase does not read.
//!   LMK "PICTURES"     - parsed, but redundant: the SUM side already carries the same measurements.
//!
//! The archived 67P fixture in tests/data/sum_lmk/real_67p has its LMK PICTURES sections trimmed out and
//! contains no SUM LIMB FITS sections, so neither is covered by the real-data tests.

//! Tolerance for the orthonormality and determinant checks on a SUM CX/CY/CZ attitude matrix.
//!
//! The SPC archive contains valid matrices with a Frobenius orthogonality error just below 8e-8,
//! due to its stored numerical precision. This tolerance admits those matrices while still
//! rejecting matrices that are materially non-rotational.
constexpr double sumCameraRotationMatrixTolerance = 1.0E-7;

//! Tolerance for the orthonormality and determinant checks on an LMK UX/UY/UZ landmark frame.
//!
//! Deliberately looser than sumCameraRotationMatrixTolerance. The worst orthonormality error across
//! the archived LMK files is 2.3e-7, so reusing the camera tolerance here would reject the majority
//! of real landmark files while still admitting nothing that a 1e-6 check would not.
constexpr double lmkLandmarkFrameTolerance = 1.0E-6;

//! One row of a SUM file's LANDMARKS section: where a named landmark was measured in this image.
//! [GIANT] Summary.landmarks_, "dictionary mapping landmark names to pixel locations (x, y)".
struct SumLandmarkObservation {
    //! [CERTAIN] Landmark name, matching the LMK file of the same name. [GIANT]: 6 characters, typically
    //! alphanumeric.
    std::string landmarkId_;

    //! [CERTAIN] Measured landmark centre in the image, (sample/column, line/row), in pixels. THESE ARE
    //! THE OBSERVATIONS the pixel_coordinates observable fits.
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );
};

//! One row of a SUM file's LIMB FITS section: a point measured on the body's limb rather than on a landmark.
//! [GIANT] Summary.limb_fits_, "dictionary mapping limb landmark names to pixel locations and
//! uncertainty (x, y, sigma)". Parsed but NOT converted to observations by this codebase.
struct SumLimbFitObservation {
    //! [CERTAIN] Name of the limb feature.
    std::string featureId_;

    //! [CERTAIN] Measured limb point, (sample/column, line/row), in pixels.
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );

    //! [UNVERIFIED] Per-point uncertainty on the limb measurement. [GIANT] confirms a third value named
    //! "sigma" but gives no unit; pixels is the only reading consistent with the neighbouring columns.
    double sigma_ = TUDAT_NAN;
};

//! Parsed contents of one SUM file, i.e. SPC's camera calibration and spacecraft state solution for one
//! image, plus the pixel measurements made in it.
struct SumImageData {
    //! [CERTAIN] Image name, first line of the file. Becomes the Camera_<imageId> reference point.
    std::string imageId_;

    //! [CERTAIN] UTC exposure time, second line of the file. [GIANT] "UTC exposure time".
    std::string utcEpochString_;

    //! [CERTAIN] NPX, NLN from the "NPX, NLN, THRSH" row: (number of columns, number of rows), pixels.
    Eigen::Vector2i imageSize_ = Eigen::Vector2i::Zero( );

    //! [CERTAIN] First THRSH value. [GIANT] min_illum: "Pixels with DN values below this threshold are
    //! ignored". Used by SPC's own correlation, not by the observation model.
    int threshold_ = 0;

    //! [CERTAIN] Second THRSH value. [GIANT] max_illum: "Pixels with DN values above this threshold are
    //! ignored", i.e. the saturation cut. Note the row is labelled "NPX, NLN, THRSH" but carries FOUR
    //! values, because THRSH is a (min, max) pair.
    int maxDn_ = 0;

    //! [CERTAIN] MMFL: camera focal length in MILLIMETRES, left unconverted. [GIANT] focal_length,
    //! "Camera focal length in millimeters".
    double focalLengthMm_ = TUDAT_NAN;

    //! [CERTAIN] CTR: principal point, (sample/column, line/row), in pixels. [GIANT] princ_point.
    Eigen::Vector2d opticalCenter_ = Eigen::Vector2d::Zero( );

    //! [CERTAIN] SCOBJ: vector from the SPACECRAFT TO THE TARGET CENTRE, expressed in the target
    //! body-fixed frame. Converted from km to METRES here. [GIANT] position_camera_to_target, "Vector
    //! from spacecraft to target in target body fixed frame".
    //!
    //! MIND THE SIGN. The spacecraft position relative to the target centre is -spacecraftObjectVector_,
    //! and a landmark's position relative to the spacecraft is VLM + SCOBJ. This direction was first
    //! established by [MEAS] (getting it backwards moved the reprojection by tens of pixels) and is now
    //! confirmed by [GIANT]. Guarded by testRealSumLmkReprojection.
    Eigen::Vector3d spacecraftObjectVector_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] CX, CY, CZ as the three ROWS: the rotation matrix from the target body-fixed frame to
    //! the camera frame. [GIANT] rotation_target_fixed_to_camera. Camera +Z is the boresight.
    //! Orthonormality is checked against sumCameraRotationMatrixTolerance.
    Eigen::Matrix3d cameraAxes_ = Eigen::Matrix3d::Constant( TUDAT_NAN );

    //! [CERTAIN] SZ: unit vector from the TARGET TO THE SUN, in the target body-fixed frame. [GIANT]
    //! direction_target_to_sun. Parsed but not used by the observation model.
    Eigen::Vector3d sunDirectionBodyFixed_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] K-MATRIX: the 2x3 camera intrinsic matrix. [GIANT] intrinsic_matrix. In the projection
    //! this codebase implements (PsfCameraProjectionModel), pixel/line = K * [x, y, x*y] + CTR with
    //! x = f * X/Z, y = f * Y/Z in the camera frame. Confirmed by [MEAS]: reprojects the archived
    //! measurements to ~1 px (testRealSumLmkReprojection).
    Eigen::Matrix< double, 2, 3 > kMatrix_ = Eigen::Matrix< double, 2, 3 >::Constant( TUDAT_NAN );

    //! [CERTAIN] Raw four-value DISTORTION row. [GIANT] near_dist_params, "Deprecated NEAR distortion
    //! parameters". SPC SUM files conventionally keep these zero and supply an actual Owen distortion
    //! model separately in INIT_LITHOS, so the converter REJECTS non-zero values rather than silently
    //! ignoring a distortion it does not model.
    Eigen::Vector4d distortionCoefficients_ = Eigen::Vector4d::Zero( );

    //! [CERTAIN] SIGMA_VSO: SPC's 1-sigma formal uncertainty on SCOBJ. [GIANT]
    //! sig_pos_camera_to_target. Converted from km to METRES here, on the assumption that it shares
    //! SCOBJ's unit - which no source states outright, though the archived values (21-34 m) are only
    //! sensible that way. One value in the file is expanded to all three components.
    //!
    //! NOT used as an observation, and it should not be: SCOBJ is SPC's own fit to the same pixels we
    //! fit, so feeding it back in would double-count the data.
    Eigen::Vector3d spacecraftObjectSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [UNVERIFIED unit] SIGMA_PTG: SPC's 1-sigma formal uncertainty on the camera pointing. [GIANT]
    //! sig_rot_target_fixed_to_camera, "Uncertainty on pointing rotation" - the quantity is documented,
    //! the UNIT is not. RADIANS is inferred from magnitude alone: the archived values are 1.0-1.9e-4,
    //! i.e. 21-40 arcsec or 5-10 px at the OSIRIS NAC focal length, which is the right size for
    //! thermoelastic pointing variation. Left unconverted, and consumed as radians when building the
    //! a-priori for the camera_pointing_correction parameter.
    Eigen::Vector3d pointingSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] The LANDMARKS section: the pixel measurements this image contributes.
    std::vector< SumLandmarkObservation > landmarkObservations_;

    //! [CERTAIN] The LIMB FITS section. Parsed, but not turned into observations by this codebase.
    //! Absent from every SUM file in the archived 67P fixture, so this path is exercised only synthetically.
    std::vector< SumLimbFitObservation > limbFitObservations_;

    //! Not an SPC field: path this record was read from.
    std::string sourceFile_;
};

//! One row of an LMK file's PICTURES section: which images see this landmark, and where.
//! Redundant with the SUM side, which already carries the same measurements, so this codebase parses it
//! but builds its observations from the SUM files. The archived 67P fixture has this section trimmed out.
struct LmkPictureObservation {
    //! [CERTAIN] Name of the image that sees this landmark.
    std::string imageId_;

    //! [CERTAIN] Landmark centre in that image, (sample/column, line/row), in pixels.
    Eigen::Vector2d pixelCoordinates_ = Eigen::Vector2d::Zero( );

    //! [UNKNOWN] A trailing marker present on some rows. Whether it means "rejected", "held fixed", or
    //! something else is not established. Parsed, never acted on.
    bool flagged_ = false;
};

//! Parsed contents of one LMK file: one landmark, i.e. the stereo-solved location and orientation of one
//! "maplet" (the small local terrain patch SPC correlates against each image).
//!
//! [WIKI] "Landmark (.LMK) files describe the locations and orientations of maplets, the pixel/line
//! locations of the maplets in images and on image limbs, and the relative locations of overlapping
//! maplets."
struct LmkLandmarkData {
    //! [CERTAIN] Landmark name. [GIANT] "This should be a 6 character string, typically alpha-numeric."
    //! Becomes the name of the body-fixed ground station this landmark is registered as, and hence the
    //! reference point of the pixel observation's transmitter link end.
    std::string landmarkId_;

    //! [UNKNOWN] Single character on the first line, after the name. 'T' throughout the archived 67P
    //! dataset. Neither [WIKI] nor [GIANT] documents it. Parsed, never interpreted.
    char typeFlag_ = '\0';

    //! [CERTAIN] SIZE: "Half the number of grid cells on a side for the corresponding Maplet minus 1"
    //! [GIANT]. So the maplet is (2 * patchSize_ + 1) cells on a side - 99 for the archived value of 49.
    //! It is a HALF-width, not a width: do not read it as the maplet's extent.
    int patchSize_ = 0;

    //! [CERTAIN] SCALE: "The ground sample distance of each grid cell for the corresponding Maplet in
    //! units of KM" [GIANT]. Left in KILOMETRES, unconverted, because the SPC label carries the unit.
    double patchScaleKm_ = TUDAT_NAN;

    //! [UNKNOWN] HORIZON: four integers, -1 -1 -1 -1 throughout the archived 67P dataset. Neither source
    //! documents them; the name suggests limb/horizon usage flags but that is a guess. Parsed, never
    //! interpreted.
    Eigen::Vector4i horizonFlags_ = Eigen::Vector4i::Zero( );

    //! [UNVERIFIED] SIGKM: a scalar uncertainty associated with the landmark position. [GIANT] sigkm,
    //! "The 1 sigma uncertainty on the landmark body fixed vector in units of km?" - THE QUESTION MARK
    //! IS GIANT'S OWN, so even that reader is unsure. Left in KILOMETRES, unconverted, because the field
    //! name carries the unit.
    //!
    //! Distinct from landmarkPositionSigma_ and NOT interchangeable with it. Across the 40 archived
    //! landmarks SIGKM spans 1.25e-3 to 5.0e-3 km and SIGMA_LMK 2.2e-4 to 7.9e-4 km, a ratio of 3.2 to
    //! 9.3 (median 4.7).
    //!
    //! [MEAS] hints at, but does not settle, which is which: SIGKM takes only THREE distinct values over
    //! all 40 landmarks, which is the signature of a per-landmark-class a-priori the operator chose,
    //! whereas SIGMA_LMK varies continuously as a solved formal uncertainty would. Suggestive, not
    //! conclusive - 40 landmarks from one arc is a small sample. This codebase does NOT use SIGKM; the
    //! landmark a-priori is built from SIGMA_LMK.
    double sigKm_ = TUDAT_NAN;

    //! [UNVERIFIED] RMSLMK: a residual RMS for this landmark's solution. [GIANT] rmslmk, "The residual
    //! RMS of the landmark body fixed vector in units of km?" - again GIANT'S OWN QUESTION MARK, and
    //! here [MEAS] actively contradicts the kilometre reading. Across the 40 archived landmarks the
    //! value spans 0.43 to 2.20 with a mean of 0.77; as kilometres that would be 430 m to 2.2 km of fit
    //! residual on a body only ~2 km across, which is impossible. Clustering around 1 instead points to
    //! a DIMENSIONLESS normalised residual (residual divided by its own sigma), which is about 1 for a
    //! well-fitted landmark. That is inference from 40 samples, not documentation.
    //!
    //! Left unconverted and unused. DO NOT treat this as a length without settling it first.
    double rmsLmk_ = TUDAT_NAN;

    //! [CERTAIN] VLM: the landmark position in the target body-fixed frame. Converted from km to METRES
    //! here. [GIANT] vlm, "The body-fixed landmark vector in units of km". This is the quantity the
    //! landmark position parameter (ground_station_position on the target body) estimates.
    Eigen::Vector3d bodyFixedPosition_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] UX, UY, UZ: the maplet local frame axes, each a unit vector expressed in the TARGET
    //! BODY-FIXED frame. They are therefore the COLUMNS of the rotation from the maplet local frame to
    //! the body-fixed frame, and localZAxis_ is the maplet plane NORMAL.
    //!
    //! [WIKI] MAPFILES states it directly: "Ux body fixed unit map axis vector", "Uy body fixed unit map
    //! axis vector", "Uz body fixed unit map normal vector". [GIANT] names the same quantity
    //! rot_map2bod, "The rotation matrix from the corresponding Maplet local frame to the body-fixed
    //! frame". [MEAS] rules out the alternative reading in which the three records are the ROWS: that
    //! would give an inward-pointing local Z for 6 of the 40 archived landmarks, whereas the columns
    //! reading gives an outward normal for all 40 (cosine with the outward radial +0.51 to +0.95).
    //!
    //! Beware of one thing the normal is NOT: it is not the line of sight. Across the archived
    //! observations it sits a mean 38 degrees away from it, so "the normal direction is the poorly
    //! determined one" is a conflation - the weak direction in an image is the line of sight.
    //!
    //! Orthonormality is checked against lmkLandmarkFrameTolerance, NOT
    //! sumCameraRotationMatrixTolerance: see that constant for why.
    Eigen::Vector3d localXAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d localYAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    Eigen::Vector3d localZAxis_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] SIGMA_LMK: per-component 1-sigma uncertainty on the landmark position, in the
    //! BODY-FIXED frame - not in the UX/UY/UZ local frame. Converted from km to METRES here. [GIANT]
    //! sigma_lmk, "The 1 sigma uncertainty on all of the components of the landmark body fixed vector in
    //! units of km", stated without the question mark that sigkm and rmslmk carry. A single value in the
    //! file is expanded to all three components.
    //!
    //! Being body-fixed is what makes the landmark a-priori a plain diagonal on the body-fixed
    //! ground_station_position parameter. Expressed in the local frame it would instead be the full
    //! block U^T diag(1/sigma^2) U.
    //!
    //! NOT independent information: like SIGMA_VSO, it is an output of SPC's own adjustment of the same
    //! pixels being fitted, so using it as a prior double-counts the data to some degree. It is still
    //! needed, because something has to regularise the landmark network against the orbit.
    Eigen::Vector3d landmarkPositionSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );

    //! [CERTAIN] The PICTURES section. Redundant with the SUM side; see LmkPictureObservation.
    std::vector< LmkPictureObservation > pictures_;

    //! Not an SPC field: path this record was read from.
    std::string sourceFile_;
};

SumImageData readSumFile( const std::string& sumFile );

std::vector< SumImageData > readSumFiles( const std::vector< std::string >& sumFiles );

LmkLandmarkData readLmkFile( const std::string& lmkFile );

std::map< std::string, LmkLandmarkData > readLmkFiles( const std::vector< std::string >& lmkFiles );

//! Name of the camera created for one SUM image: "Camera_" followed by the image ID.
std::string getSumImageCameraName( const std::string& imageId );

//! Observation epoch of a SUM image, in seconds since J2000 in the UTC time scale.
double getSumImageObservationTime( const SumImageData& image );

//! Rotation from the target body-fixed frame to the camera frame, given by the SUM CX/CY/CZ rows.
Eigen::Quaterniond getSumImageRotationFromBodyFixedToCameraFrame( const SumImageData& image );

//! Check SUM/LMK data before conversion and return the images to convert.
/*!
 *  Throws if an image has unusable camera data, or if two image IDs give the same camera name. If
 *  skipObservationsWithMissingLandmarks is true, landmarks without an LMK file are dropped from the images (and images
 *  left without landmarks are dropped). Otherwise a missing landmark is an error. If validateSpacecraftObjectGeometry
 *  is true, every landmark must lie in front of the camera when placed with the SUM SCOBJ vector.
 *  Repeated images and image/landmark measurements are merged; conflicting metadata or pixels are rejected.
 */
std::vector< SumImageData > getSumImagesToConvert( const std::vector< SumImageData >& sumImages,
                                                   const std::map< std::string, LmkLandmarkData >& landmarks,
                                                   const bool skipObservationsWithMissingLandmarks = false,
                                                   const bool validateSpacecraftObjectGeometry = true );

//! Rotation from a landmark's UX/UY/UZ local frame to the target body-fixed frame.
/*!
 *  UX/UY/UZ are the columns of this rotation, and UZ is the maplet surface normal. Throws if the axes are not finite or
 *  do not form a rotation within lmkLandmarkFrameTolerance. The result is snapped to the nearest exact rotation.
 */
Eigen::Matrix3d getLmkLandmarkRotationFromLocalToBodyFixedFrame( const LmkLandmarkData& landmark );

//! Create the supplementary data for SUM/LMK observations.
/*!
 *  Each landmark becomes a reference point on the target body, carrying its SIGMA_LMK position uncertainty and its
 *  UX/UY/UZ local frame when the LMK file gives them. Each image becomes a camera on the receiver body,
 *  pointed relative to the target body-fixed frame and carrying the SUM SIGMA_PTG pointing uncertainty.
 */
std::vector< std::shared_ptr< data::TrackingSupplementaryData > > getSumLmkTrackingSupplementaryData(
        const std::vector< SumImageData >& sumImages,
        const std::map< std::string, LmkLandmarkData >& landmarks,
        const std::string& targetBodyName,
        const std::string& receiverBodyName );

inline data::PlainLinkDefinition getSumLmkLinkEnds( const std::string& targetBodyName,
                                                    const std::string& landmarkId,
                                                    const std::string& receiverBodyName,
                                                    const std::string& cameraName )
{
    return data::PlainLinkDefinition( { std::make_pair( std::make_pair( targetBodyName, landmarkId ), "transmitter" ),
                                        std::make_pair( std::make_pair( receiverBodyName, cameraName ), "receiver" ) } );
}

//! Convert parsed SUM/LMK data to pixel-coordinate tracking data and supplementary data.
/*!
 *  Every (image, landmark) pair becomes one pixel-coordinate observation, sent by the landmark on the target body and
 *  received by the image's camera on the receiver body. Epochs are kept in UTC.
 */
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
convertSumLmkData( const std::vector< SumImageData >& sumImages,
                   const std::map< std::string, LmkLandmarkData >& landmarks,
                   const std::string& targetBodyName,
                   const std::string& receiverBodyName,
                   const bool skipObservationsWithMissingLandmarks = false,
                   const bool validateSpacecraftObjectGeometry = true )
{
    const std::vector< SumImageData > imagesToConvert =
            getSumImagesToConvert( sumImages, landmarks, skipObservationsWithMissingLandmarks, validateSpacecraftObjectGeometry );

    std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataSets;
    for( const SumImageData& image : imagesToConvert )
    {
        const TimeType observationTime = static_cast< TimeType >( getSumImageObservationTime( image ) );
        const std::string cameraName = getSumImageCameraName( image.imageId_ );
        for( const SumLandmarkObservation& observation : image.landmarkObservations_ )
        {
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observable( 2 );
            observable( 0 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 0 ) );
            observable( 1 ) = static_cast< ObservationScalarType >( observation.pixelCoordinates_( 1 ) );

            trackingDataSets.push_back( std::make_shared< data::TrackingData< ObservationScalarType, TimeType > >(
                    "PixelCoordinates",
                    getSumLmkLinkEnds( targetBodyName, observation.landmarkId_, receiverBodyName, cameraName ),
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( 1, observable ),
                    std::vector< TimeType >( 1, observationTime ),
                    "receiver",
                    "UTC" ) );
        }
    }

    return std::make_pair( trackingDataSets,
                           getSumLmkTrackingSupplementaryData( imagesToConvert, landmarks, targetBodyName, receiverBodyName ) );
}

//! Read SUM and LMK files and convert them to pixel-coordinate tracking data and supplementary data.
template< typename ObservationScalarType = double, typename TimeType = Time >
std::pair< std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > >,
           std::vector< std::shared_ptr< data::TrackingSupplementaryData > > >
readSumLmkData( const std::vector< std::string >& sumFiles,
                const std::vector< std::string >& lmkFiles,
                const std::string& targetBodyName,
                const std::string& receiverBodyName,
                const bool skipObservationsWithMissingLandmarks = false,
                const bool validateSpacecraftObjectGeometry = true )
{
    return convertSumLmkData< ObservationScalarType, TimeType >( readSumFiles( sumFiles ),
                                                                 readLmkFiles( lmkFiles ),
                                                                 targetBodyName,
                                                                 receiverBodyName,
                                                                 skipObservationsWithMissingLandmarks,
                                                                 validateSpacecraftObjectGeometry );
}

}  // namespace sum_lmk

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READ_SUM_LMK_FILES_H
