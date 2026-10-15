#ifndef TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H
#define TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H

#include <optional>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace reference_frames
{

//! Transformation between terrestrial reference frames at one epoch: position_out = translation + matrix * position_in.
//! The matrix combines the (linearized) rotation and the scale, so it is not exactly orthogonal.
struct TerrestrialFrameAffineTransformation {
    //! Translation (metres) and combined rotation and scale matrix (dimensionless).
    Eigen::Vector3d translation = Eigen::Vector3d::Zero( );
    Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity( );

    //! Time derivatives of the translation (metres per second) and of the matrix (per second).
    Eigen::Vector3d translationRate = Eigen::Vector3d::Zero( );
    Eigen::Matrix3d matrixRate = Eigen::Matrix3d::Zero( );

    //! Apply the translation and matrix to a Cartesian position in metres.
    Eigen::Vector3d transformPosition( const Eigen::Vector3d& position ) const;

    //! Include parameter derivatives when transforming position and velocity.
    Eigen::Vector6d transformState( const Eigen::Vector6d& state ) const;

    //! Invert the full matrix and differentiate that inverse.
    TerrestrialFrameAffineTransformation inverse( ) const;

    //! Apply the supplied transformation first, including both matrix and translation derivatives.
    TerrestrialFrameAffineTransformation composeAfter( const TerrestrialFrameAffineTransformation& before ) const;
};

//! Evaluate published Helmert parameters with an explicit year and years-per-input-second rate.
TerrestrialFrameAffineTransformation evaluateLinearizedHelmert( const input_output::TerrestrialFrameTransformationParameters& parameters,
                                                                double decimalYear,
                                                                double decimalYearRateWithRespectToInputSeconds );

//! Linearized Helmert transformations with parameter rates from a preloaded terrestrial catalogue.
class HelmertTransformationModel
{
public:
    //! Validate and retain the original catalogue and its declared rotation convention.
    explicit HelmertTransformationModel( const input_output::TerrestrialFrameTransformationData& data );

    //! Return translation, matrix and their derivatives; equal nonempty realization labels give identity.
    //! The rate is years per input second. Distinct realizations require known hub records.
    TerrestrialFrameAffineTransformation getTransformation( const std::string& sourceRealization,
                                                            const std::string& targetRealization,
                                                            double decimalYear,
                                                            double decimalYearRateWithRespectToInputSeconds ) const;

    //! Transform position at the supplied decimal year; no time derivative is needed.
    Eigen::Vector3d transformPosition( const Eigen::Vector3d& position,
                                       const std::string& sourceRealization,
                                       const std::string& targetRealization,
                                       double decimalYear ) const;

    //! Transform position and velocity using an explicit year rate in the state's input seconds.
    Eigen::Vector6d transformState( const Eigen::Vector6d& state,
                                    const std::string& sourceRealization,
                                    const std::string& targetRealization,
                                    double decimalYear,
                                    double decimalYearRateWithRespectToInputSeconds ) const;

    //! Expose the retained neutral catalogue and its original provenance.
    const input_output::TerrestrialFrameTransformationData& getData( ) const
    {
        return data_;
    }

private:
    //! Share source-to-hub-to-target evaluation; absent rates are used only for position.
    TerrestrialFrameAffineTransformation evaluateTransformation( const std::string& sourceRealization,
                                                                 const std::string& targetRealization,
                                                                 double decimalYear,
                                                                 std::optional< double > yearRate ) const;

    //! Published catalogue identities and SI-normalized parameters.
    input_output::TerrestrialFrameTransformationData data_;
};

}  // namespace reference_frames
}  // namespace tudat
#endif
