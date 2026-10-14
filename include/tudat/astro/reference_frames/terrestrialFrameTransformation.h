#ifndef TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H
#define TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H

#include <optional>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace reference_frames
{

//! Linearized Helmert map; the matrix includes scale and need not be orthogonal.
struct TerrestrialFrameAffineTransformation {
    //! Metre-valued translation and dimensionless rotation/scale matrix.
    Eigen::Vector3d translation = Eigen::Vector3d::Zero( );

    Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity( );

    //! Derivatives in metres per input second and inverse input seconds.
    Eigen::Vector3d translationRate = Eigen::Vector3d::Zero( );

    Eigen::Matrix3d matrixRate = Eigen::Matrix3d::Zero( );

    //! Apply the affine map to a Cartesian position in metres.
    Eigen::Vector3d transformPosition( const Eigen::Vector3d& position ) const;

    //! Include parameter derivatives when transforming position and velocity.
    Eigen::Vector6d transformState( const Eigen::Vector6d& state ) const;

    //! Invert the full matrix and differentiate that inverse.
    TerrestrialFrameAffineTransformation inverse( ) const;

    //! Apply the supplied map first, including both sets of affine derivatives.
    TerrestrialFrameAffineTransformation composeAfter( const TerrestrialFrameAffineTransformation& before ) const;

};

//! Evaluate published Helmert parameters with an explicit year and years-per-input-second rate.
TerrestrialFrameAffineTransformation evaluateLinearizedHelmert(
        const input_output::TerrestrialFrameTransformationParameters& parameters,
        double decimalYear,
        double decimalYearRateWithRespectToInputSeconds );

//! Apply a preloaded terrestrial catalogue with the caller's explicit time convention.
class TerrestrialFrameTransformationModel
{
public:
    //! Validate and retain the original catalogue and its declared rotation convention.
    explicit TerrestrialFrameTransformationModel( const input_output::TerrestrialFrameTransformationData& data );

    //! Return affine values and derivatives; equal nonempty realization labels give identity.
    //! The rate is years per input second. Distinct realizations require known hub records.
    TerrestrialFrameAffineTransformation getTransformation( const std::string& sourceRealization,
                                                            const std::string& targetRealization,
                                                            double decimalYear,
                                                            double decimalYearRateWithRespectToInputSeconds ) const;

    //! Transform position from a declared decimal year without requiring or inventing a rate.
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
