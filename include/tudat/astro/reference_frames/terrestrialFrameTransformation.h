#ifndef TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H
#define TUDAT_TERRESTRIAL_FRAME_TRANSFORMATION_H
#include <functional>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/basicTypedefs.h"
namespace tudat
{
namespace reference_frames
{
//! Published parameter coordinate and its derivative relative to the input state's seconds.
//! Unit-normalizing a per-year rate does not define this physical time mapping.
struct TerrestrialFrameEvaluationEpoch {
    double decimalYear;
    double decimalYearRateWithRespectToInputSeconds;
};
//! Linearized Helmert affine transformation; its matrix includes scale and is not an orthogonal rotation.
struct TerrestrialFrameAffineTransformation {
    Eigen::Vector3d translation = Eigen::Vector3d::Zero( );
    Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity( );
    Eigen::Vector3d translationRate = Eigen::Vector3d::Zero( );
    Eigen::Matrix3d matrixRate = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d transformPosition( const Eigen::Vector3d& position ) const;
    Eigen::Vector6d transformState( const Eigen::Vector6d& state ) const;
    TerrestrialFrameAffineTransformation inverse( ) const;
    //! Apply this transformation after the supplied transformation, including all affine derivative terms.
    TerrestrialFrameAffineTransformation composeAfter( const TerrestrialFrameAffineTransformation& before ) const;
};
TerrestrialFrameAffineTransformation evaluateLinearizedHelmert( const input_output::TerrestrialFrameTransformationParameters& parameters,
                                                                const TerrestrialFrameEvaluationEpoch& epoch );
enum class ItrfTransformationCatalogue { ignItrf2014, ignItrf2020 };
//! Select pinned conventional input through its neutral reader, without redefining the table's identity.
input_output::TerrestrialFrameTransformationData getConventionalItrfTransformationData( ItrfTransformationCatalogue catalogue );
//! Native model owns an explicitly supplied physical-time mapping; evaluation never opens files.
class TerrestrialFrameTransformationModel
{
public:
    using Clock = std::function< TerrestrialFrameEvaluationEpoch( const Time& ) >;
    explicit TerrestrialFrameTransformationModel( const input_output::TerrestrialFrameTransformationData& data, Clock clock = {} );
    //! Equal nonempty realization labels give identity before catalogue lookup, including unknown labels.
    //! This identity behavior does not validate catalogue membership; distinct labels require known hub records.
    TerrestrialFrameAffineTransformation getTransformation( const std::string& sourceRealization,
                                                            const std::string& targetRealization,
                                                            const TerrestrialFrameEvaluationEpoch& epoch ) const;
    TerrestrialFrameAffineTransformation getTransformationAtTime( const std::string& sourceRealization,
                                                                  const std::string& targetRealization,
                                                                  const Time& inputTime ) const;
    Eigen::Vector3d transformPosition( const Eigen::Vector3d& position,
                                       const std::string& sourceRealization,
                                       const std::string& targetRealization,
                                       const TerrestrialFrameEvaluationEpoch& epoch ) const;
    //! Value-only calendar/declared coordinate evaluation; no physical rate is needed or invented.
    Eigen::Vector3d transformPositionAtDecimalYear( const Eigen::Vector3d& position,
                                                    const std::string& sourceRealization,
                                                    const std::string& targetRealization,
                                                    double decimalYear ) const;
    Eigen::Vector6d transformState( const Eigen::Vector6d& state,
                                    const std::string& sourceRealization,
                                    const std::string& targetRealization,
                                    const TerrestrialFrameEvaluationEpoch& epoch ) const;
    Eigen::Vector3d transformPositionAtTime( const Eigen::Vector3d& position,
                                             const std::string& sourceRealization,
                                             const std::string& targetRealization,
                                             const Time& inputTime ) const;
    Eigen::Vector6d transformStateAtTime( const Eigen::Vector6d& state,
                                          const std::string& sourceRealization,
                                          const std::string& targetRealization,
                                          const Time& inputTime ) const;
    const input_output::TerrestrialFrameTransformationData& getData( ) const
    {
        return data_;
    }

private:
    input_output::TerrestrialFrameTransformationData data_;
    Clock clock_;
};
//! Compatibility coordinate: uniform Julian-year seconds since J2000; not calendar decimal year or TDB conversion.
TerrestrialFrameEvaluationEpoch getLegacyItrfEvaluationEpoch( double secondsSinceJ2000 );
}  // namespace reference_frames
}  // namespace tudat
#endif
