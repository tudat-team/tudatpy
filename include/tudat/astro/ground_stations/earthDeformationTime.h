#ifndef TUDAT_EARTH_DEFORMATION_TIME_H
#define TUDAT_EARTH_DEFORMATION_TIME_H

#include <memory>
#include <optional>
#include "tudat/basics/timeType.h"

namespace tudat
{
namespace earth_orientation
{
class TerrestrialTimeScaleConverter;
}
namespace basic_astrodynamics
{
//! One physical epoch with distinct environment TDB and terrestrial clocks.
/*!
 * The station environment supplies TDB seconds since J2000. Conversion is lazy:
 * generic/custom callbacks can consume the unchanged TDB scalar without invoking
 * a terrestrial clock. Native Earth models use TT or continuous UT1 explicitly.
 * The geocentric TT/TAI clock needs no EOP data; UT1 requires a supplied converter.
 * Original double inputs retain their type for canonical endpoint-image handling;
 * original Time inputs remain precise and are never snapped to an endpoint.
 */
class EarthDeformationEpoch
{
public:
    using Converter = earth_orientation::TerrestrialTimeScaleConverter;
    EarthDeformationEpoch( const Time& tdb, const std::shared_ptr< Converter >& converter = nullptr ):
        tdb_( tdb ), originalDouble_( false ), converter_( converter )
    {}
    EarthDeformationEpoch( double tdb, const std::shared_ptr< Converter >& converter = nullptr ):
        tdb_( tdb ), originalDouble_( true ), converter_( converter )
    {}
    const Time& getTdb( ) const
    {
        return tdb_;
    }
    double getTdbSeconds( ) const
    {
        return static_cast< double >( tdb_ );
    }
    const Time& getTai( ) const;
    const Time& getTt( ) const;
    Time getUt1( ) const;
    //! Evaluate UT1 with the requesting model's own converter and its EOP policy.
    Time getUt1( const std::shared_ptr< Converter >& converter ) const;
    //! Derivative of the geocentric terrestrial clock with respect to TDB.
    double getTtRateWithRespectToTdb( ) const;
    bool originatedAsDouble( ) const
    {
        return originalDouble_;
    }
    EarthDeformationEpoch withTdb( const Time& tdb ) const
    {
        return EarthDeformationEpoch( tdb, converter_ );
    }
    //! Replace the clock while preserving the original epoch and its precision.
    EarthDeformationEpoch withTimeConverter( const std::shared_ptr< Converter >& converter ) const
    {
        return originalDouble_ ? EarthDeformationEpoch( getTdbSeconds( ), converter ) : EarthDeformationEpoch( tdb_, converter );
    }

private:
    const std::shared_ptr< Converter >& clock( ) const;
    Time tdb_;
    bool originalDouble_;
    std::shared_ptr< Converter > converter_;
    mutable std::optional< Time > tai_;
    mutable std::optional< Time > tt_;
};

//! Validate a closed TT interval for the current Earth models' scalar-UTC phase convention.
/*!
 * Reject ambiguous scalar UTC samples and UTC offset steps, including steps at either
 * endpoint. Smooth pre-1972 UTC drift is permitted. This is a limitation of the current
 * UTC phase representation, not a physical tide singularity. Custom argument callbacks
 * with a different time basis are responsible for their own interval policy.
 */
void validateEarthDeformationUtcInterval( const double lowerTt, const double upperTt );

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif
