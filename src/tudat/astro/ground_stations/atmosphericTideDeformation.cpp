#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace basic_astrodynamics
{
namespace
{
template< typename Convention >
Convention resolve( Convention source, Convention declared, Convention unspecified, const char* name )
{
    if( source != unspecified && declared != unspecified && source != declared )
    {
        throw std::runtime_error( std::string( "Conflicting atmospheric tide " ) + name + " declarations" );
    }
    return declared == unspecified ? source : declared;
}
input_output::AtmosphericTideMetadata stationMetadata( const input_output::AtmosphericTideMetadata& station,
                                                       const input_output::AtmosphericTideMetadata& common )
{
    using namespace input_output;
    auto result = station;
    result.displacementBasis =
            resolve( station.displacementBasis, common.displacementBasis, AtmosphericTideDisplacementBasis::unspecified, "basis" );
    result.referenceOrigin =
            resolve( station.referenceOrigin, common.referenceOrigin, AtmosphericLoadingReferenceOrigin::unspecified, "origin" );
    result.sourceProfile = resolve( station.sourceProfile, common.sourceProfile, AtmosphericTideSourceProfile::unspecified, "profile" );
    if( !station.modelName.empty( ) && !common.modelName.empty( ) && station.modelName != common.modelName )
    {
        throw std::runtime_error( "Conflicting atmospheric tide model declarations" );
    }
    if( result.modelName.empty( ) )
    {
        result.modelName = common.modelName;
    }
    if( !station.greensFunction.empty( ) && !common.greensFunction.empty( ) && station.greensFunction != common.greensFunction )
    {
        throw std::runtime_error( "Conflicting atmospheric Green-function declarations" );
    }
    if( result.greensFunction.empty( ) )
    {
        result.greensFunction = common.greensFunction;
    }
    return result;
}
}  // namespace

void AtmosphericTideDeformation::validateCommonSettings( ) const
{
    using namespace input_output;
    if( settings_.lookupLatitudeConvention != AtmosphericTideLatitudeConvention::unspecified &&
        settings_.lookupLatitudeConvention != AtmosphericTideLatitudeConvention::geocentric &&
        settings_.lookupLatitudeConvention != AtmosphericTideLatitudeConvention::geodetic )
    {
        throw std::runtime_error( "Unknown atmospheric lookup latitude convention" );
    }
    if( settings_.displacementBasis != AtmosphericTideDisplacementBasis::unspecified &&
        settings_.displacementBasis != AtmosphericTideDisplacementBasis::geocentricRadialNorthEast &&
        settings_.displacementBasis != AtmosphericTideDisplacementBasis::geodeticUpNorthEast )
    {
        throw std::runtime_error( "Unknown atmospheric displacement basis" );
    }
}
void AtmosphericTideDeformation::validateShapeForBasis( input_output::AtmosphericTideDisplacementBasis basis, bool geodeticLookup ) const
{
    if( basis == input_output::AtmosphericTideDisplacementBasis::unspecified )
    {
        throw std::runtime_error( "Atmospheric loading requires an explicit displacement basis" );
    }
    if( basis == input_output::AtmosphericTideDisplacementBasis::geodeticUpNorthEast || geodeticLookup )
    {
        if( !bodyShape_ )
        {
            throw std::runtime_error( "Geodetic atmospheric lookup/basis requires a configured body shape" );
        }
        if( !std::dynamic_pointer_cast< OblateSpheroidBodyShapeModel >( bodyShape_ ) &&
            !std::dynamic_pointer_cast< SphericalBodyShapeModel >( bodyShape_ ) )
        {
            throw std::runtime_error( "Unsupported configured shape for geodetic atmospheric conversion" );
        }
    }
}
void AtmosphericTideDeformation::validateOriginConversion( const input_output::AtmosphericTideMetadata& metadata ) const
{
    using namespace input_output;
    if( !settings_.convertToCentreOfMass )
    {
        return;
    }
    if( metadata.sourceProfile != AtmosphericTideSourceProfile::iers2010Rp03NoIbFarrell || metadata.modelName != "Ray-Ponte 2003" ||
        metadata.greensFunction != "Farrell elastic Green functions" )
    {
        throw std::runtime_error( "IERS atmospheric origin translation requires the explicitly identified RP03 no-IB Farrell product" );
    }
    if( metadata.referenceOrigin != AtmosphericLoadingReferenceOrigin::solidEarthCentre &&
        metadata.referenceOrigin != AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass )
    {
        throw std::runtime_error( "IERS atmospheric origin translation requires a declared CE or CM source origin" );
    }
}

AtmosphericTideDeformation::AtmosphericTideDeformation(
        const input_output::AtmosphericTideGrid& grid,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
        const AtmosphericTideEvaluationSettings& settings,
        const std::shared_ptr< BodyShapeModel >& bodyShape ):
    BodyDeformationModel( converter ), settings_( settings ), bodyShape_( bodyShape )
{
    using namespace input_output;
    if( !converter ) throw std::runtime_error( "Atmospheric tide deformation requires its own shared Earth time converter" );
    validateCommonSettings( );
    validateAtmosphericTideGrid( grid );
    resolvedGridMetadata_ = grid.metadata;
    resolvedGridMetadata_.latitudeConvention = resolve( grid.metadata.latitudeConvention,
                                                        settings.lookupLatitudeConvention,
                                                        AtmosphericTideLatitudeConvention::unspecified,
                                                        "lookup latitude" );
    resolvedGridMetadata_.displacementBasis =
            resolve( grid.metadata.displacementBasis, settings.displacementBasis, AtmosphericTideDisplacementBasis::unspecified, "basis" );

    // RP03 grids: default to the geocentric convention when neither the source nor the settings declare one.
    if( grid.metadata.modelName == "Ray-Ponte 2003" && grid.metadata.greensFunction == "Farrell elastic Green functions" )
    {
        if( resolvedGridMetadata_.latitudeConvention == AtmosphericTideLatitudeConvention::unspecified )
        {
            resolvedGridMetadata_.latitudeConvention = AtmosphericTideLatitudeConvention::geocentric;
        }
        if( resolvedGridMetadata_.displacementBasis == AtmosphericTideDisplacementBasis::unspecified )
        {
            resolvedGridMetadata_.displacementBasis = AtmosphericTideDisplacementBasis::geocentricRadialNorthEast;
        }
    }
    if( resolvedGridMetadata_.latitudeConvention == AtmosphericTideLatitudeConvention::unspecified )
    {
        throw std::runtime_error( "Automatic atmospheric grid lookup requires an explicit latitude convention" );
    }
    validateShapeForBasis( resolvedGridMetadata_.displacementBasis,
                           resolvedGridMetadata_.latitudeConvention == AtmosphericTideLatitudeConvention::geodetic );
    validateOriginConversion( resolvedGridMetadata_ );
    gridInterpolator_ = std::make_shared< AtmosphericTideGridInterpolator >( grid );
}
AtmosphericTideDeformation::AtmosphericTideDeformation(
        const input_output::AtmosphericTideStationData& stations,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter,
        const AtmosphericTideEvaluationSettings& settings,
        const std::shared_ptr< BodyShapeModel >& bodyShape ):
    BodyDeformationModel( converter ), stationData_( stations ), settings_( settings ), bodyShape_( bodyShape )
{
    using namespace input_output;
    if( !converter ) throw std::runtime_error( "Atmospheric tide deformation requires its own shared Earth time converter" );
    validateCommonSettings( );
    validateAtmosphericTideStationData( stations );
    for( const auto& entry : stations.stations )
    {
        auto metadata = stationMetadata( entry.second.metadata, stations.metadata );
        metadata.displacementBasis =
                resolve( metadata.displacementBasis, settings.displacementBasis, AtmosphericTideDisplacementBasis::unspecified, "basis" );
        validateShapeForBasis( metadata.displacementBasis );
        validateOriginConversion( metadata );
        resolvedStationMetadata_.emplace( entry.first, metadata );
    }
}

//! Evaluate displacement alone so valid positions do not require an available UT1 rate.
Eigen::Vector3d AtmosphericTideDeformation::evaluateDisplacement( const Time& tdb,
                                                                  const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluateCartesianLoading( tdb, station, false ).head< 3 >( );
}

//! Evaluate Cartesian position and analytic rate with one coefficient, geometry and clock calculation.
Eigen::Vector6d AtmosphericTideDeformation::evaluateDisplacementState(
        const Time& tdb,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluateCartesianLoading( tdb, station, true );
}

//! Synthesize metre displacements and optional UT1 derivatives, then rotate and apply the requested origin correction.
Eigen::Vector6d AtmosphericTideDeformation::evaluateCartesianLoading( const Time& tdb,
                                                                      const std::shared_ptr< ground_stations::GroundStationState >& station,
                                                                      bool withRate )
{
    using namespace input_output;
    if( !station )
    {
        throw std::runtime_error( "Missing atmospheric loading station" );
    }

    // Nominal body-fixed geometry and angular coordinates in radians.
    const Eigen::Vector3d position = station->getNominalCartesianPosition( );

    if( !position.allFinite( ) || position.norm( ) == 0 )
    {
        throw std::runtime_error( "Invalid atmospheric station position" );
    }
    const double longitude = std::atan2( position.y( ), position.x( ) );

    const double geocentricLatitude = std::atan2( position.z( ), std::hypot( position.x( ), position.y( ) ) );

    // Provider declarations and metre coefficients, resolved for the selected station or grid.
    AtmosphericTideMetadata metadata;

    Eigen::Matrix< double, 3, 4 > coefficients;

    if( usesStationCoefficients( ) )
    {
        const auto& record = stationData_.getStation( station->getSiteId( ) );

        metadata = resolvedStationMetadata_.at( station->getSiteId( ) );
        coefficients = record.coefficientsMetres;
    }
    else
    {
        metadata = resolvedGridMetadata_;
    }
    double geodeticLatitude = geocentricLatitude;

    if( metadata.displacementBasis == AtmosphericTideDisplacementBasis::geodeticUpNorthEast ||
        ( !usesStationCoefficients( ) && metadata.latitudeConvention == AtmosphericTideLatitudeConvention::geodetic ) )
    {
        // Use the configured body shape, independently of how this Cartesian point was represented by the station.
        geodeticLatitude =
                coordinate_conversions::convertPositionElements(
                        position, coordinate_conversions::cartesian_position, coordinate_conversions::geodetic_position, bodyShape_ )
                        .y( );
    }
    if( !usesStationCoefficients( ) )
    {
        const double latitude =
                metadata.latitudeConvention == AtmosphericTideLatitudeConvention::geodetic ? geodeticLatitude : geocentricLatitude;

        coefficients =
                gridInterpolator_->interpolate( latitude * 180 / mathematical_constants::PI, longitude * 180 / mathematical_constants::PI );
    }

    // Obtain one UT1 phase; only a state request asks for its physical derivative.
    const Time tai = timeConverter_->getTaiTime( tdb_scale, tdb );

    double ut1, ut1RateWrtTdb = 1;

    if( withRate )
    {
        const auto ut1WithRate = timeConverter_->getUt1FromTaiWithRate( tai );

        ut1 = static_cast< double >( ut1WithRate.ut1 );
        ut1RateWrtTdb = ut1WithRate.dUt1DTai * timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
    }
    else
    {
        ut1 = static_cast< double >( timeConverter_->getUt1FromTai( tai ) );
    }

    // The same nominal component basis applies to the displacement and its derivative.
    const double basisLatitude =
            metadata.displacementBasis == AtmosphericTideDisplacementBasis::geodeticUpNorthEast ? geodeticLatitude : geocentricLatitude;

    const auto rotation =
            reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, basisLatitude );

    const Eigen::Vector3d local = synthesizeAtmosphericTideLoading( coefficients, ut1 );

    // Source R/N/E (or U/N/E) components become the existing helper's E/N/U axes.
    const Eigen::Vector3d enu( local.z( ), local.y( ), local.x( ) );

    Eigen::Vector6d state = Eigen::Vector6d::Zero( );

    state.head< 3 >( ) = rotation * enu;
    if( withRate )
    {
        // Local analytic velocity in metres per UT1 second, before the clock derivative is applied.
        const Eigen::Vector3d localRate = differentiateAtmosphericTideLoadingWrtUt1( coefficients, ut1 );

        state.tail< 3 >( ) = rotation * Eigen::Vector3d( localRate.z( ), localRate.y( ), localRate.x( ) );
    }
    if( settings_.convertToCentreOfMass && metadata.referenceOrigin == AtmosphericLoadingReferenceOrigin::solidEarthCentre )
    {
        state.head< 3 >( ) += calculateIers2010Rp03AtmosphericGeocentreTranslation( ut1 );
        if( withRate )
        {
            state.tail< 3 >( ) += differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( ut1 );
        }
    }
    if( withRate )
    {
        // Convert the full local-plus-origin velocity to metres per TDB second exactly once.
        state.tail< 3 >( ) *= ut1RateWrtTdb;
    }
    return state;
}

input_output::AtmosphericLoadingReferenceOrigin AtmosphericTideDeformation::getOutputReferenceOrigin( const std::string& identifier ) const
{
    if( usesStationCoefficients( ) )
    {
        stationData_.getStation( identifier );
    }
    if( settings_.convertToCentreOfMass )
    {
        return input_output::AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass;
    }
    return usesStationCoefficients( ) ? resolvedStationMetadata_.at( identifier ).referenceOrigin : resolvedGridMetadata_.referenceOrigin;
}
}  // namespace basic_astrodynamics
}  // namespace tudat
