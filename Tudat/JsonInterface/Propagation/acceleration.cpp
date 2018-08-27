/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Propagation/acceleration.h"
#include "Tudat/JsonInterface/Propagation/thrust.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    if ( ! accelerationSettings )
    {
        return;
    }
    using namespace basic_astrodynamics;
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration;

    const AvailableAcceleration accelerationType = accelerationSettings->accelerationType_;
    jsonObject[ K::type ] = accelerationType;

    switch ( accelerationType ) {
    case undefined_acceleration:
    case aerodynamic:
    case cannon_ball_radiation_pressure:
    {
        return;
    }
    case point_mass_gravity:
    case third_body_point_mass_gravity:
    {
        jsonObject[ K::type ] = point_mass_gravity;
        return;
    }
    case spherical_harmonic_gravity:
    case third_body_spherical_harmonic_gravity:
    {
        jsonObject[ K::type ] = spherical_harmonic_gravity;
        std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicAccelerationSettings
                = std::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
        assertNonnullptrPointer( sphericalHarmonicAccelerationSettings );
        jsonObject[ K::maximumDegree ] = sphericalHarmonicAccelerationSettings->maximumDegree_;
        jsonObject[ K::maximumOrder ] = sphericalHarmonicAccelerationSettings->maximumOrder_;
        return;
    }
    case mutual_spherical_harmonic_gravity:
    case third_body_mutual_spherical_harmonic_gravity:
    {
        jsonObject[ K::type ] = mutual_spherical_harmonic_gravity;
        std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicAccelerationSettings
                = std::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >(
                    accelerationSettings );
        assertNonnullptrPointer( mutualSphericalHarmonicAccelerationSettings );
        jsonObject[ K::maximumDegreeOfBodyExertingAcceleration ] =
                mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyExertingAcceleration_;
        jsonObject[ K::maximumOrderOfBodyExertingAcceleration ] =
                mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyExertingAcceleration_;
        jsonObject[ K::maximumDegreeOfBodyUndergoingAcceleration ] =
                mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyUndergoingAcceleration_;
        jsonObject[ K::maximumOrderOfBodyUndergoingAcceleration ] =
                mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyUndergoingAcceleration_;
        jsonObject[ K::maximumDegreeOfCentralBody ] =
                mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfCentralBody_;
        jsonObject[ K::maximumOrderOfCentralBody ] =
                mutualSphericalHarmonicAccelerationSettings->maximumOrderOfCentralBody_;
        return;
    }
    case thrust_acceleration:
    {
        std::shared_ptr< ThrustAccelerationSettings > thrustAccelerationSettings
                = std::dynamic_pointer_cast< ThrustAccelerationSettings >( accelerationSettings );
        assertNonnullptrPointer( thrustAccelerationSettings );
        jsonObject = thrustAccelerationSettings;
        return;
    }
    case relativistic_correction_acceleration:
    {
        std::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationCorrectionSettings
                = std::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >( accelerationSettings );
        assertNonnullptrPointer( relativisticAccelerationCorrectionSettings );
        jsonObject[ K::calculateSchwarzschildCorrection ] =
                relativisticAccelerationCorrectionSettings->calculateSchwarzschildCorrection_;
        jsonObject[ K::calculateLenseThirringCorrection ] =
                relativisticAccelerationCorrectionSettings->calculateLenseThirringCorrection_;
        jsonObject[ K::calculateDeSitterCorrection ] =
                relativisticAccelerationCorrectionSettings->calculateDeSitterCorrection_;
        jsonObject[ K::primaryBody ] =
                relativisticAccelerationCorrectionSettings->primaryBody_;
        jsonObject[ K::centralBodyAngularMomentum ] =
                relativisticAccelerationCorrectionSettings->centralBodyAngularMomentum_;
        return;
    }
    case empirical_acceleration:
    {
        std::shared_ptr< EmpiricalAccelerationSettings > empiricalAccelerationSettings
                = std::dynamic_pointer_cast< EmpiricalAccelerationSettings >( accelerationSettings );
        assertNonnullptrPointer( empiricalAccelerationSettings );
        jsonObject[ K::constantAcceleration ] = empiricalAccelerationSettings->constantAcceleration_;
        jsonObject[ K::sineAcceleration ] = empiricalAccelerationSettings->sineAcceleration_;
        jsonObject[ K::cosineAcceleration ] = empiricalAccelerationSettings->cosineAcceleration_;
        return;
    }
    default:
        handleUnimplementedEnumValue( accelerationType, accelerationTypes, unsupportedAccelerationTypes );
    }
}

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    using namespace json_interface;
    using namespace basic_astrodynamics;
    using K = Keys::Propagator::Acceleration;

    // Get acceleration type
    const AvailableAcceleration accelerationType = getValue< AvailableAcceleration >( jsonObject, K::type );

    switch ( accelerationType ) {
    case undefined_acceleration:
    case point_mass_gravity:
    case aerodynamic:
    case cannon_ball_radiation_pressure:
    {
        accelerationSettings = std::make_shared< AccelerationSettings >( accelerationType );
        return;
    }
    case spherical_harmonic_gravity:
    {
        accelerationSettings = std::make_shared< SphericalHarmonicAccelerationSettings >(
                    getValue< int >( jsonObject, K::maximumDegree ),
                    getValue< int >( jsonObject, K::maximumOrder ) );
        return;
    }
    case mutual_spherical_harmonic_gravity:
    {
        MutualSphericalHarmonicAccelerationSettings defaults( -1, -1, -1, -1 );
        accelerationSettings = std::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                    getValue< int >( jsonObject, K::maximumDegreeOfBodyExertingAcceleration ),
                    getValue< int >( jsonObject, K::maximumOrderOfBodyExertingAcceleration ),
                    getValue< int >( jsonObject, K::maximumDegreeOfBodyUndergoingAcceleration ),
                    getValue< int >( jsonObject, K::maximumOrderOfBodyUndergoingAcceleration ),
                    getValue( jsonObject, K::maximumDegreeOfCentralBody, defaults.maximumDegreeOfCentralBody_ ),
                    getValue( jsonObject, K::maximumOrderOfCentralBody, defaults.maximumOrderOfCentralBody_ ) );
        return;
    }
    case thrust_acceleration:
    {
        accelerationSettings = getAs< std::shared_ptr< ThrustAccelerationSettings > >( jsonObject );
        return;
    }
    case relativistic_correction_acceleration:
    {
        RelativisticAccelerationCorrectionSettings defaults;
        accelerationSettings = std::make_shared< RelativisticAccelerationCorrectionSettings >(
                    getValue( jsonObject, K::calculateSchwarzschildCorrection,
                              defaults.calculateSchwarzschildCorrection_ ),
                    getValue( jsonObject, K::calculateLenseThirringCorrection,
                              defaults.calculateLenseThirringCorrection_ ),
                    getValue( jsonObject, K::calculateDeSitterCorrection, defaults.calculateDeSitterCorrection_ ),
                    getValue( jsonObject, K::primaryBody, defaults.primaryBody_ ),
                    getValue( jsonObject, K::centralBodyAngularMomentum, defaults.centralBodyAngularMomentum_ ) );
        return;
    }
    case empirical_acceleration:
    {
        EmpiricalAccelerationSettings defaults;
        accelerationSettings = std::make_shared< EmpiricalAccelerationSettings >(
                    getValue( jsonObject, K::constantAcceleration, defaults.constantAcceleration_ ),
                    getValue( jsonObject, K::sineAcceleration, defaults.sineAcceleration_ ),
                    getValue( jsonObject, K::cosineAcceleration, defaults.cosineAcceleration_ ) );
        return;
    }
    default:
    {
        if ( accelerationType == third_body_point_mass_gravity ||
             accelerationType == third_body_spherical_harmonic_gravity ||
             accelerationType == third_body_mutual_spherical_harmonic_gravity )
        {
            std::cerr << "Whether a body will cause a third-body acceleration is determined internally "
                      << "by Tudat based on the propagation settings." << std::endl;
        }
        handleUnimplementedEnumValue( accelerationType, accelerationTypes, unsupportedAccelerationTypes );
    }
    }
}

} // namespace simulation_setup

} // namespace tudat

namespace boost
{

template<>
tudat::basic_astrodynamics::EmpiricalAccelerationComponents lexical_cast(const std::string & s) {
    return tudat::basic_astrodynamics::empiricalAccelerationComponentTypesInverse.at( s );
}

template<>
std::string lexical_cast(const tudat::basic_astrodynamics::EmpiricalAccelerationComponents & component )
{
    return tudat::basic_astrodynamics::empiricalAccelerationComponentTypes.at( component );
}

}
