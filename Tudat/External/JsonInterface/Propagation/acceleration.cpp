/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "acceleration.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Convert `AvailableAcceleration` to `json`.
void to_json( json& jsonObject, const AvailableAcceleration& availableAcceleration )
{
    jsonObject = json_interface::stringFromEnum( availableAcceleration, availableAccelerations );
}

//! Convert `json` to `AvailableAcceleration`.
void from_json( const json& jsonObject, AvailableAcceleration& availableAcceleration )
{
    availableAcceleration = json_interface::enumFromString( jsonObject.get< std::string >( ), availableAccelerations );
}

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    if ( accelerationSettings )
    {
        using namespace basic_astrodynamics;
        using namespace json_interface;
        using Keys = Keys::Acceleration;

        const AvailableAcceleration accelerationType = accelerationSettings->accelerationType_;

        switch ( accelerationType ) {
        case undefined_acceleration:
        case aerodynamic:
        case cannon_ball_radiation_pressure:
        {
            jsonObject[ Keys::type ] = accelerationType;
            return;
        }
        case point_mass_gravity:
        case third_body_point_mass_gravity:
        {
            jsonObject[ Keys::type ] = point_mass_gravity;
            return;
        }
        case spherical_harmonic_gravity:
        case third_body_spherical_harmonic_gravity:
        {
            jsonObject[ Keys::type ] = spherical_harmonic_gravity;
            boost::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicAccelerationSettings
                    = boost::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
            jsonObject[ Keys::maximumDegree ] = sphericalHarmonicAccelerationSettings->maximumDegree_;
            jsonObject[ Keys::maximumOrder ] = sphericalHarmonicAccelerationSettings->maximumOrder_;
            return;
        }
        case mutual_spherical_harmonic_gravity:
        case third_body_mutual_spherical_harmonic_gravity:
        {
            jsonObject[ Keys::type ] = mutual_spherical_harmonic_gravity;
            boost::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicAccelerationSettings
                    = boost::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >(
                        accelerationSettings );
            jsonObject[ Keys::maximumDegreeOfBodyExertingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyExertingAcceleration_;
            jsonObject[ Keys::maximumOrderOfBodyExertingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyExertingAcceleration_;
            jsonObject[ Keys::maximumDegreeOfBodyUndergoingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyUndergoingAcceleration_;
            jsonObject[ Keys::maximumOrderOfBodyUndergoingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyUndergoingAcceleration_;
            jsonObject[ Keys::maximumDegreeOfCentralBody ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfCentralBody_;
            jsonObject[ Keys::maximumOrderOfCentralBody ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfCentralBody_;
            return;
        }
        case relativistic_correction_acceleration:
        {
            boost::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationCorrectionSettings
                    = boost::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >( accelerationSettings );
            jsonObject[ Keys::calculateSchwarzschildCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateSchwarzschildCorrection_;
            jsonObject[ Keys::calculateLenseThirringCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateLenseThirringCorrection_;
            jsonObject[ Keys::calculateDeSitterCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateDeSitterCorrection_;
            jsonObject[ Keys::primaryBody ] =
                    relativisticAccelerationCorrectionSettings->primaryBody_;
            jsonObject[ Keys::centralBodyAngularMomentum ] =
                    relativisticAccelerationCorrectionSettings->centralBodyAngularMomentum_;
            return;
        }
        case empirical_acceleration:
        {
            boost::shared_ptr< EmpiricalAccelerationSettings > empiricalAccelerationSettings
                    = boost::dynamic_pointer_cast< EmpiricalAccelerationSettings >( accelerationSettings );
            jsonObject[ Keys::constantAcceleration ] = empiricalAccelerationSettings->constantAcceleration_;
            jsonObject[ Keys::sineAcceleration ] = empiricalAccelerationSettings->sineAcceleration_;
            jsonObject[ Keys::cosineAcceleration ] = empiricalAccelerationSettings->cosineAcceleration_;
            return;
        }
        default:
            throw std::runtime_error( stringFromEnum( accelerationType, availableAccelerations )
                                      + " not supported by json_interface." );
        }
    }
}

//! Convert `json` to `AccelerationSettings`.
void from_json( const json& jsonObject, boost::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    accelerationSettings = json_interface::createAccelerationSettings( jsonObject );
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::AccelerationSettings > createAccelerationSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace basic_astrodynamics;
    using namespace simulation_setup;
    using Keys = Keys::Acceleration;

    // Get acceleration type
    const AvailableAcceleration accelerationType = getValue< AvailableAcceleration >( settings, keyTree + Keys::type );

    if ( accelerationType == third_body_point_mass_gravity ||
         accelerationType == third_body_spherical_harmonic_gravity ||
         accelerationType == third_body_mutual_spherical_harmonic_gravity )
    {
        std::cerr << "Whether a body will cause a third-body acceleration is determined internally "
                  << "by Tudat based on the propagation settings.\nRemove \"thirdBody\" from \""
                  << stringFromEnum( accelerationType, availableAccelerations )
                  << "\" at key " << keyTree
                  << " to silence this warning." << std::endl;
    }

    switch ( accelerationType ) {
    case undefined_acceleration:
    case aerodynamic:
    case cannon_ball_radiation_pressure:
        return boost::make_shared< AccelerationSettings >( accelerationType );
    case point_mass_gravity:
    case third_body_point_mass_gravity:
        return boost::make_shared< AccelerationSettings >( point_mass_gravity );
    case spherical_harmonic_gravity:
    case third_body_spherical_harmonic_gravity:
        return boost::make_shared< SphericalHarmonicAccelerationSettings >(
                    getNumeric< int >( settings, keyTree + Keys::maximumDegree ),
                    getNumeric< int >( settings, keyTree + Keys::maximumOrder ) );
    case mutual_spherical_harmonic_gravity:
    case third_body_mutual_spherical_harmonic_gravity:
    {
        MutualSphericalHarmonicAccelerationSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
        return boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                    getNumeric< int >( settings, keyTree + Keys::maximumDegreeOfBodyExertingAcceleration ),
                    getNumeric< int >( settings, keyTree + Keys::maximumOrderOfBodyExertingAcceleration ),
                    getNumeric< int >( settings, keyTree + Keys::maximumDegreeOfBodyUndergoingAcceleration ),
                    getNumeric< int >( settings, keyTree + Keys::maximumOrderOfBodyUndergoingAcceleration ),
                    getNumeric( settings, keyTree + Keys::maximumDegreeOfCentralBody,
                                defaults.maximumDegreeOfCentralBody_ ),
                    getNumeric( settings, keyTree + Keys::maximumOrderOfCentralBody,
                                defaults.maximumOrderOfCentralBody_ ) );
    }
    case relativistic_correction_acceleration:
    {
        RelativisticAccelerationCorrectionSettings defaults;
        return boost::make_shared< RelativisticAccelerationCorrectionSettings >(
                    getValue( settings, keyTree + Keys::calculateSchwarzschildCorrection,
                              defaults.calculateSchwarzschildCorrection_ ),
                    getValue( settings, keyTree + Keys::calculateLenseThirringCorrection,
                              defaults.calculateLenseThirringCorrection_ ),
                    getValue( settings, keyTree + Keys::calculateDeSitterCorrection,
                              defaults.calculateDeSitterCorrection_ ),
                    getValue( settings, keyTree + Keys::primaryBody,
                              defaults.primaryBody_ ),
                    getValue( settings, keyTree + Keys::centralBodyAngularMomentum,
                              defaults.centralBodyAngularMomentum_ ) );
    }
    case empirical_acceleration:
    {
        EmpiricalAccelerationSettings defaults;
        return boost::make_shared< EmpiricalAccelerationSettings >(
                    getValue( settings, keyTree + Keys::constantAcceleration, defaults.constantAcceleration_ ),
                    getValue( settings, keyTree + Keys::sineAcceleration, defaults.sineAcceleration_ ),
                    getValue( settings, keyTree + Keys::cosineAcceleration, defaults.cosineAcceleration_ ) );
    }
    default:
        throw std::runtime_error( stringFromEnum( accelerationType, availableAccelerations )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
