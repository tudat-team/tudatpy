/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TORQUESETTINGS_H
#define TUDAT_TORQUESETTINGS_H

#include "tudat/io/serialization/core.h"
#include "tudat/io/serialization/file_io_declarations.h"

#include "tudat/astro/basic_astro/torqueModelTypes.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"

namespace tudat
{

namespace simulation_setup
{

// Class for providing settings for torque model.
/*!
 *  Class for providing settings for torque model. This class is a functional (base) class for
 *  settings of torque models that  require no information in addition to their type.
 *  Classes defining settings for torque models requiring additional information must be
 *  derived from this class.
 *  Bodies exerting and undergong torque are set externally from this class.
 *  This class can be used for the easy setup of torque models
 *  (see createTorqueModels.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */

//! @get_docstring(TorqueSettings.__docstring__)
class TorqueSettings
{
public:
    // Constructor, sets type of torque.
    /*!
     *  Constructor, sets type of torque.
     *  \param torqueType Type of torque from AvailableTorque enum.
     */
    TorqueSettings( const basic_astrodynamics::AvailableTorque torqueType ): torqueType_( torqueType ) {}

    // Destructor
    virtual ~TorqueSettings( ) {}

    // Type of torque that is to be created.
    basic_astrodynamics::AvailableTorque torqueType_;

    // Used for serialization testing
    bool operator==( const TorqueSettings& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const TorqueSettings& rhs ) const
    {
        return !equals( rhs );
    }

    //! Save torque settings to a JSON file
    TUDAT_DECLARE_FILE_IO_POLYMORPHIC( TorqueSettings )

protected:
    // Default constructor for serialization
    TorqueSettings( ): torqueType_( basic_astrodynamics::underfined_torque ) {}

    // Each derived class should implement this function such that it returns true if a deserialized object is
    // equal to the original object.
    virtual bool equals( const TorqueSettings& rhs ) const
    {
        return torqueType_ == rhs.torqueType_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( torqueType_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( torqueType_ ) );
    }
};

// Class to define settings for a spherical harmonic gravitational torque exerted by a point mass.
//! @get_docstring(SphericalHarmonicTorqueSettings.__docstring__)
class SphericalHarmonicTorqueSettings : public TorqueSettings
{
public:
    // Constructor
    /*!
     * Constructor
     * \param maximumDegree Maximum degree to which gravity field of body undergoing torque is to be exerted
     * \param maximumOrder Maximum order to which gravity field of body undergoing torque is to be exerted
     */
    SphericalHarmonicTorqueSettings( const int maximumDegree, const int maximumOrder ):
        TorqueSettings( basic_astrodynamics::spherical_harmonic_gravitational_torque ), maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder )
    {}

    // Maximum degree to which gravity field of body undergoing torque is to be exerted
    int maximumDegree_;

    // Maximum order to which gravity field of body undergoing torque is to be exerted
    int maximumOrder_;

protected:
    // Default constructor for serialization
    SphericalHarmonicTorqueSettings( ): maximumDegree_( 0 ), maximumOrder_( 0 ) {}

    bool equals( const TorqueSettings& rhs ) const override
    {
        const auto* derived = dynamic_cast< const SphericalHarmonicTorqueSettings* >( &rhs );
        if( !derived ) return false;
        if( !TorqueSettings::equals( rhs ) ) return false;
        return maximumDegree_ == derived->maximumDegree_ && maximumOrder_ == derived->maximumOrder_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
        ar( CEREAL_NVP( maximumDegree_ ) );
        ar( CEREAL_NVP( maximumOrder_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
        ar( CEREAL_NVP( maximumDegree_ ) );
        ar( CEREAL_NVP( maximumOrder_ ) );
    }
};

class FullTwoBodySphericalHarmonicTorqueSettings : public TorqueSettings
{
public:
    FullTwoBodySphericalHarmonicTorqueSettings(
            const std::shared_ptr< AccelerationSettings > fullTwoBodySphericalHarmonicAccelerationSettings ):
        TorqueSettings( basic_astrodynamics::full_two_body_spherical_harmonic_gravitational_torque ),
        fullTwoBodySphericalHarmonicAccelerationSettings_( fullTwoBodySphericalHarmonicAccelerationSettings )
    {}

    std::shared_ptr< AccelerationSettings > fullTwoBodySphericalHarmonicAccelerationSettings_;

protected:
    // Default constructor for serialization
    FullTwoBodySphericalHarmonicTorqueSettings( ): fullTwoBodySphericalHarmonicAccelerationSettings_( nullptr ) {}

    bool equals( const TorqueSettings& rhs ) const override
    {
        const auto* derived = dynamic_cast< const FullTwoBodySphericalHarmonicTorqueSettings* >( &rhs );
        if( derived == nullptr || !TorqueSettings::equals( rhs ) )
        {
            return false;
        }
        return ( fullTwoBodySphericalHarmonicAccelerationSettings_ == derived->fullTwoBodySphericalHarmonicAccelerationSettings_ ) ||
                ( fullTwoBodySphericalHarmonicAccelerationSettings_ != nullptr &&
                  derived->fullTwoBodySphericalHarmonicAccelerationSettings_ != nullptr &&
                  *fullTwoBodySphericalHarmonicAccelerationSettings_ == *derived->fullTwoBodySphericalHarmonicAccelerationSettings_ );
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
        ar( CEREAL_NVP( fullTwoBodySphericalHarmonicAccelerationSettings_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
        ar( CEREAL_NVP( fullTwoBodySphericalHarmonicAccelerationSettings_ ) );
    }
};

class FourthDegreeFullTwoBodyGravitationalTorqueSettings : public TorqueSettings
{
public:
    FourthDegreeFullTwoBodyGravitationalTorqueSettings( ):
        TorqueSettings( basic_astrodynamics::fourth_degree_full_two_body_gravitational_torque )
    {}

protected:
    bool equals( const TorqueSettings& rhs ) const override
    {
        return dynamic_cast< const FourthDegreeFullTwoBodyGravitationalTorqueSettings* >( &rhs ) != nullptr &&
                TorqueSettings::equals( rhs );
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< TorqueSettings >( this ) );
    }
};

inline Eigen::Vector3d applyTorqueScalingFunction( const std::function< Eigen::Vector3d( const double ) > torqueFunction,
                                                   const std::function< double( const double ) > scalingFunction,
                                                   const double time )
{
    return torqueFunction( time ) * scalingFunction( time );
}

//! @get_docstring(CustomTorqueSettings.__docstring__)
class CustomTorqueSettings : public TorqueSettings
{
public:
    CustomTorqueSettings( const std::function< Eigen::Vector3d( const double ) > torqueFunction ):
        TorqueSettings( basic_astrodynamics::custom_torque ), torqueFunction_( torqueFunction )
    {}

    CustomTorqueSettings( const std::function< Eigen::Vector3d( const double ) > torqueFunction,
                          std::function< double( const double ) > scalingFunction ):
        TorqueSettings( basic_astrodynamics::custom_torque ),
        torqueFunction_( std::bind( &applyTorqueScalingFunction, torqueFunction, scalingFunction, std::placeholders::_1 ) )
    {}

    std::function< Eigen::Vector3d( const double ) > torqueFunction_;

protected:
    // Default constructor for serialization
    CustomTorqueSettings( ): torqueFunction_( nullptr ) {}

    bool equals( const TorqueSettings& rhs ) const override
    {
        const auto* derived = dynamic_cast< const CustomTorqueSettings* >( &rhs );
        if( !derived ) return false;
        if( !TorqueSettings::equals( rhs ) ) return false;
        // std::function cannot be compared; always return false for custom torque
        return false;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        static_cast< void >( ar );
        throw std::runtime_error(
                "CustomTorqueSettings cannot be serialized: std::function member 'torqueFunction_' is not serializable." );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        static_cast< void >( ar );
        throw std::runtime_error(
                "CustomTorqueSettings cannot be serialized: std::function member 'torqueFunction_' is not serializable." );
    }
};

//! @get_docstring(aerodynamicTorque)
inline std::shared_ptr< TorqueSettings > aerodynamicTorque( )
{
    return std::make_shared< TorqueSettings >( basic_astrodynamics::aerodynamic_torque );
}

inline std::shared_ptr< TorqueSettings > radiationPressureTorque( )
{
    return std::make_shared< TorqueSettings >( basic_astrodynamics::radiation_pressure_torque );
}

//! @get_docstring(secondDegreeGravitationalTorque)
inline std::shared_ptr< TorqueSettings > secondDegreeGravitationalTorque( )
{
    return std::make_shared< TorqueSettings >( basic_astrodynamics::second_order_gravitational_torque );
}

//! @get_docstring(sphericalHarmonicGravitationalTorque)
inline std::shared_ptr< TorqueSettings > sphericalHarmonicGravitationalTorque( const int maximumDegree, const int maximumOrder )
{
    return std::make_shared< SphericalHarmonicTorqueSettings >( maximumDegree, maximumOrder );
}

inline std::shared_ptr< TorqueSettings > fullTwoBodySphericalHarmonicGravitationalTorque( const int maximumDegreeOfBodyUndergoingTorque,
                                                                                          const int maximumOrderOfBodyUndergoingTorque,
                                                                                          const int maximumDegreeOfBodyExertingTorque,
                                                                                          const int maximumOrderOfBodyExertingTorque )
{
    return std::make_shared< FullTwoBodySphericalHarmonicTorqueSettings >(
            fullTwoBodySphericalHarmonicAcceleration( maximumDegreeOfBodyUndergoingTorque,
                                                      maximumOrderOfBodyUndergoingTorque,
                                                      maximumDegreeOfBodyExertingTorque,
                                                      maximumOrderOfBodyExertingTorque ) );
}

inline std::shared_ptr< TorqueSettings > fullTwoBodySphericalHarmonicGravitationalTorque(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    return std::make_shared< FullTwoBodySphericalHarmonicTorqueSettings >(
            std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( coefficientCombinationsToUse ) );
}

inline std::shared_ptr< TorqueSettings > fourthDegreeFullTwoBodyGravitationalTorque( )
{
    return std::make_shared< FourthDegreeFullTwoBodyGravitationalTorqueSettings >( );
}

inline std::shared_ptr< TorqueSettings > dissipativeTorque( )
{
    return std::make_shared< TorqueSettings >( basic_astrodynamics::dissipative_torque );
}

//! @get_docstring(customTorqueSettings)
inline std::shared_ptr< TorqueSettings > customTorqueSettings( const std::function< Eigen::Vector3d( const double ) > torqueFunction,
                                                               const std::function< double( const double ) > scalingFunction = nullptr )
{
    if( scalingFunction == nullptr )
    {
        return std::make_shared< CustomTorqueSettings >( torqueFunction );
    }
    else
    {
        return std::make_shared< CustomTorqueSettings >( torqueFunction, scalingFunction );
    }
}

typedef std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< TorqueSettings > > > > SelectedTorqueMap;

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_TORQUESETTINGS_H
