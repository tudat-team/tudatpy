/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATORTYPE_H
#define TUDAT_PROPAGATORTYPE_H

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>

#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

struct PropagatorType {
    PropagatorType( const TranslationalPropagatorType translationalPropagatorType ):
        translationalPropagatorType_( translationalPropagatorType )
    {}

    PropagatorType( const RotationalPropagatorType rotationalPropagatorType ): rotationalPropagatorType_( rotationalPropagatorType ) {}

    PropagatorType( ): otherPropagator_( true ) {}

    PropagatorType( const int customStateSize ): customStateSize_( customStateSize ) {}

    TranslationalPropagatorType translationalPropagatorType_ = undefined_translational_propagator;
    RotationalPropagatorType rotationalPropagatorType_ = undefined_rotational_propagator;
    bool otherPropagator_ = false;
    int customStateSize_ = 0;

private:
    friend class cereal::access;
    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( translationalPropagatorType_ ) );
        ar( CEREAL_NVP( rotationalPropagatorType_ ) );
        ar( CEREAL_NVP( otherPropagator_ ) );
        ar( CEREAL_NVP( customStateSize_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( translationalPropagatorType_ ) );
        ar( CEREAL_NVP( rotationalPropagatorType_ ) );
        ar( CEREAL_NVP( otherPropagator_ ) );
        ar( CEREAL_NVP( customStateSize_ ) );
    }

public:
    //! Equality comparison for PropagatorType
    bool operator==( const PropagatorType& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const PropagatorType& rhs ) const
    {
        return !( *this == rhs );
    }

    //! Equality comparison via equals method
    bool equals( const PropagatorType& rhs ) const
    {
        return translationalPropagatorType_ == rhs.translationalPropagatorType_ &&
                rotationalPropagatorType_ == rhs.rotationalPropagatorType_ && otherPropagator_ == rhs.otherPropagator_ &&
                customStateSize_ == rhs.customStateSize_;
    }

private:
};

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_PROPAGATORTYPE_H
