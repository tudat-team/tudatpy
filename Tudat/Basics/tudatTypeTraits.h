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

#ifndef TUDAT_TYPE_TRAITS_H
#define TUDAT_TYPE_TRAITS_H

#include <type_traits>

#include "Tudat/Basics/timeType.h"

namespace tudat
{

/// Implementation of \p std::enable_if_t.
/**
 * Implementation of \p std::enable_if_t, from C++14. See: http://en.cppreference.com/w/cpp/types/enable_if.
 */
template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

namespace is_eigen_matrix_detail
{
template< typename T >
std::true_type test( const Eigen::MatrixBase< T >* );

std::false_type test( ... );
}

template< typename T >
struct is_eigen_matrix: public decltype( is_eigen_matrix_detail::test( std::declval< T* >( ) ) )
{ };

//template< typename T,
//          enable_if_t< is_eigen_matrix< T >::value, int > = 0 >
//struct is_eigen_vector
//{
//    static const int value = ( T::ColsAtCompileTime == 1 );
//};

//template< typename T,
//          enable_if_t< !is_eigen_matrix< T >::value, int > = 0 >
//struct is_eigen_vector
//{
//    static const int value = false;
//};

template< typename T >
struct is_state_scalar {
  static const bool value = false;
};

template< >
struct is_state_scalar< double > {
  static const bool value = true;
};

template< >
struct is_state_scalar< long double > {
  static const bool value = true;
};

template< typename T >
struct is_time_type {
  static const bool value = false;
};

template< >
struct is_time_type< double > {
  static const bool value = true;
};

template< >
struct is_time_type< Time > {
  static const bool value = true;
};

template< typename StateScalarType, typename TimeType >
struct is_state_scalar_and_time_type {
  static const bool value = ( is_time_type< TimeType >::value && is_state_scalar< StateScalarType >::value );
};

template< typename T >
struct is_direct_gravity_partial {
    static const bool value = false;
};

namespace acceleration_partials
{
class CentralGravitationPartial;
}

template< >
struct is_direct_gravity_partial< acceleration_partials::CentralGravitationPartial > {
    static const bool value = true;
};

namespace acceleration_partials
{
class SphericalHarmonicsGravityPartial;
}

template< >
struct is_direct_gravity_partial< acceleration_partials::SphericalHarmonicsGravityPartial > {
    static const bool value = true;
};

namespace acceleration_partials
{
class MutualSphericalHarmonicsGravityPartial;
}

template< >
struct is_direct_gravity_partial< acceleration_partials::MutualSphericalHarmonicsGravityPartial > {
    static const bool value = true;
};

template< typename T >
struct is_direct_gravity_acceleration {
    static const bool value = false;
};

namespace gravitation
{
template< typename ReturnType >
class CentralGravitationalAccelerationModel;
}

template< >
struct is_direct_gravity_acceleration< gravitation::CentralGravitationalAccelerationModel< Eigen::Vector3d > > {
    static const bool value = true;
};

namespace gravitation
{
class SphericalHarmonicsGravitationalAccelerationModel;
}

template< >
struct is_direct_gravity_acceleration< gravitation::SphericalHarmonicsGravitationalAccelerationModel > {
    static const bool value = true;
};

namespace gravitation
{
class MutualSphericalHarmonicsGravitationalAccelerationModel;
}

template< >
struct is_direct_gravity_acceleration< gravitation::MutualSphericalHarmonicsGravitationalAccelerationModel > {
    static const bool value = true;
};



} // namespace tudat

#endif // TUDAT_TYPE_TRAITS_H
