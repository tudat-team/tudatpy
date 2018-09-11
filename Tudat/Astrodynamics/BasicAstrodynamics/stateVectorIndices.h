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

#ifndef TUDAT_STATE_INDICES_H
#define TUDAT_STATE_INDICES_H

namespace tudat
{

namespace orbital_element_conversions
{

//! Cartesian elements indices.
enum CartesianElementIndices
{
    xCartesianPositionIndex = 0,
    yCartesianPositionIndex = 1,
    zCartesianPositionIndex = 2,
    xCartesianVelocityIndex = 3,
    yCartesianVelocityIndex = 4,
    zCartesianVelocityIndex = 5
};

//! Keplerian elements indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex = 0,
    eccentricityIndex = 1,
    inclinationIndex = 2,
    argumentOfPeriapsisIndex = 3,
    longitudeOfAscendingNodeIndex = 4,
    trueAnomalyIndex = 5,
    semiLatusRectumIndex = 0
};

//! Modified equinoctial element vector indices.
enum ModifiedEquinoctialElementVectorIndices
{
    semiParameterIndex = 0,
    fElementIndex = 1,
    gElementIndex = 2,
    hElementIndex = 3,
    kElementIndex = 4,
    trueLongitudeIndex = 5
};

//! Spherical orbital state element indices
enum SphericalOrbitalStateElementIndices
{
    radiusIndex = 0,
    latitudeIndex = 1,
    longitudeIndex = 2,
    speedIndex = 3,
    flightPathIndex = 4,
    headingAngleIndex = 5
};

//! Unified state model with quaternions indices.
enum UnifiedStateModelQuaternionsElementIndices
{
    CHodographUSM7Index = 0,
    Rf1HodographUSM7Index = 1,
    Rf2HodographUSM7Index = 2,
    etaUSM7Index = 3,
    epsilon1USM7Index = 4,
    epsilon2USM7Index = 5,
    epsilon3USM7Index = 6
};

//! Unified state model with modified Rodrigues parameters indices.
enum UnifiedStateModelModifiedRodriguesParametersElementIndices
{
    CHodographUSM6Index = 0,
    Rf1HodographUSM6Index = 1,
    Rf2HodographUSM6Index = 2,
    sigma1USM6Index = 3,
    sigma2USM6Index = 4,
    sigma3USM6Index = 5,
    shadowFlagUSM6Index = 6
};

//! Unified state model with exponential map indices.
enum UnifiedStateModelExponentialMapElementIndices
{
    CHodographUSMEMIndex = 0,
    Rf1HodographUSMEMIndex = 1,
    Rf2HodographUSMEMIndex = 2,
    e1USMEMIndex = 3,
    e2USMEMIndex = 4,
    e3USMEMIndex = 5,
    shadowFlagUSMEMIndex = 6
};

//! Quaternions indices.
enum QuaternionsElementIndices
{
    etaQuaternionIndex = 0,
    epsilon1QuaternionIndex = 1,
    epsilon2QuaternionIndex = 2,
    epsilon3QuaternionIndex = 3
};

//! Modified Rodrigues parameters indices.
enum ModifiedRodriguesParametersElementIndices
{
    sigma1ModifiedRodriguesParametersIndex = 0,
    sigma2ModifiedRodriguesParametersIndex = 1,
    sigma3ModifiedRodriguesParametersIndex = 2,
    shadowFlagModifiedRodriguesParametersIndex = 3
};

//! Exponential map indices.
enum ExponentialMapElementIndices
{
    e1ExponentialMapIndex = 0,
    e2ExponentialMapIndex = 1,
    e3ExponentialMapIndex = 2,
    shadowFlagExponentialMapIndex = 3
};

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_STATE_INDICES_H
