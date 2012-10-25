#ifndef TUDAT_STATE_INDICES_H
#define TUDAT_STATE_INDICES_H

namespace tudat
{
namespace basic_astrodynamics
{

//! Keplerian elements indices.
enum CartesianElementIndices
{
    xCartesianPositionIndex,
    yCartesianPositionIndex,
    zCartesianPositionIndex,
    xCartesianVelocityIndex,
    yCartesianVelocityIndex,
    zCartesianVelocityIndex
};

//! Keplerian elements indices.
enum KeplerianElementIndices
{
    semiMajorAxisIndex,
    eccentricityIndex,
    inclinationIndex,
    argumentOfPeriapsisIndex,
    longitudeOfAscendingNodeIndex,
    trueAnomalyIndex
};

} // namespace basic_astrodynamics
} // namespace tudat

#endif // TUDAT_STATE_INDICES_H
