/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110617    F.M. Engelen      First creation of code.
 *      110822    D.Dirkx           Removed pointer to double member, removed location of center
 *                                  of mass, minor changes.
 *
 *    References
 *
 */

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"

namespace tudat
{
namespace astrodynamics
{
namespace moment_models
{

//! Compute the aerodynamic moment in same reference frame as input coefficients.
Eigen::MatrixXd computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const double referenceLength,
                                          const Eigen::Vector3d& momentCoefficients )
{
    return dynamicPressure * referenceArea * referenceLength * momentCoefficients;
}

//! Calculates the aerodynamic moment in same reference frame as input coefficients.
Eigen::MatrixXd computeAerodynamicMoment(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface )
{
    return computeAerodynamicMoment( dynamicPressure,
                                     coefficientInterface.getReferenceArea( ),
                                     coefficientInterface.getReferenceLength( ),
                                     coefficientInterface.getCurrentMomentCoefficients( ) );
}

} // namespace moment_models
} // namespace astrodynamics
} // namespace tudat
