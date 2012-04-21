/*   Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110822    D. Dirkx          Removed pointer to double member, minor changes.
 *      120209    D. Dirkx          Changed force models to free functions, moved implementation
 *                                  to source file.
 *      120324    K. Kumar          Added astrodynamics namespace layer; updated function input
 *                                  arguments.
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"

namespace tudat
{
namespace astrodynamics
{
namespace force_models
{

//! Compute the aerodynamic force in same reference frame as input coefficients.
Eigen::VectorXd computeAerodynamicForce( const double dynamicPressure,
                                         const double referenceArea,
                                         const Eigen::Vector3d& aerodynamicCoefficients )
{
    return dynamicPressure * referenceArea * aerodynamicCoefficients;
}

//! Compute the aerodynamic force in same reference frame as input coefficients.
Eigen::MatrixXd computeAerodynamicForce(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface )
{
    return computeAerodynamicForce( dynamicPressure,
                                    coefficientInterface.getReferenceArea( ),
                                    coefficientInterface.getCurrentForceCoefficients( ) );
}

} // namespace force_models
} // namespace astrodynamics
} // namespace tudat
