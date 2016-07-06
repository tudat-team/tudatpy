#ifndef RELATIVISTICLIGHTTIMECORRECTIONS_H
#define RELATIVISTICLIGHTTIMECORRECTIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include <cmath>
#include <vector>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace relativity
{

//! Function to calculate first order relativistic light time correction due to a gravitating point mass.
/*!
 *  Function to calculate first order relativistic light time correction due to a gravitating point mass, according to Eq. (11.17) of 2010 IERS
 *  conventions.
 *  \param bodyGravitationalParameter Gravitational parameter of gravitating body.
 *  \param transmitterPosition Position of origin of electromagnetic signal (at time of transmission).
 *  \param receiverPosition Position of target of electromagentic signal (at time of reception)
 *  \param ppnParameterGamma Parametric post-Newtonian parameter gamma, a measure for the space-time curvature due to a unit rest mass (1.0 in GR)
 *  \return Light time correction (in seconds) due to the gravitating point mass.
 */
double calculateFirstOrderLightTimeCorrectionFromCentralBody( const double bodyGravitationalParameter,
                                                              const Eigen::Vector3d& transmitterPosition,
                                                              const Eigen::Vector3d& receiverPosition,
                                                              const Eigen::Vector3d& centralBodyPosition,
                                                              const double ppnParameterGamma = 1.0 );
}

}
#endif // RELATIVISTICLIGHTTIMECORRECTIONS_H
