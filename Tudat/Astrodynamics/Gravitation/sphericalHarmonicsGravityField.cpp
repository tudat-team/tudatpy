/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      101117    K. Kumar          File created.
 *      101214    K. Kumar          Updated getGradientOfPotential( ) and getLaplacianOfPotential( ).
 *      101215    K. Kumar          Simplified getGradientOfPotential( ) and
 *                                  getLaplacianOfPotential( ).
 *      101216    K. Kumar          Updated functions to use position of origin for relative
 *                                  position.
 *      110106    K. Kumar          Added set/get functions for degree and order of expansion.
 *      110202    K. Kumar          Updated code to make use of the CartesianPositionElements
 *                                  class.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      110805    K. Kumar          Added predefined functionality with WGS-72 and WGS-84 predefined
 *                                  predefined Earth gravity fields.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      141020    D. Dirkx          Change of architecture.
 *
 *    References
 *      Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. Revisiting Spacetrack Report #3:
 *          Rev 1, Proceedings of the AIAA/AAS Astrodynamics Specialist Conference. Keystone, CO,
 *          2006.
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
namespace tudat
{

namespace gravitation
{

//! Function to calculate the gravitational potential from a spherical harmonic field expansion.
double calculateSphericalHarmonicGravitationalPotential(
        const Eigen::Vector3d& bodyFixedPosition, const double gravitationalParameter,
        const double referenceRadius,
        const Eigen::MatrixXd& cosineCoefficients,
        const Eigen::MatrixXd& sineCoefficients,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const int minimumumDegree,
        const int minimumumOrder )
{
    // Initialize (distance/reference radius)^n (n=ratioToPowerDegree)
    double ratioToPowerDegree = 1.0;
    double radiusRatio = referenceRadius / bodyFixedPosition.norm( );

    // Declare local variables used in calculation
    double legendrePolynomial = 0.0;
    double singleDegreeTerm = 0.0;

    // Determine body fixed spherical position of body udnergoing acceleration.
    Eigen::Vector3d sphericalPositon =
            coordinate_conversions::convertCartesianToSpherical( bodyFixedPosition );
    double latitude = mathematical_constants::PI / 2.0 - sphericalPositon.y( );
    double longitude = sphericalPositon.z( );

    double potential = 0.0;
    int startDegree = 0;

    // Initialize value of potential to 1 (C_{0,0})
    if( minimumumDegree == 0 )
    {
        potential = 1.0;
        startDegree = 1;
    }
    else
    {
        startDegree = minimumumDegree;
        ratioToPowerDegree *= basic_mathematics::raiseToIntegerPower< double >( radiusRatio, startDegree - 1 );
    }

    basic_mathematics::LegendreCache& legendreCacheReference = *sphericalHarmonicsCache->getLegendreCache( );
    legendreCacheReference.update( std::sin( latitude ) );

    // Iterate over all degrees
    for( int degree = startDegree; degree < cosineCoefficients.rows( ); degree++ )
    {
        singleDegreeTerm = 0.0;

        // Iterate over all orders in current degree for which coefficients are provided.
        for( int order = minimumumOrder; ( order < cosineCoefficients.cols( ) &&
                                           order <= degree ); order++ )
        {
            // Calculate legendre polynomial (geodesy-normalized) at current degree and order
            legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomial(
                        degree, order, legendreCacheReference );

            // Calculate contribution to potential from current degree and order
            singleDegreeTerm += legendrePolynomial * ( cosineCoefficients( degree, order ) *
                                                       std::cos( order * longitude ) +
                                                       sineCoefficients( degree, order ) *
                                                       std::sin( order * longitude ) );
        }

        // Add potential contributions from current degree to toal value.
        ratioToPowerDegree *= radiusRatio;
        potential += singleDegreeTerm * ratioToPowerDegree;
    }

    // Multiply by central term and return
    return potential * gravitationalParameter / bodyFixedPosition.norm( );
}
}
}
