/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      130806    J. Geul           File created.
 *
 *    References
 *
 *    Notes  
 *        
 *      Previously the cartesian stateDerivative was calculated by combining the second half 
 *      of the state with the accelerations. This has the benefit of handling any dimensional 
 *      vector. It does however make unnoted assumptions of the order of the vector elements. 
 *      Therefore a fixed size method is implemented below. Different dimensional problems
 *      must therefore use different methods, or stick to a 6D state scheme.
 *
 */

#ifndef TUDAT_STATEDERIVATIVE_MAP_CARTESIAN_H
#define TUDAT_STATEDERIVATIVE_MAP_CARTESIAN_H

#include <cmath>
#include <Eigen/Core>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>

namespace tudat
{
namespace state_derivative_models
{

    //! XYZ accelerations to Cartesian state derivatives.
    /*!
     * Uses the accelerations in a 6D state, such that the system is
     * reduced to a first-order ODE.
     *
     * \param state of the spacecraft in Cartesian elements.
     * \param acceleration in inertial XYZ frame of all accelerations [m/s^2].
     */
    Eigen::Matrix< double, 6, 1 > stateDerivativeMapCartesian( 
	const Eigen::Matrix< double, 6, 1 >& state, 
	const Eigen::Vector3d& acceleration )
    {

	// Index of Cartesian Coordinates
	using basic_astrodynamics::xCartesianPositionIndex;
	using basic_astrodynamics::yCartesianPositionIndex;
	using basic_astrodynamics::zCartesianPositionIndex;
	using basic_astrodynamics::xCartesianVelocityIndex;
	using basic_astrodynamics::yCartesianVelocityIndex;
	using basic_astrodynamics::zCartesianVelocityIndex;

	// Indices of Cartesian Acceleration
	using basic_astrodynamics::xCartesianAccelerationIndex;
	using basic_astrodynamics::yCartesianAccelerationIndex;
	using basic_astrodynamics::zCartesianAccelerationIndex;

	// Declare Orbital state derivative of the same size as Orbital state.
	Eigen::Matrix< double, 6, 1 > stateDerivative = state * 0;

	stateDerivative( xCartesianPositionIndex ) = state( xCartesianVelocityIndex );
	stateDerivative( yCartesianPositionIndex ) = state( yCartesianVelocityIndex );
	stateDerivative( zCartesianPositionIndex ) = state( zCartesianVelocityIndex );
	stateDerivative( xCartesianVelocityIndex ) = acceleration( xCartesianAccelerationIndex );
	stateDerivative( yCartesianVelocityIndex ) = acceleration( yCartesianAccelerationIndex );
	stateDerivative( zCartesianVelocityIndex ) = acceleration( zCartesianAccelerationIndex );

	return stateDerivative;

    }

} // namespace state_derivative_models
} // namespace tudat


#endif // TUDAT_STATEDERIVATIVE_MAP_CARTESIAN
