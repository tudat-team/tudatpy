/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      101111    E. Iorfida        Creation of code.
 *      101111    E. Iorfida        Implementation of all the equations up to the Newton method.
 *      101117    E. Iorfida        Velocities computations added.
 *      101126    E. Iorfida        Get/set codes deleted.
 *      101206    E. Iorfida        LambertTargetingElements class deleted,
 *                                  added setInitialState, modified punctuation. Set single
 *                                  variables, change variables names in more understandable ones.
 *      101209    E. Iorfida        Corrected some coding errors.
 *      101213    E. Iorfida        Deleted lambertAngle, added numberOfRevolution, modified
 *                                  implementation.
 *      101214    E. Iorfida        Implementation only for the case with numberOfRevolution = 0.
 *      110113    E. Iorfida        Added necessary elements to build pointer-to-member-function
 *                                  to RootFinderAlgorithms and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to use the last
 *                                  version of Newton-Raphson code.
 *      110126    E. Iorfida        Initialized member functions.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfdeletVelocityVector'
 *                                  became 'velocity'. Requested references to specific formulas. Also
 *                                  corrected 'tangential' to 'transverse'. Simplified computation
 *                                  of radial unit vector. Corrected computation of transverse
 *                                  heliocentric velocity.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names (from
 *                                  heliocentric, to inertial). Added patch for negative case of
 *                                  initialLambertGuess_. Added equations references.
 *      110206    E. Iorfida        Added unique function for Newton-Raphson method. Added
 *                                  computeAbsoluteValue to the initialLambertGuess_ for
 *                                  non-converging cases.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      110418    E. Iorfida        Added a new normal plane that take into account the case of two
 *                                  parallel position vector (with a relative angle of 180
 *                                  degrees). Better defined the pointers to the output
 *                                  CartesianVelocityElements.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120620    T. Secretin       Removed inheritance from TrajectoryDesignMethod class.
 *                                  Turned into base class for Lambert Targeters. Previous code
 *                                  adapted and moved to LambertTargeterGooding.cpp.
 *      120704    P. Musegaas       Moved getInertialVelocityVectors to the header.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"

namespace tudat
{
namespace mission_segments
{

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, LambertTargeter& lambertTargeter )
{
    stream << "The position vector at departure is set to: "
           << lambertTargeter.cartesianPositionAtDeparture
           << "The position vector at arrival is set to: "
           << lambertTargeter.cartesianPositionAtArrival
           << "The velocity vector at departure is computed as: "
           << lambertTargeter.getInertialVelocityAtDeparture( )
           << "The velocity vector at departure is computed as: "
           << lambertTargeter.getInertialVelocityAtArrival( ) << std::endl;

    // Return stream.
    return stream;
}

} // namespace mission_segments
} // namespace tudat
