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
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; corrected Doxygen comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor.
 *
 *    References
 *      E.H. Hirschel and C. Weiland, Selected Aerothermodynamic Design Problems
 *          of Hypersonic Flight Vehicles (chapter 5), Springer/AIAA, 2009.
 *      D. Dirkx, Continuous Shape Optimization of Entry Vehicles, MSc thesis,
 *          Delft University of Technology, 2011 (Unpublished).
 *
 */

#ifndef TUDAT_CAPSULE_H
#define TUDAT_CAPSULE_H

#include <iostream>

#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"

namespace tudat
{
namespace mathematics
{
namespace geometric_shapes
{

//! Capsule class.
/*!
 * Class that defines a Capsule shape, like that of the Apollo or ARD (Advanced
 * Re-entry Demonstrator) capsule. The class derives from the
 * CompositeSurfaceGeometry class, its constituent surface geometries are
 * a SphereSegment object for the nose region, a ConicalFrustum
 * object for the afterbody, a SphereSegment object for the
 * rear cap, and a Torus object for the connection between the nose and torus
 * parts. The setCapsule( ) function needs to be called after each of the
 * parameters has been set to initialize the capsule.
 */
class Capsule: public CompositeSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     * Default constructor, initializes single and composite surface lists.
     * \param noseRadius Radius of curvature of capsule nose.
     * \param middleRadius Maximum radius of capsule
     * \param rearLength Length of frustum section
     * \param rearAngle Cone half angle of frustum part
     * \param sideRadius Radius of curvature of capsule shoulder.
     */
    Capsule( const double noseRadius,
             const double middleRadius,
             const double rearLength,
             const double rearAngle,
             const double sideRadius );

    //! Get nose radius.
    /*!
     * Returns the nose radius of the capsule.
     * \return Nose radius.
     */
    double getNoseRadius( ) { return noseRadius_; }

    //! Get middle radius.
    /*!
     * Returns the middle radius of the capsule.
     * \return Middle radius.
     */
    double getMiddleRadius( ) { return middleRadius_; }

    //! Get rear length.
    /*!
     * Returns the rear length of the capsule.
     * \return Rear length.
     */
    double getRearLength( ) { return rearLength_; }

    //! Get rear angle.
    /*!
     * Returns the rear angle.
     * \return Rear angle.
     */
    double getRearAngle( ) { return rearAngle_; }

    //! Get side radius.
    /*!
     * Returns the side radius.
     * \return side radius.
     */
    double getSideRadius( ) { return sideRadius_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * and number defining parameters.
     */
    friend std::ostream &operator<<( std::ostream &stream, Capsule& capsule );

protected:

private:

    //! Middle radius.
    /*!
     * Middle radius.
     */
    double middleRadius_;

    //! Nose radius.
    /*!
     * Nose radius.
     */
    double noseRadius_;

    //! Rear length.
    /*!
     * Rear length.
     */
    double rearLength_;

    //! Side radius.
    /*!
     * Side radius.
     */
    double sideRadius_;

    //! Rear angle.
    /*!
     * Rear angle.
     */
    double rearAngle_;
};


} // namespace geometric_shapes
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_CAPSULE_H
