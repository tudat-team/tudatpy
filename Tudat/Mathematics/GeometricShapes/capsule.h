/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      E.H. Hirschel and C. Weiland, Selected Aerothermodynamic Design Problems of Hypersonic
 *          Flight Vehicles (chapter 5), Springer/AIAA, 2009.
 *      D. Dirkx, Continuous Shape Optimization of Entry Vehicles, MSc thesis, Delft University
 *          of Technology, 2011 (Unpublished).
 *
 */

#ifndef TUDAT_CAPSULE_H
#define TUDAT_CAPSULE_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"

namespace tudat
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
     * \param stream Stream to which info is to be printed.
     * \param capsule Capsule of which info is to be printed.
     * \return Stream with printed info.
     */
    friend std::ostream &operator << ( std::ostream &stream, Capsule& capsule );

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

//! Typedef for shared-pointer to Capsule object.
typedef boost::shared_ptr< Capsule > CapsulePointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_CAPSULE_H
