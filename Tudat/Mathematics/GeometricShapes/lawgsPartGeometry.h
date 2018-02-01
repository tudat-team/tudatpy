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
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard (LaWGS) format, NASA
 *          TECHNICAL MEMORANDUM 85767.
 *
 */

#ifndef TUDAT_LAWGS_PART_GEOMETRY_H
#define TUDAT_LAWGS_PART_GEOMETRY_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/GeometricShapes/quadrilateralMeshedSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Class to define a mesh accoring to Lawgs standards.
/*!
 * Class to define a surface mesh accoring to the Langley Wireframe Geometry
 * standard, seet reference.
 */
class LawgsPartGeometry : public QuadrilateralMeshedSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LawgsPartGeometry( ) : name_( "" ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~LawgsPartGeometry( ) { }

    //! Create a mesh surface on a single-surface geometry.
    /*!
     * Creates a mesh surface on a single-surface geometry. A meshed surface based on the given
     * surface geometry is created.
     * \param originalSurface surface from which a mesh is to be created.
     * \param numberOfLines Number of points to be sampled from 1st
     *          independent variable.
     * \param numberOfPoints Number of points to be sampled from 2nd
     *          independent variable.
     */
    void setMesh( boost::shared_ptr< SingleSurfaceGeometry > originalSurface,
                  int numberOfLines, int numberOfPoints );

    //! Copy constructor.
    /*!
     * Copy constructor to deep-copy contents of a LawgsPartGeomtry object to a new one.
     * \param partToCopy LawgsPartGeometry to copy
     */
     LawgsPartGeometry( const LawgsPartGeometry& partToCopy );

    //! Function to retrieve surface point.
    /*!
     * This function retrieves a surface point from on a panel. Because four
     * points define a panel in this mesh format, the collection of panels
     * will not be necessarilly watertight, making this function non-trivial,
     * since the four points defining the panel will not necessarily lie in
     * the same plane.
     * \param independentVariable1 Independent variable indicating on (or
     * between) which lines to retrieve a point.
     * \param independentVariable2 Independent variable indicating on (or
     * between) which points to retrieve a point.
     * \return point on mesh panel.
     */
     virtual Eigen::VectorXd getSurfacePoint( const double independentVariable1,
                                              const double independentVariable2 );

    //! Get surface derivative (currently not implemented).
    /*!
     *  Currently unavailable function to return surface derivative.
     *  \param u NOTE: Function unavailable.
     *  \param v NOTE: Function unavailable.
     *  \param uDerivative NOTE: Function unavailable.
     *  \param vDerivative NOTE: Function unavailable.
     *  \return  NOTE: Function unavailable.
     */
    virtual Eigen::VectorXd getSurfaceDerivative( const double u, const double v,
                                                  const int uDerivative, const int vDerivative );

    //! Set name of a Lawgs part.
    /*!
     *  Sets the name of a Lawgs part.
     *  \param name New name of a Lawgs part.
     */
    void setName( const std::string& name ) { name_ = name; }

    //! Get parameter.
    /*!
     *  Function not currently impemented for Lawgs.
     *  \param i Function not currently impemented for Lawgs.
     *  \return Function not currently impemented for Lawgs.
     */
    virtual double getParameter( const int i );

    //! Set parameter.
    /*!
     * Sets parameter.
     *  Function not currently impemented for Lawgs.
     *  \param parameterIndex Function not currently impemented for Lawgs.
     *  \param value Function not currently impemented for Lawgs.
     *  \return Function not currently impemented for Lawgs.
     */
    virtual void setParameter( const int parameterIndex, const double value );

    //! Get part name.
    /*!
     *  Returns part name.
     *  \return Part name.
     */
    std::string getName( ) { return name_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the number of
     * lines and points, and the name of the part.
     * \param stream Stream object.
     * \param lawgsPartGeometry Lawgs part geometry.
     * \return Stream object.
     */
    friend std::ostream& operator << ( std::ostream& stream, LawgsPartGeometry& lawgsPartGeometry );

protected:

    //! Part name.
    /*!
     * Part name.
     */
    std::string name_;

private:
};

//! Typedef for shared-pointer to LawgsPartGeometry object.
typedef boost::shared_ptr< LawgsPartGeometry > LawgsPartGeometryPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_LAWGS_PART_GEOMETRY_H
