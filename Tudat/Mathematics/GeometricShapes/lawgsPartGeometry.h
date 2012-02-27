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
 *      101125    D. Dirkx          First version of file.
 *      110127    D. Dirkx          Finalized for code check.
 *      110206    J. Melman         Minor formatting issues.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard
 *          (LaWGS) format, NASA TECHNICAL MEMORANDUM 85767.
 *
 */

#ifndef TUDAT_LAWGS_PART_GEOMETRY_H
#define TUDAT_LAWGS_PART_GEOMETRY_H

#include <Eigen/Core>
#include <iostream>
#include "Tudat/Mathematics/GeometricShapes/quadrilateralMeshedSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{

//! Class to defined a mesh accoring to Lawgs standards.
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
    void setMesh( SingleSurfaceGeometry* originalSurface, int numberOfLines, int numberOfPoints );

    //! Copy constructor.
    /*!
     * Copy constructor to deep-copy contents of a LawgsPartGeomtry object to a new one.
     */
     LawgsPartGeometry( const LawgsPartGeometry& PartIn );

    //! Function to retrieve surface point.
    /*!
     * This function retrieves a surface point from on a panel. Because four
     * points define a panel in this mesh format, the collection of panels
     * will not be necessarilly watertight, making this function non-trivial,
     * since the four points defining the panle will not necessarily lie in
     * the same plane.
     * \param independentVariable1 Independent variable indicating on (or
     * between) which lines to retrieve a point.
     * \param independentVariable2 Independent variable indicating on (or
     * between) which points to retrieve a point.
     * \return point on mesh panel.
     */
     virtual Eigen::VectorXd getSurfacePoint( double independentVariable1,
                                              double independentVariable2 );

    //! Get surface derivative (currently not implemented).
    /*!
     * Currently unavailable function to return surface derivative.
     */
    virtual Eigen::VectorXd getSurfaceDerivative( double u, double v,
                                                  int uDerivative, int vDerivative );

    //! Set name of a Lawgs part.
    /*!
     * Sets the name of a Lawgs part.
     */
    void setName( const std::string& name ) { name_ = name; }

    //! Get parameter.
    /*!
     * Returns parameter.
     */
    virtual double getParameter( int i );

    //! Set parameter.
    /*!
     * Sets parameter.
     */
    virtual void setParameter( int parameterIndex, double value );

    //! Get part name.
    /*!
     * Returns part name.
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
    friend std::ostream& operator<<( std::ostream& stream, LawgsPartGeometry& lawgsPartGeometry );

protected:

    //! Part name.
    /*!
     * Part name.
     */
    std::string name_;

private:
};

} // namespace tudat

#endif // TUDAT_LAWGS_PART_GEOMETRY_H
