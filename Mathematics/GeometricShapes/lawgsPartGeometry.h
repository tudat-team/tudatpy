/*! \file lawgsPartGeometry.h
 *    This file contains the definition of the LawgsPartGeometry class,
 *    Lawgs is short for Langley Wireframe Geometry Standard and is a geometry
 *    format developed by NASA (see reference) using quadrilateral surface
 *    panels.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 6 February, 2011
 *
 *    References
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard
 *          (LaWGS) format, NASA TECHNICAL MEMORANDUM 85767.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 */

#ifndef LAWGSPARTGEOMETRY_H
#define LAWGSPARTGEOMETRY_H

// Include statements.
#include "quadrilateralMeshedSurfaceGeometry.h"
#include "singleSurfaceGeometry.h"

//! Class to defined a mesh accoring to Lawgs standards.
/*!
 *  Class to define a surface mesh accoring to the Langley Wireframe Geometry
 *  standard, seet reference.
 */
class LawgsPartGeometry : public QuadrilateralMeshedSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    LawgsPartGeometry();

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~LawgsPartGeometry();

    //! Function to create a mesh surface on a single-surface geometry.
    /*!
     *  Function to create a mesh surface on a single-surface geometry. A meshed
     *  surface based on the given surface geometry is created.
     *  \param originalSurface surface from which a mesh is to be created.
     *  \param numberOfLines Number of points to be sampled from 1st
     *  independent variable.
     *  \param numberOfPoints Number of points to be sampled from 2nd
     *  independent variable.
     */
    void setMesh( SingleSurfaceGeometry* originalSurface,
                  int numberOfLines,
                  int numberOfPoints );

    //! Copy constructor.
    /*!
     *  Copy constructor to deep-copy contents of a LawgsPartGeomtry
     *  object to a new one.
     */
     LawgsPartGeometry( const LawgsPartGeometry& PartIn );

    //! Function to retrieve surface point.
    /*!
     *  This function retrieves a surface point from on a panel. Because four
     *  points define a panel in this mesh format, the collection of panels
     *  will not be necessarilly watertight, making this function non-trivial,
     *  since the four points defining the panle will not necessarily lie in
     *  the same plane.
     *  \param independentVariable1 Independent variable indicating on (or
     *  between) which lines to retrieve a point.
     *  \param independentVariable2 Independent variable indicating on (or
     *  between) which points to retrieve a point.
     *  \return point on mesh panel.
     */
    virtual VectorXd getSurfacePoint( const double& independentVariable1,
                                      const double& independentVariable2 );

    //! Function to retrieve surface derivative (currently not implemented).
    /*!
     *  Currently unavailable function to return surface derivative.
     */
    virtual VectorXd getSurfaceDerivative( const double& u ,
                                           const double& v ,
                                           const int& uDerivative ,
                                           const int& vDerivative );

    //! Function to set the name of a Lawgs part.
    /*!
     *  Function to set the name of a Lawgs part.
     */
    void setName( const std::string& name );

    //! Non-functional function to retrieve parameter.
    /*!
     *  Non functional function of base class, could be implemented in future
     *  modification.
     */
    virtual double getParameter( const int& i );

    //! Non-functional function to set parameter.
    /*!
     *  Non functional function of base class, could be implemented in future
     *  modification.
     */
    virtual void setParameter( const int &parameterIndex, const double &value );

    //! Function to de-allocate arrays in object.
    /*!
     * Function to de-allocate arrays in object.
     */
    void deAllocateArrays();

    //! Function to retrieve part name.
    /*!
     *  Function to retrieve part name.
     */
    std::string getName( );

    //! Overloaded ostream to print class information.
    /*!
     *  Overloads ostream to print class information, prints the number of
     *  lines and points, and the name of the part.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     LawgsPartGeometry& geometry );

protected:

    //! Part name.
    /*!
     *  Part name.
     */
    std::string name_;

};

#endif // LAWGSPARTGEOMETRY_H

// End of file.
