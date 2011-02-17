/*! \file quadrilateralMeshedSurfaceGeometry.h
 *    This file contains the definition of the quadrilateral meshed surface
 *    geometry class.
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
 *      An example of a heritage code which uses such a mesh is found in:
 *          The Mark IV Supersonic-Hypersonic Arbitrary Body Program, Volume
 *          II-Program Formulation, Douglas Aircraft Company, AFFDL-TR-73-159,
 *          Volume II.
 *
 *    Notes
 *      This class uses pointers to pointers to denote two-dimensional arrays
 *      for heritage reasons, instead of the more modern C++ vector or map
 *      types.
 *
 *      The numberOfLines_ and numberOfPoints_ member variables denote the
 *      number of mesh points. The number of panels in the mesh will be
 *      numberOfLines_ - 1 by numberOfPoints_ -1.
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

#ifndef QUADRILATERALMESHEDSURFACEGEOMETRY_H
#define QUADRILATERALMESHEDSURFACEGEOMETRY_H

// Include statements.
#include "singleSurfaceGeometry.h"

//! Class for quadrilateral meshed surface geometry.
/*!
 *  Base class for quadrilateral meshed surface geometry.
 */
class QuadrilateralMeshedSurfaceGeometry : public SingleSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    QuadrilateralMeshedSurfaceGeometry( );

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~QuadrilateralMeshedSurfaceGeometry( );

    //! Function to calculate panel characteristics.
    /*!
     *  This function calculates the normal, centroid and area of panels in mesh.
     */
    void performPanelCalculations( );

    //! Function to retrieve a surface point.
    /*!
     *  Function to retrieve a surface point.
     *  \param lineIndex Line on which to retrieve point.
     *  \param pointIndex Point from line to retrieve.
     *  \return Surface point.
     */
    Vector3d getMeshPoint( const int& lineIndex, const int& pointIndex );

    //! Function to retrieve a panel area.
    /*!
     *  Function to retrieve a panel area.
     *  \param lineIndex Line of panels from which to retrieve.
     *  \param pointIndex Point from line from which to retrieve.
     *  \return Panel area
     */
    double getPanelArea( const int& lineIndex, const int& pointIndex );

    //! Function to retrieve a panel centroid.
    /*!
     *  Function to retrieve a panel centroid.
     *  \param lineIndex Line of panels from which to retrieve.
     *  \param pointIndex Point from line from which to retrieve.
     *  \return Panel centroid.
     */
    Vector3d getPanelCentroid( const int& lineIndex, const int& pointIndex );

    //! Function to retrieve an outward panel surface normal
    /*!
     *  Function to retrieve an outward panel surface normal.
     *  \param lineIndex Line of panels from which to retrieve.
     *  \param pointIndex Point from line from which to retrieve.
     *  \return Outward panel surface normal.
     */
    Vector3d getPanelSurfaceNormal( const int& lineIndex, const int& pointIndex );

    //! Function to retrieve number of lines.
    /*!
     *  Function to retrieve number of lines.
     *  \return Number of lines on mesh.
     */
    int getNumberOfLines( );

    //! Function to retrieve number of points.
    /*!
     *  Function to retrieve number of points.
     *  \return Number of points on mesh.
     */
    int getNumberOfPoints( );

    //! Returns the total area of the mesh.
    /*!
     * Returns the total area of the mesh.
     * \return Total mesh area.
     */
    double getTotalArea( );

    //! Sets reversal operator.
    /*!
     *  Function set if the panels in the mesh are to be inverted.
     *  \param isMeshInverted Boolean to denote whether the mesh panels are to
     *  be inverted.
     */
    void setReversalOperator( bool isMeshInverted );

    //! Returns a boolean denoting if the mesh is inverted.
    /*!
     *  Returns a boolean denoting if the mesh is inverted.
     *  \return Boolean which is true if mesh is inverted, false if not.
     */
     bool getReversalOperator( );

     //! Overload ostream to print class information.
     /*!
      *  Overloads ostream to print class information, prints the number of
      *  lines and points, and the name of the part.
      */
     friend std::ostream& operator<<( std::ostream& stream,
                                      QuadrilateralMeshedSurfaceGeometry&
                                      geometry );

protected:

    //! Number of lines in mesh.
    /*!
     *  Number of lines in mesh.
     */
    int numberOfLines_;

    //! Number of points in mesh.
    /*!
     *  Number of points in mesh.
     */
    int numberOfPoints_;

    //! Variable to denote whether mesh orientation is to be inverted.
    /*!
     *  Variable to denote whether mesh orientation is to be inverted, -1 if
     *  it is inverted, 1 if not.
     */
    int reversalOperator_;

    //! Mesh points.
    /*!
     *  2-Dimensional array containing mesh point locations.
     */
    Vector3d** meshPoints_;

    //! Panel centroids.
    /*!
     *  2-Dimensional array containing panel centroid locations.
     */
    Vector3d** panelCentroids_;

    //! Panel surface normals.
    /*!
     *  2-Dimensional array containing outward panel surface normal vectors.
     */
    Vector3d** panelSurfaceNormals_;

    //! Panel doubles.
    /*!
     *  2-Dimensional array containing panel areas.
     */
    double** panelAreas_;

    //! Total mesh surface area/
    /*!
     *  Total mesh surface area, contains the sum of all areas in panelAreas_.
     */
    double totalArea_;
};

#endif // QUADRILATERALMESHEDSURFACEGEOMETRY_H

// End of file.
