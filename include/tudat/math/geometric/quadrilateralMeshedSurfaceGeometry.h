/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      An example of a heritage code which uses such a mesh is found in:
 *          The Mark IV Supersonic-Hypersonic Arbitrary Body Program, Volume
 *          II-Program Formulation, Douglas Aircraft Company, AFFDL-TR-73-159, Volume II.
 *
 *    Notes
 *      The numberOfLines_ and numberOfPoints_ member variables denote the number of mesh points.
 *      The number of panels in the mesh will be numberOfLines_ - 1 by numberOfPoints_ - 1.
 *
 */

#ifndef TUDAT_QUADRILATERAL_MESHED_SURFACE_GEOMETRY_H
#define TUDAT_QUADRILATERAL_MESHED_SURFACE_GEOMETRY_H

#include <boost/multi_array.hpp>

#include <Eigen/Core>

#include "tudat/math/geometric/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Class for quadrilateral meshed surface geometry.
/*!
 * Base class for quadrilateral meshed surface geometry.
 */
class QuadrilateralMeshedSurfaceGeometry : public SingleSurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    QuadrilateralMeshedSurfaceGeometry( )
        : numberOfLines_( -0 ),
          numberOfPoints_( -0 ),
          reversalOperator_( 1 ),
          totalArea_( -0.0 )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~QuadrilateralMeshedSurfaceGeometry( ) { }

    //! Get surface point.
    /*!
     * Returns a surface point.
     * \param lineIndex Line on which to retrieve point.
     * \param pointIndex Point from line to retrieve.
     * \return Surface point.
     */
    Eigen::Vector3d getMeshPoint( const int lineIndex, const int pointIndex )
    {
        return meshPoints_[ lineIndex ][ pointIndex ];
    }

    //! Get panel area.
    /*!
     * Returns a panel area.
     * \param lineIndex Line of panels from which to retrieve.
     * \param pointIndex Point from line from which to retrieve.
     * \return Panel area
     */
    double getPanelArea( const int lineIndex, const int pointIndex )
    {
        return panelAreas_[ lineIndex ][ pointIndex ];
    }

    //! Get panel centroid.
    /*!
     * Returns a panel centroid.
     * \param lineIndex Line of panels from which to retrieve.
     * \param pointIndex Point from line from which to retrieve.
     * \return Panel centroid.
     */
    Eigen::Vector3d getPanelCentroid( const int lineIndex, const int pointIndex )
    {
        return panelCentroids_[ lineIndex ][ pointIndex ];
    }

    //! Get outward panel surface normal.
    /*!
     * Returns an outward panel surface normal.
     * \param lineIndex Line of panels from which to retrieve.
     * \param pointIndex Point from line from which to retrieve.
     * \return Outward panel surface normal.
     */
    Eigen::Vector3d getPanelSurfaceNormal( const int lineIndex, const int pointIndex )
    {
        return panelSurfaceNormals_[ lineIndex ][ pointIndex ];
    }

    //! Get number of lines.
    /*!
     * Returns number of lines.
     * \return Number of lines on mesh.
     */
    int getNumberOfLines( ) { return numberOfLines_; }

    //! Get number of points.
    /*!
     * Returns number of points.
     * \return Number of points on mesh.
     */
    int getNumberOfPoints( ) { return numberOfPoints_; }

    //! Get total area of the mesh.
    /*!
     * Returns the total area of the mesh.
     * \return Total mesh area.
     */
    double getTotalArea( ) { return totalArea_; }

    //! Set reversal operator.
    /*!
     * Sets if the panels in the mesh are to be inverted.
     * \param isMeshInverted Boolean to denote whether the mesh panels are to
     * be inverted.
     */
    void setReversalOperator( const bool isMeshInverted );

    //! Get boolean denoting if the mesh is inverted.
    /*!
     * Returns a boolean denoting if the mesh is inverted.
     * \return Boolean which is true if mesh is inverted, false if not.
     */
     bool getReversalOperator( );

     //! Overload ostream to print class information.
     /*!
      * Overloads ostream to print class information, prints the number of
      * lines and points, and the name of the part.
      * \param stream Stream object.
      * \param quadrilateralMeshedSurfaceGeometry Quadrilateral meshed surface
      *         geometry.
      * \return Stream object.
      */
     friend std::ostream& operator << ( std::ostream& stream,
                                      QuadrilateralMeshedSurfaceGeometry&
                                      quadrilateralMeshedSurfaceGeometry );

     boost::multi_array< Eigen::Vector3d, 2 > getMeshPoints( )
     {
         return meshPoints_;
     }

     //! Panel centroids.
     /*!
      * 2-Dimensional array containing panel centroid locations.
      */
     boost::multi_array< Eigen::Vector3d, 2 > getPanelCentroids( )
     {
         return panelCentroids_;
     }

     //! Panel centroids.
     /*!
      * 2-Dimensional array containing panel centroid locations.
      */
     boost::multi_array< Eigen::Vector3d, 2 > getPanelSurfaceNormals( )
     {
         return panelSurfaceNormals_;
     }

     void clear( )
     {
         boost::array< int, 2 > numberOfPointsPerIndependentVariables;
         numberOfPointsPerIndependentVariables[ 0 ] = 0;
         numberOfPointsPerIndependentVariables[ 1 ] = 0;

         meshPoints_.resize( numberOfPointsPerIndependentVariables );
         panelCentroids_.resize( numberOfPointsPerIndependentVariables );
         panelSurfaceNormals_.resize( numberOfPointsPerIndependentVariables );
         panelAreas_.resize( numberOfPointsPerIndependentVariables );
     }

protected:

     //! Calculate panel characteristics.
     /*!
      * Calculates the normal, centroid and area of panels in mesh.
      */
     void performPanelCalculations( );

    //! Number of lines in mesh.
    /*!
     * Number of lines in mesh.
     */
    int numberOfLines_;

    //! Number of points in mesh.
    /*!
     * Number of points in mesh.
     */
    int numberOfPoints_;

    //! Variable to denote whether mesh orientation is to be inverted.
    /*!
     * Variable to denote whether mesh orientation is to be inverted, -1 if
     * it is inverted, 1 if not.
     */
    int reversalOperator_;

    //! Mesh points.
    /*!
     * 2-Dimensional array containing mesh point locations.
     */
    boost::multi_array< Eigen::Vector3d, 2 > meshPoints_;

    //! Panel centroids.
    /*!
     * 2-Dimensional array containing panel centroid locations.
     */
    boost::multi_array< Eigen::Vector3d, 2 > panelCentroids_;

    //! Panel surface normals.
    /*!
     * 2-Dimensional array containing outward panel surface normal vectors.
     */
    boost::multi_array< Eigen::Vector3d, 2 > panelSurfaceNormals_;

    //! Panel doubles.
    /*!
     * 2-Dimensional array containing panel areas.
     */
    boost::multi_array< double, 2 > panelAreas_;

    //! Total mesh surface area/
    /*!
     * Total mesh surface area, contains the sum of all areas in panelAreas_.
     */
    double totalArea_;

private:
};

//! Typedef for shared-pointer to QuadrilateralMeshedSurfaceGeometry object.
typedef std::shared_ptr< QuadrilateralMeshedSurfaceGeometry >
QuadrilateralMeshedSurfaceGeometryPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_QUADRILATERAL_MESHED_SURFACE_GEOMETRY_H
