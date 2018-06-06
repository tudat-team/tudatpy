/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
#define TUDAT_COMPOSITE_SURFACE_GEOMETRY_H

#include <vector>

#include <memory>

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Composite surface geometry class.
/*!
 * Class for a composite surface geometry, which can be composed of any
 * number of single surface or composite surface geometries, to be set by the
 * user. This class is a container class for surface geometries.
 *
 * SingleSurfaceGeometry objects should be retrieved from this class ( or
 * recursively from a member of CompositeSurfaceGeometry objects ) to perform
 * geometric operations.
 */
class CompositeSurfaceGeometry: public SurfaceGeometry
{
public:

    //! Typedef for shared-pointer to CompositeSurfaceGeometry object.
    typedef std::shared_ptr< CompositeSurfaceGeometry > CompositeSurfaceGeometryPointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CompositeSurfaceGeometry( ) { }

    //! Constructor, sets constituent surface geometries of objects.
    /*!
     *  Constructor, sets constituent surface geometries of objects.
     *  \param singleSurfaceGeometryList vector of pointers to SingleSurfaceGeometries
     *  \param compositeSurfaceGeometryList vector of pointers to CompositeSurfaceGeometries
     */
    CompositeSurfaceGeometry( std::vector< std::shared_ptr< SingleSurfaceGeometry > >
                              singleSurfaceGeometryList,
                              std::vector< std::shared_ptr< CompositeSurfaceGeometry > >
                              compositeSurfaceGeometryList );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~CompositeSurfaceGeometry( ) { }

    //! Get pointer to stored SingleSurfaceGeometry object.
    /*!
     * Returns a pointer to a SingleSurfaceGeometry object stored in
     * singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ where the object to be
     *          retrieved is stored.
     * \return Pointer to SingleSurfaceGeometry object at index location in
     *          singleSurfaceGeometryList_.
     */
    std::shared_ptr< SingleSurfaceGeometry > getSingleSurfaceGeometry( const unsigned int index )
    {
        return singleSurfaceGeometryList_[ index ];
    }

    //! Get pointer to a stored CompositeSurfaceGeometry object.
    /*!
     * Returns a pointer to a CompositeSurfaceGeometry object stored in
     * compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ where the object to
     *          be retrieved is stored.
     * \return Pointer to CompositeSurfaceGeometry object at index location in
     *          compositeSurfaceGeometryList_.
     */
    std::shared_ptr< CompositeSurfaceGeometry > getCompositeSurfaceGeometry(
            const unsigned int index )
    {
        return compositeSurfaceGeometryList_[ index ];
    }

    //! Get number of single surface geometries.
    /*!
     * Returns the number of single surface geometries stored in
     * singleSurfaceGeometryList_.
     * \return Number of SingleSurfaceGeometry objects stored in class.
     */
    unsigned int getNumberOfSingleSurfaceGeometries( )
    {
        return numberOfSingleSurfaceGeometries_;
    }

    //! Get number of composite surface geometries.
    /*!
     * Returns the number of composite surface geometries stored in
     * compositeSurfaceGeometryList_.
     * \return Number of CompositeSurfaceGeometry objects stored in class.
     */
    unsigned int getNumberOfCompositeSurfaceGeometries( )
    {
        return numberOfCompositeSurfaceGeometries_;
    }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * and number of constituent SingleSurfaceGeometry and
     * compositeSurfaceGeometry objects.
     * \param stream Stream object.
     * \param compositeSurfaceGeometry compositeSurfaceGeometry object.
     * \return Stream object.
     */
    friend std::ostream &operator << ( std::ostream &stream,
                                     CompositeSurfaceGeometry& compositeSurfaceGeometry );

protected:

    //! Set pointer to SingleSurfaceGeometry object.
    /*!
     * Sets a pointer to a SingleSurfaceGeometry object in singleSurfaceGeometryList_.
     * \param singleSurfaceGeometry Pointer to SingleSurfaceGeometry
     *           object which is to be stored in singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ at which the surface is to be set.
     */
    void setSingleSurfaceGeometry( std::shared_ptr< SingleSurfaceGeometry >
                                   singleSurfaceGeometry, const unsigned int index )
    {
        singleSurfaceGeometryList_[ index ] = singleSurfaceGeometry;
    }

    //! Set pointer to a CompositeSurfaceGeometry object.
    /*!
     * Sets a pointer to a CompositeSurfaceGeometry in
     * compositeSurfaceGeometryList_.
     * \param compositeSurfaceGeometry Pointer to
     *          CompositeSurfaceGeometry object which is to be stored in
     *          compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ at which the surface is to be set.
     */
    void setCompositeSurfaceGeometry( std::shared_ptr< CompositeSurfaceGeometry>
                                      compositeSurfaceGeometry, const unsigned int index )
    {
        compositeSurfaceGeometryList_[ index ] = compositeSurfaceGeometry;
    }

    //! Set number of single surface geometries.
    /*!
     * Sets and allocates number of SingleSurfaceGeometry objects.
     * \param numberOfSingleSurfaceGeometries Number of SingleSurfaceGeometry objects stored.
     */
    void setNumberOfSingleSurfaceGeometries( const unsigned int numberOfSingleSurfaceGeometries )
    {
        numberOfSingleSurfaceGeometries_ = numberOfSingleSurfaceGeometries;
        singleSurfaceGeometryList_.resize( numberOfSingleSurfaceGeometries_ );
    }

    //! Set number of composite surface geometries.
    /*!
     * Sets and allocates number of CompositeSurfaceGeometry objects.
     * \param numberOfCompositeSurfaceGeometries Number of
     *          CompositeSurfaceGeometry objects stored.
     */
    void setNumberOfCompositeSurfaceGeometries( const unsigned int
                                                numberOfCompositeSurfaceGeometries )
    {
        numberOfCompositeSurfaceGeometries_ = numberOfCompositeSurfaceGeometries;
        compositeSurfaceGeometryList_.resize( numberOfCompositeSurfaceGeometries_ );
    }

    //! Size of singleSurfaceGeometryList_.
    /*!
     *  Size of singleSurfaceGeometryList_.
     */
    unsigned int numberOfSingleSurfaceGeometries_;

    //! Size of compositeSurfaceGeometryList_.
    /*!
     *  Size of compositeSurfaceGeometryList_.
     */
    unsigned int numberOfCompositeSurfaceGeometries_;

    //! Array of pointers to SingleSurfaceGeometries.
    /*!
     *  Array of pointers to SingleSurfaceGeometries.
     */
    std::vector< SingleSurfaceGeometryPointer > singleSurfaceGeometryList_;

    //! Array of pointers to CompositeSurfaceGeometries.
    /*!
     *  Array of pointers to CompositeSurfaceGeometries.
     */
    std::vector< CompositeSurfaceGeometryPointer > compositeSurfaceGeometryList_;

private:
};

//! Typedef for shared-pointer to CompositeSurfaceGeometry object.
typedef CompositeSurfaceGeometry::CompositeSurfaceGeometryPointer CompositeSurfaceGeometryPointer;

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
