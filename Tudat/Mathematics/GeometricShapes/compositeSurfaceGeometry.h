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
 *      102511    D. Dirkx          First version of file.
 *      110119    K. Kumar          Minor comments changes; path updated; Doxygen comments updated;
 *                                  updated function arguments to use references and unsigned ints;
 *                                  added "End of file" comment; minor changes to layout.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Currently, pointer arrays are used to stored SingleSurfaceGeometry and
// CompositeSurfaceGeometry classes. In future, it might be fruitful to
// consider the use of standard STL containers, such as vector and map
// instead.
// 

#ifndef TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
#define TUDAT_COMPOSITE_SURFACE_GEOMETRY_H

#include <vector>

#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/surfaceGeometry.h"

namespace tudat
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

    //! Default constructor.
    /*!
     * Default constructor, sets number of single and composite geometries to zero.
     */
    CompositeSurfaceGeometry( ) : numberOfSingleSurfaceGeometries_( 0 ),
        numberOfCompositeSurfaceGeometries_( 0 ), singleSurfaceGeometryList_( NULL ),
        compositeSurfaceGeometryList_( NULL ) { }

    //! Default destructor.
    /*!
     * Default destructor, deallocates surface lists and resets them to NULL.
     */
    virtual ~CompositeSurfaceGeometry( ) { }

    //! Set pointer to SingleSurfaceGeometry object.
    /*!
     * Sets a pointer to a SingleSurfaceGeometry object in singleSurfaceGeometryList_.
     * \param pointerToSingleSurfaceGeometry Pointer to SingleSurfaceGeometry
     *           object which is to be stored in singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ at which the surface is to be set.
     */
    void setSingleSurfaceGeometry( SingleSurfaceGeometry* pointerToSingleSurfaceGeometry,
                                   const unsigned int& index )
    { singleSurfaceGeometryList_[ index ] = pointerToSingleSurfaceGeometry; }

    //! Set pointer to a CompositeSurfaceGeometry object.
    /*!
     * Sets a pointer to a CompositeSurfaceGeometry in
     * compositeSurfaceGeometryList_.
     * \param pointerToCompositeSurfaceGeometry Pointer to
     *          CompositeSurfaceGeometry object which is to be stored in
     *          compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ at which the surface is to be set.
     */
    void setCompositeSurfaceGeometry( CompositeSurfaceGeometry* pointerToCompositeSurfaceGeometry,
                                      const unsigned int& index )
    { compositeSurfaceGeometryList_[ index ] = pointerToCompositeSurfaceGeometry; }

    //! Set number of single surface geometries.
    /*!
     * Sets and allocates number of SingleSurfaceGeometry objects.
     * \param numberOfSingleSurfaceGeometries Number of SingleSurfaceGeometry objects stored.
     */
    void setNumberOfSingleSurfaceGeometries( const unsigned int& numberOfSingleSurfaceGeometries )
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
    void setNumberOfCompositeSurfaceGeometries( const unsigned int&
                                                numberOfCompositeSurfaceGeometries )
    {
        numberOfCompositeSurfaceGeometries_ = numberOfCompositeSurfaceGeometries;
        compositeSurfaceGeometryList_.resize( numberOfCompositeSurfaceGeometries_ );
    }

    //! Get pointer to stored SingleSurfaceGeometry object.
    /*!
     * Returns a pointer to a SingleSurfaceGeometry object stored in
     * singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ where the object to be
     *          retrieved is stored.
     * \return Pointer to SingleSurfaceGeometry object at index location in
     *          singleSurfaceGeometryList_.
     */
    SingleSurfaceGeometry* getSingleSurfaceGeometry( const unsigned int& index ) 
    { return singleSurfaceGeometryList_[ index ]; }

    //! Get pointer to a stored CompositeSurfaceGeometry object.
    /*!
     * Returns a pointer to a CompositeSurfaceGeometry object stored in
     * compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ where the object to
     *          be retrieved is stored.
     * \return Pointer to CompositeSurfaceGeometry object at index location in
     *          compositeSurfaceGeometryList_.
     */
    CompositeSurfaceGeometry* getCompositeSurfaceGeometry( const unsigned int& index ) 
    { return compositeSurfaceGeometryList_[index]; }

    //! Get number of single surface geometries.
    /*!
     * Returns the number of single surface geometries stored in
     * singleSurfaceGeometryList_.
     * \return Number of SingleSurfaceGeometry objects stored in class.
     */
    unsigned int& getNumberOfSingleSurfaceGeometries( )
    { return numberOfSingleSurfaceGeometries_; }

    //! Get number of composite surface geometries.
    /*!
     * Returns the number of composite surface geometries stored in
     * compositeSurfaceGeometryList_.
     * \return Number of CompositeSurfaceGeometry objects stored in class.
     */
    unsigned int& getNumberOfCompositeSurfaceGeometries( )
    { return numberOfCompositeSurfaceGeometries_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the class type,
     * and number of constituent SingleSurfaceGeometry and
     * compositeSurfaceGeometry objects.
     * \param stream Stream object.
     * \param compositeSurfaceGeometry compositeSurfaceGeometry object.
     * \return Stream object.
     */
    friend std::ostream &operator<<( std::ostream &stream,
                                     CompositeSurfaceGeometry& compositeSurfaceGeometry );

protected:

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
    std::vector< SingleSurfaceGeometry* > singleSurfaceGeometryList_;

    //! Array of pointers to CompositeSurfaceGeometries.
    /*!
     *  Array of pointers to CompositeSurfaceGeometries.
     */
    std::vector< CompositeSurfaceGeometry* > compositeSurfaceGeometryList_;

private:
};

} // namespace tudat

#endif // TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
