/*! \file compositeSurfaceGeometry.h
 *    This file contains the definition of the CompositeSurfaceGeometry class
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *      Currently, pointer arrays are used to stored SingleSurfaceGeometry and
 *      CompositeSurfaceGeometry classes. In future, it might be fruitful to
 *      consider the use of standard STL containers, such as vector and map
 *      instead.
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
 *      102511    D. Dirkx          First version of file.
 *      110119    K. Kumar          Minor comments changes; path updated;
 *                                  Doxygen comments updated; updated function
 *                                  arguments to use references and unsigned
 *                                  ints; added "End of file" comment; minor
 *                                  changes to layout.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef COMPOSITESURFACEGEOMETRY_H
#define COMPOSITESURFACEGEOMETRY_H

// Include statements.
#include "singleSurfaceGeometry.h"
#include "surfaceGeometry.h"

//! Composite surface geometry class.
/*!
 * Class for a composite surface geometry, which can be composed of any
 * number of single surface or composite surface geometries, to be set by the
 * user. This class is a container class for surface geometries.
 *
 * SingleSurfaceGeometry objects should be retreived from this class ( or
 * recursively from a member of CompositeSurfaceGeometry objects ) to perform
 * geometric operations.
 */
class CompositeSurfaceGeometry: public SurfaceGeometry
{
public:

    //! Default constructor.
    /*!
     * Default constructor, sets number of single and composite geometries to
     * zero.
     */
    CompositeSurfaceGeometry( );

    //! Default destructor.
    /*!
     * Default destructor, deallocates surface lists and resets them to NULL.
     */
    virtual ~CompositeSurfaceGeometry();

    //! Set pointer to SingleSurfaceGeometry object.
    /*!
     * Sets a pointer to a SingleSurfaceGeometry object in
     * singleSurfaceGeometryList_.
     * \param pointerToSingleSurfaceGeometry Pointer to SingleSurfaceGeometry
     *           object which is to be stored in singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ at which the surface is
     *          to be set.
     */
    void setSingleSurfaceGeometry( SingleSurfaceGeometry*
                                   pointerToSingleSurfaceGeometry,
                                   const unsigned int& index );

    //! Set pointer to a CompositeSurfaceGeometry object.
    /*!
     * Sets a pointer to a CompositeSurfaceGeometry in
     * compositeSurfaceGeometryList_.
     * \param pointerToCompositeSurfaceGeometry Pointer to
     *          CompositeSurfaceGeometry object which is to be stored in
     *          compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ at which the surface
     *          is to be set.
     */
    void setCompositeSurfaceGeometry( CompositeSurfaceGeometry*
                                      pointerToCompositeSurfaceGeometry,
                                      const unsigned int& index );

    //! Set number of single surface geometries.
    /*!
     * Sets and allocates number of SingleSurfaceGeometry objects.
     * \param numberOfSingleSurfaceGeometries Number of SingleSurfaceGeometry
     *          objects stored.
     */
    void setNumberOfSingleSurfaceGeometries( const unsigned int&
                                             numberOfSingleSurfaceGeometries );

    //! Set number of composite surface geometries.
    /*!
     * Sets and allocates number of CompositeSurfaceGeometry objects.
     * \param numberOfCompositeSurfaceGeometries Number of
     *          CompositeSurfaceGeometry objects stored.
     */
    void setNumberOfCompositeSurfaceGeometries(
            const unsigned int& numberOfCompositeSurfaceGeometries );

    //! Get pointer to stored SingleSurfaceGeometry object.
    /*!
     * Returns a pointer to a SingleSurfaceGeometry object stored in
     * singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ where the object to be
     *          retrieved is stored.
     * \return Pointer to SingleSurfaceGeometry object at index location in
     *          singleSurfaceGeometryList_.
     */
    SingleSurfaceGeometry* getSingleSurfaceGeometry( const unsigned int&
                                                     index );

    //! Get pointer to a stored CompositeSurfaceGeometry object.
    /*!
     * Returns a pointer to a CompositeSurfaceGeometry object stored in
     * compositeSurfaceGeometryList_.
     * \param index Index of compositeSurfaceGeometryList_ where the object to
     *          be retrieved is stored.
     * \return Pointer to CompositeSurfaceGeometry object at index location in
     *          compositeSurfaceGeometryList_.
     */
    CompositeSurfaceGeometry* getCompositeSurfaceGeometry( const unsigned int&
                                                           index );

    //! Get number of single surface geometries.
    /*!
     * Returns the number of single surface geometries stored in
     * singleSurfaceGeometryList_.
     * \return Number of SingleSurfaceGeometry objects stored in class.
     */
    unsigned int& getNumberOfSingleSurfaceGeometries( );

    //! Get number of composite surface geometries.
    /*!
     * Returns the number of composite surface geometries stored in
     * compositeSurfaceGeometryList_.
     * \return Number of CompositeSurfaceGeometry objects stored in class.
     */
    unsigned int& getNumberOfCompositeSurfaceGeometries( );

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
                                     CompositeSurfaceGeometry&
                                     compositeSurfaceGeometry );

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
    SingleSurfaceGeometry** singleSurfaceGeometryList_;

    //! Array of pointers to CompositeSurfaceGeometries.
    /*!
     *  Array of pointers to CompositeSurfaceGeometries.
     */
    CompositeSurfaceGeometry** compositeSurfaceGeometryList_;

private:
};

#endif // COMPOSITESURFACEGEOMETRY_H

// End of file.
