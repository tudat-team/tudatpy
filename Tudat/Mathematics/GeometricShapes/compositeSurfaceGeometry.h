/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110119    K. Kumar          Minor comments changes; path updated; Doxygen comments updated;
 *                                  updated function arguments to use references and unsigned ints;
 *                                  added "End of file" comment; minor changes to layout.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
#define TUDAT_COMPOSITE_SURFACE_GEOMETRY_H

#include <vector>

#include <boost/shared_ptr.hpp>

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
    CompositeSurfaceGeometry( std::vector< boost::shared_ptr< SingleSurfaceGeometry > >
                              singleSurfaceGeometryList,
                              std::vector< boost::shared_ptr< CompositeSurfaceGeometry > >
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
    boost::shared_ptr< SingleSurfaceGeometry > getSingleSurfaceGeometry( const unsigned int index )
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
    boost::shared_ptr< CompositeSurfaceGeometry > getCompositeSurfaceGeometry(
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
    friend std::ostream &operator<<( std::ostream &stream,
                                     CompositeSurfaceGeometry& compositeSurfaceGeometry );

protected:

    //! Set pointer to SingleSurfaceGeometry object.
    /*!
     * Sets a pointer to a SingleSurfaceGeometry object in singleSurfaceGeometryList_.
     * \param singleSurfaceGeometry Pointer to SingleSurfaceGeometry
     *           object which is to be stored in singleSurfaceGeometryList_.
     * \param index Index of singleSurfaceGeometryList_ at which the surface is to be set.
     */
    void setSingleSurfaceGeometry( boost::shared_ptr< SingleSurfaceGeometry >
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
    void setCompositeSurfaceGeometry( boost::shared_ptr< CompositeSurfaceGeometry>
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
    std::vector< boost::shared_ptr< SingleSurfaceGeometry > > singleSurfaceGeometryList_;

    //! Array of pointers to CompositeSurfaceGeometries.
    /*!
     *  Array of pointers to CompositeSurfaceGeometries.
     */
    std::vector< boost::shared_ptr< CompositeSurfaceGeometry > > compositeSurfaceGeometryList_;

private:
};

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_COMPOSITE_SURFACE_GEOMETRY_H
