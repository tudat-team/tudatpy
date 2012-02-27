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
 *      110119    K. Kumar          First creation of code.
 *
 *    References
 *
 */

#ifndef TUDAT_ROOT_FINDER_BASE_H
#define TUDAT_ROOT_FINDER_BASE_H

namespace tudat
{

//! An abstract base class for RootFinder.
/*!
 * An abstract base class for the RootFinder class.
 */
class RootFinderBase
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~RootFinderBase( ) { }

    //! Compute mathematical function value.
    /*!
     * Computes the value of the mathematical function used for the root-finder
     * algorithm being used.
     * \param inputValue Input value.
     */
    virtual double computeFunction( double& inputValue ) = 0;

protected:

private:
};

} // namespace tudat

#endif // TUDAT_ROOT_FINDER_BASE_H
