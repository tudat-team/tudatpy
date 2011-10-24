/*! \file rootFinderBase.h
 *    This header file contains an abstract base class for the RootFinder class included in Tudat.
 *
 *    Path              : /Mathematics/RootFindingMethods/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 19 January, 2011
 *    Last modified     : 19 January, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 */

#ifndef ROOTFINDERBASE_H
#define ROOTFINDERBASE_H

//! An abstract base class for RootFinder.
/*!
 * This is an abstract base class for the RootFinder class.
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

#endif // ROOTFINDERBASE_H

// End of file.
