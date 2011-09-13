/*! \file randomNumberGenerator.h
 *    This header file contains a base class for all random number generators
 *    in Tudat.
 *
 *    Path              : /Mathematics/RandomNumberGenerators/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 7 October, 2010
 *    Last modified     : 17 May, 2011
 *
 *    References

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
 *      101007    K. Kumar          First creation of code.
 *      101015    K. Kumar          Added RandomNumbers structure.
 *      101020    K. Kumar          Consolidated code into single class and
 *                                  completed code comments. Object can now
 *                                  only be created by explicit seeding.
 *      101025    K. Kumar          Updated code based on D. Dirkx's codecheck
 *                                  comments.
 *      110117    K. Kumar          Corrected path.
 *      110121    K. Kumar          Added comment to "Notes".
 *      110517    K. Kumar          Converted file into generic base class.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

//! Random number generator base class.
/*!
 * Base class for random number generators in Tudat.
 */
class RandomNumberGenerator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RandomNumberGenerator( ) {}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~RandomNumberGenerator( )    {}

protected:

private:
};

#endif // RANDOMNUMBERGENERATOR_H

// End of file.
